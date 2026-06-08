// Dear ImGui: standalone example application for SDL2 + Metal
// (SDL is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)

// Learn about Dear ImGui:
// - FAQ                  https://dearimgui.com/faq
// - Getting Started      https://dearimgui.com/getting-started
// - Documentation        https://dearimgui.com/docs (same as your local docs/ folder).
// - Introduction, links and more at the top of imgui.cpp

#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "../frontend/imgui-renderer.h"
#include "../frontend/metal_imgui_backend.h"
#include "../frontend/gif_recorder.h"
#include <stdio.h>
#include <SDL.h>

#import <Metal/Metal.h>
#import <QuartzCore/QuartzCore.h>

int main(int, char**)
{

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details. If you like the default font but want it to scale better, consider using the 'ProggyVector' from the same author!
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //style.FontSizeBase = 20.0f;
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf");
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf");
    //IM_ASSERT(font != nullptr);

    // Setup SDL
    // (Some versions of SDL before <2.0.10 appears to have performance/stalling issues on a minority of Windows systems,
    // depending on whether SDL_INIT_GAMECONTROLLER is enabled or disabled.. updating to latest version of SDL is recommended!)
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
    {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    // Inform SDL that we will be using metal for rendering. Without this hint initialization of metal renderer may fail.
    SDL_SetHint(SDL_HINT_RENDER_DRIVER, "metal");

    // Enable native IME.
    SDL_SetHint(SDL_HINT_IME_SHOW_UI, "1");

    SDL_Window* window = SDL_CreateWindow("Dear ImGui SDL+Metal example", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 800, SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    if (window == nullptr)
    {
        printf("Error creating window: %s\n", SDL_GetError());
        return -2;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (renderer == nullptr)
    {
        printf("Error creating renderer: %s\n", SDL_GetError());
        return -3;
    }

    // Setup Platform/Renderer backends
    CAMetalLayer* layer = (__bridge CAMetalLayer*)SDL_RenderGetMetalLayer(renderer);
    layer.pixelFormat = MTLPixelFormatBGRA8Unorm;
    layer.framebufferOnly = NO; // allow reading back the drawable for GIF capture
    MetalImguiBackend imguiBackend(layer, window);
    ImGuiRenderer imguiRenderer(imguiBackend);
    imguiRenderer.initialize();
    imguiRenderer.io();

    id<MTLCommandQueue> commandQueue = [layer.device newCommandQueue];
    MTLRenderPassDescriptor* renderPassDescriptor = [MTLRenderPassDescriptor new];

    // GIF recording state. The recorder reacts to transitions of state.gif_recording
    // (toggled by the header button / Cmd+Shift+G in the renderer). While recording we
    // blit the drawable into a reusable staging buffer, throttled to kFps.
    GifRecorder recorder;
    bool wasRecording = false;
    const int kFps = 15;
    uint32_t lastCaptureTicks = 0;
    id<MTLBuffer> stagingBuffer = nullptr;
    NSUInteger stagingCapacity = 0;

    // Main loop
    bool done = false;
    while (!done)
    {
        @autoreleasepool
        {
            // Poll and handle events (inputs, window resize, etc.)
            // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
            // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
            // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
            // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
            SDL_Event event;
            while (SDL_PollEvent(&event))
            {
                ImGui_ImplSDL2_ProcessEvent(&event);
                if (event.type == SDL_QUIT)
                    done = true;
                if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window))
                    done = true;
            }

            int width, height;
            SDL_GetRendererOutputSize(renderer, &width, &height);
            layer.drawableSize = CGSizeMake(width, height);
            id<CAMetalDrawable> drawable = [layer nextDrawable];

            ImGuiRenderer::FrameState& state = imguiRenderer.frameState();

            id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
            renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(state.clear_color[0] * state.clear_color[3], state.clear_color[1] * state.clear_color[3], state.clear_color[2] * state.clear_color[3], state.clear_color[3]);
            renderPassDescriptor.colorAttachments[0].texture = drawable.texture;
            renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
            renderPassDescriptor.colorAttachments[0].storeAction = MTLStoreActionStore;
            id <MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
            [renderEncoder pushDebugGroup:@"ImGui demo"];

            imguiBackend.setFrameResources((__bridge void*)commandBuffer,
                                           renderPassDescriptor,
                                           (__bridge void*)renderEncoder);

            // Render ImGui frame (using default UI if no custom draw function is provided)
            imguiRenderer.renderFrame(nullptr);

            [renderEncoder popDebugGroup];
            [renderEncoder endEncoding];

            // --- GIF recording: react to start/stop, then capture this frame ---
            bool wantRecording = state.gif_recording;
            if (wantRecording && !wasRecording)
            {
                recorder.start(width, height, kFps);
                state.gif_status = recorder.status();
                if (!recorder.isRecording())
                {
                    // start failed (e.g. ffmpeg missing) — don't leave the button stuck red
                    state.gif_recording = false;
                }
                lastCaptureTicks = 0;
            }
            else if (!wantRecording && wasRecording)
            {
                recorder.stop();
                state.gif_status = recorder.status();
            }
            wasRecording = state.gif_recording;

            bool capturedThisFrame = false;
            if (recorder.isRecording())
            {
                uint32_t now = SDL_GetTicks();
                if (lastCaptureTicks == 0 || (now - lastCaptureTicks) >= (uint32_t)(1000 / kFps))
                {
                    lastCaptureTicks = now;
                    NSUInteger needed = (NSUInteger)width * (NSUInteger)height * 4u;
                    if (!stagingBuffer || stagingCapacity < needed)
                    {
                        stagingBuffer = [layer.device newBufferWithLength:needed
                                                                 options:MTLResourceStorageModeShared];
                        stagingCapacity = needed;
                    }
                    id<MTLBlitCommandEncoder> blit = [commandBuffer blitCommandEncoder];
                    [blit copyFromTexture:drawable.texture
                              sourceSlice:0
                              sourceLevel:0
                             sourceOrigin:MTLOriginMake(0, 0, 0)
                               sourceSize:MTLSizeMake(width, height, 1)
                                 toBuffer:stagingBuffer
                        destinationOffset:0
                   destinationBytesPerRow:(NSUInteger)width * 4u
                 destinationBytesPerImage:needed];
                    [blit endEncoding];
                    capturedThisFrame = true;
                }
            }

            [commandBuffer presentDrawable:drawable];
            [commandBuffer commit];

            if (capturedThisFrame)
            {
                [commandBuffer waitUntilCompleted];
                recorder.addFrame((const uint8_t*)stagingBuffer.contents, width, height);
            }
        }
    }

    // Cleanup
    if (recorder.isRecording())
    {
        recorder.stop();
    }
    imguiRenderer.shutdown();

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
