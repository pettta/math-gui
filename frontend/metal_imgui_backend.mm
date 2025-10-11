#include "metal_imgui_backend.h"

#import <QuartzCore/QuartzCore.h>
#include <SDL.h>

#include "imgui_impl_metal.h"
#include "imgui_impl_sdl2.h"

MetalImguiBackend::MetalImguiBackend(CAMetalLayer* layer, SDL_Window* window)
    : layer_(layer), window_(window)
{
}

void MetalImguiBackend::setFrameResources(void* commandBuffer,
                                          MTLRenderPassDescriptor* renderPassDescriptor,
                                          void* renderEncoder)
{
#ifdef __OBJC__
    commandBuffer_ = (__bridge id<MTLCommandBuffer>)commandBuffer;
    renderEncoder_ = (__bridge id<MTLRenderCommandEncoder>)renderEncoder;
#else
    commandBuffer_ = commandBuffer;
    renderEncoder_ = renderEncoder;
#endif
    renderPassDescriptor_ = renderPassDescriptor;
}

void MetalImguiBackend::initializeBackend()
{
#ifdef __OBJC__
    ImGui_ImplMetal_Init(layer_.device);
    ImGui_ImplSDL2_InitForMetal(window_);
#endif
}

void MetalImguiBackend::newFrame()
{
#ifdef __OBJC__
    ImGui_ImplMetal_NewFrame(renderPassDescriptor_);
    ImGui_ImplSDL2_NewFrame();
#endif
}

void MetalImguiBackend::renderDrawData(ImDrawData* drawData)
{
#ifdef __OBJC__
    if (!commandBuffer_ || !renderEncoder_)
    {
        return;
    }

    ImGui_ImplMetal_RenderDrawData(drawData, commandBuffer_, renderEncoder_);
#endif
}

void MetalImguiBackend::shutdownBackend()
{
#ifdef __OBJC__
    ImGui_ImplMetal_Shutdown();
    ImGui_ImplSDL2_Shutdown();
#endif
}
