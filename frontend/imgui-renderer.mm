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
    commandBuffer_ = (__bridge id<MTLCommandBuffer>)commandBuffer;
    renderPassDescriptor_ = renderPassDescriptor;
    renderEncoder_ = (__bridge id<MTLRenderCommandEncoder>)renderEncoder;
}

void MetalImguiBackend::initializeBackend()
{
    ImGui_ImplMetal_Init(layer_.device);
    ImGui_ImplSDL2_InitForMetal(window_);
}

void MetalImguiBackend::newFrame()
{
    ImGui_ImplMetal_NewFrame(renderPassDescriptor_);
    ImGui_ImplSDL2_NewFrame();
}

void MetalImguiBackend::renderDrawData(ImDrawData* drawData)
{
    if (!commandBuffer_ || !renderEncoder_)
    {
        return;
    }

    ImGui_ImplMetal_RenderDrawData(drawData, commandBuffer_, renderEncoder_);
}

void MetalImguiBackend::shutdownBackend()
{
    ImGui_ImplMetal_Shutdown();
    ImGui_ImplSDL2_Shutdown();
}
