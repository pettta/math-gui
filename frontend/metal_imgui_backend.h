#pragma once

#include "imgui-renderer.h"

#ifdef __OBJC__
@class CAMetalLayer;
@class MTLRenderPassDescriptor;
@protocol CAMetalDrawable;
@protocol MTLCommandBuffer;
@protocol MTLRenderCommandEncoder;
#else
typedef struct objc_object CAMetalLayer;
typedef struct objc_object MTLRenderPassDescriptor;
typedef struct objc_object MTLCommandBuffer;
typedef struct objc_object MTLRenderCommandEncoder;
#endif

struct SDL_Window;

class MetalImguiBackend final : public ImGuiBackend {
public:
    MetalImguiBackend(CAMetalLayer* layer, SDL_Window* window);

    void setFrameResources(void* commandBuffer,
                           MTLRenderPassDescriptor* renderPassDescriptor,
                           void* renderEncoder);

    void initializeBackend() override;
    void newFrame() override;
    void renderDrawData(ImDrawData* drawData) override;
    void shutdownBackend() override;

private:
    CAMetalLayer* layer_;
    SDL_Window* window_;

#ifdef __OBJC__
    id<MTLCommandBuffer> commandBuffer_ = nullptr;
    id<MTLRenderCommandEncoder> renderEncoder_ = nullptr;
#else
    void* commandBuffer_ = nullptr;
    void* renderEncoder_ = nullptr;
#endif
    MTLRenderPassDescriptor* renderPassDescriptor_ = nullptr;
};
