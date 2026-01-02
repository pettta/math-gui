#pragma once

#include "imgui-renderer.h"

#include <functional>

#include <vulkan/vulkan.h>

struct SDL_Window;

// Backend that bridges Dear ImGui with SDL2 + Vulkan using dynamic rendering.
class VulkanImguiBackend final : public ImGuiBackend {
public:
    VulkanImguiBackend(SDL_Window* window,
                       VkInstance instance,
                       VkPhysicalDevice physicalDevice,
                       VkDevice device,
                       uint32_t graphicsQueueFamily,
                       VkQueue graphicsQueue,
                       VkFormat swapchainImageFormat,
                       uint32_t minImageCount,
                       uint32_t imageCount);

    void setFrameResources(VkCommandBuffer commandBuffer,
                           VkImageView targetImageView,
                           VkExtent2D extent);

    void updateSwapchainInfo(VkFormat swapchainImageFormat,
                             uint32_t imageCount,
                             uint32_t minImageCount);

    void initializeBackend() override;
    void newFrame() override;
    void renderDrawData(ImDrawData* drawData) override;
    void shutdownBackend() override;

private:
    SDL_Window* window_;
    VkInstance instance_;
    VkPhysicalDevice physicalDevice_;
    VkDevice device_;
    uint32_t graphicsQueueFamily_;
    VkQueue graphicsQueue_;
    VkFormat swapchainImageFormat_;
    uint32_t imageCount_;
    uint32_t minImageCount_;
    VkDescriptorPool descriptorPool_{VK_NULL_HANDLE};

    VkCommandBuffer commandBuffer_{VK_NULL_HANDLE};
    VkImageView targetImageView_{VK_NULL_HANDLE};
    VkExtent2D renderExtent_{};
    bool initialized_{false};
};
