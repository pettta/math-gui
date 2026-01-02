#include "vulkan_imgui_backend.h"

#include <SDL.h>

#include "imgui_impl_sdl2.h"
#include "imgui_impl_vulkan.h"

#include <vk_initializers.h>
#include <vk_types.h>

#include <stdexcept>
#include <utility>

VulkanImguiBackend::VulkanImguiBackend(SDL_Window* window,
                                       VkInstance instance,
                                       VkPhysicalDevice physicalDevice,
                                       VkDevice device,
                                       uint32_t graphicsQueueFamily,
                                       VkQueue graphicsQueue,
                                       VkFormat swapchainImageFormat,
                                       uint32_t minImageCount,
                                       uint32_t imageCount)
    : window_(window)
    , instance_(instance)
    , physicalDevice_(physicalDevice)
    , device_(device)
    , graphicsQueueFamily_(graphicsQueueFamily)
    , graphicsQueue_(graphicsQueue)
    , swapchainImageFormat_(swapchainImageFormat)
    , imageCount_(imageCount)
    , minImageCount_(minImageCount)
{
}

void VulkanImguiBackend::setFrameResources(VkCommandBuffer commandBuffer,
                                           VkImageView targetImageView,
                                           VkExtent2D extent)
{
    commandBuffer_ = commandBuffer;
    targetImageView_ = targetImageView;
    renderExtent_ = extent;
}

void VulkanImguiBackend::updateSwapchainInfo(VkFormat swapchainImageFormat,
                                             uint32_t imageCount,
                                             uint32_t minImageCount)
{
    swapchainImageFormat_ = swapchainImageFormat;
    imageCount_ = imageCount;
    minImageCount_ = minImageCount;
    if (initialized_)
    {
        ImGui_ImplVulkan_SetMinImageCount(minImageCount_);
    }
}

void VulkanImguiBackend::initializeBackend()
{
    if (initialized_)
    {
        return;
    }

    VkDescriptorPoolSize poolSizes[] = {
        {VK_DESCRIPTOR_TYPE_SAMPLER, 1000},
        {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1000},
        {VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_TEXEL_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_TEXEL_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 1000},
        {VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 1000},
    };

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    poolInfo.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    poolInfo.maxSets = 1000;
    poolInfo.poolSizeCount = static_cast<uint32_t>(std::size(poolSizes));
    poolInfo.pPoolSizes = poolSizes;

    VK_CHECK(vkCreateDescriptorPool(device_, &poolInfo, nullptr, &descriptorPool_));

    ImGui_ImplSDL2_InitForVulkan(window_);

    ImGui_ImplVulkan_InitInfo initInfo{};
    initInfo.Instance = instance_;
    initInfo.PhysicalDevice = physicalDevice_;
    initInfo.Device = device_;
    initInfo.QueueFamily = graphicsQueueFamily_;
    initInfo.Queue = graphicsQueue_;
    initInfo.PipelineCache = VK_NULL_HANDLE;
    initInfo.DescriptorPool = descriptorPool_;
    initInfo.Subpass = 0;
    initInfo.MinImageCount = minImageCount_;
    initInfo.ImageCount = imageCount_;
    initInfo.MSAASamples = VK_SAMPLE_COUNT_1_BIT;
    initInfo.UseDynamicRendering = true;
    initInfo.Allocator = nullptr;

    VkPipelineRenderingCreateInfo pipelineInfo{};
    pipelineInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO;
    pipelineInfo.colorAttachmentCount = 1;
    pipelineInfo.pColorAttachmentFormats = &swapchainImageFormat_;
    initInfo.PipelineRenderingCreateInfo = pipelineInfo;

    if (!ImGui_ImplVulkan_Init(&initInfo))
    {
        throw std::runtime_error("Failed to initialize ImGui Vulkan backend.");
    }

    initialized_ = true;
}

void VulkanImguiBackend::newFrame()
{
    if (!initialized_)
    {
        return;
    }

    ImGui_ImplVulkan_NewFrame();
    ImGui_ImplSDL2_NewFrame();
}

void VulkanImguiBackend::renderDrawData(ImDrawData* drawData)
{
    if (!initialized_ || !drawData || drawData->CmdListsCount == 0)
    {
        commandBuffer_ = VK_NULL_HANDLE;
        targetImageView_ = VK_NULL_HANDLE;
        return;
    }

    if (commandBuffer_ == VK_NULL_HANDLE || targetImageView_ == VK_NULL_HANDLE || renderExtent_.width == 0 || renderExtent_.height == 0)
    {
        commandBuffer_ = VK_NULL_HANDLE;
        targetImageView_ = VK_NULL_HANDLE;
        return;
    }

    VkRenderingAttachmentInfo colorAttachment = vkinit::attachment_info(targetImageView_, nullptr, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL);
    VkRenderingInfo renderInfo = vkinit::rendering_info(renderExtent_, &colorAttachment, nullptr);

    vkCmdBeginRendering(commandBuffer_, &renderInfo);
    ImGui_ImplVulkan_RenderDrawData(drawData, commandBuffer_);
    vkCmdEndRendering(commandBuffer_);

    commandBuffer_ = VK_NULL_HANDLE;
    targetImageView_ = VK_NULL_HANDLE;
}

void VulkanImguiBackend::shutdownBackend()
{
    if (!initialized_)
    {
        return;
    }

    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplSDL2_Shutdown();

    if (descriptorPool_ != VK_NULL_HANDLE)
    {
        vkDestroyDescriptorPool(device_, descriptorPool_, nullptr);
        descriptorPool_ = VK_NULL_HANDLE;
    }

    commandBuffer_ = VK_NULL_HANDLE;
    targetImageView_ = VK_NULL_HANDLE;
    initialized_ = false;
}

