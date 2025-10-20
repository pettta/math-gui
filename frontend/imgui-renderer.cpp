#include "imgui-renderer.h"

#include "plugins/probability-plugin.h"
#include "plugins/topology-plugin.h" 

#include <stdexcept>
#include <array>
#include <filesystem>
#include <iostream>
#include <system_error>

ImGuiRenderer::ImGuiRenderer(ImGuiBackend& backend)
    : backend_(backend)
{
}

ImGuiRenderer::~ImGuiRenderer()
{
    if (initialized_)
    {
        shutdown();
    }
}

void ImGuiRenderer::initialize()
{
    if (initialized_)
    {
        return;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImPlot3D::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;

    ImGui::StyleColorsDark();
    ImFont* builtin_font = io.Fonts->AddFontDefault();


    ImFont* noto_font = nullptr;
    std::filesystem::path loaded_noto_font_path;
    noto_font = io.Fonts->AddFontFromFileTTF("../frontend/utils/fonts/NotoSansMath-Regular.ttf", 18.0f);


    // Set default font BEFORE backend initialization
    io.FontDefault = noto_font ? noto_font : builtin_font;

    if (!noto_font)
    {
        std::cout << "[ImGuiRenderer] Failed to load Noto Sans font. Using built-in default font." << std::endl;
    }
    backend_.initializeBackend();
    initialized_ = true;
}

ImGuiIO& ImGuiRenderer::io()
{
    if (!initialized_)
    {
        throw std::runtime_error("ImGuiRenderer::io called before initialize()");
    }
    return ImGui::GetIO();
}

ImGuiRenderer::FrameState& ImGuiRenderer::frameState()
{
    return initialState_;
}

const ImGuiRenderer::FrameState& ImGuiRenderer::frameState() const
{
    return initialState_;
}

void ImGuiRenderer::beginFrame()
{
    if (!initialized_)
    {
        initialize();
    }

    backend_.newFrame();
    ImGui::NewFrame();
}

void ImGuiRenderer::businessLogic(FrameState& state)
{
    // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
    if (state.show_demo_window){
        ImGui::ShowDemoWindow(&state.show_demo_window);
        ImPlot::ShowDemoWindow(&state.show_demo_window);
        ImPlot3D::ShowDemoWindow(&state.show_demo_window);
    }
        

    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
    {

        ImGui::Begin("Hello, world!");      
        ImGui::Checkbox("Demo Window", &state.show_demo_window);      
        ImGui::Checkbox("Linear Algebra Window", &state.show_linear_algebra_window);
        ImGui::Checkbox("Probability Window", &state.show_probability_window);
        ImGui::Checkbox("Topology Window", &state.show_topology_window); 
        ImGui::ColorEdit3("clear color", state.clear_color);          // Edit 3 floats representing a color

        ImGuiIO& io = ImGui::GetIO(); 
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
        ImGui::End();
    }

    if (state.show_linear_algebra_window)
    {
        // TODO IMPLEMENT 
        ImGui::Begin("Linear Algebra Window", &state.show_linear_algebra_window);  
        ImGui::Text("Hello from linear algebra window!");
        if (ImGui::Button("Close Me"))
            state.show_linear_algebra_window = false;
        ImGui::End();
    }

    if (state.show_topology_window) 
    {
        math_gui::plugins::RenderTopologyWindow(state); 
    }

    if (state.show_probability_window)
    {
        math_gui::plugins::RenderProbabilityWindow(state);
    }
}

void ImGuiRenderer::endFrame()
{
    if (!initialized_)
    {
        return;
    }

    ImGui::Render();
    backend_.renderDrawData(ImGui::GetDrawData());
}

void ImGuiRenderer::renderFrame(const std::function<void(FrameState&)>& drawUI)
{
    beginFrame();
    if (drawUI)
    {
        drawUI(initialState_);
    }
    else
    {
        businessLogic(initialState_);
    }
    endFrame();
}

void ImGuiRenderer::shutdown()
{
    if (!initialized_)
    {
        return;
    }

    backend_.shutdownBackend();
    ImGui::DestroyContext();
    ImPlot::DestroyContext();
    ImPlot3D::DestroyContext();
    initialized_ = false;
}
