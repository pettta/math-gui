#include "imgui-renderer.h"

#include <stdexcept>

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
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;

    ImGui::StyleColorsDark();

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
    if (state.show_demo_window)
        ImGui::ShowDemoWindow(&state.show_demo_window);

    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
    {
        static float f = 0.0f;
        static int counter = 0;

        ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

        ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
    ImGui::Checkbox("Demo Window", &state.show_demo_window);      // Edit bools storing our window open/close state
    ImGui::Checkbox("Another Window", &state.show_another_window);

    ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
    ImGui::ColorEdit3("clear color", state.clear_color); // Edit 3 floats representing a color

        if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
            counter++;
        ImGui::SameLine();
        ImGui::Text("counter = %d", counter);

        ImGuiIO& io = ImGui::GetIO(); 
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
        ImGui::End();
    }

    // 3. Show another simple window.
    if (state.show_another_window)
    {
        ImGui::Begin("Another Window", &state.show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
        ImGui::Text("Hello from another window!");
        if (ImGui::Button("Close Me"))
            state.show_another_window = false;
        ImGui::End();
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
    initialized_ = false;
}
