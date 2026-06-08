#include "analysis-plugin.h"

#include "imgui.h"

namespace math_gui::plugins
{
void RenderAnalysisWindow(ImGuiRenderer::FrameState& state)
{
    ImGui::Begin("Analysis Window", &state.show_analysis_window);
    ImGui::Text("Hello from analysis window!");
    if (ImGui::Button("Close Me"))
        state.show_analysis_window = false;
    ImGui::End();
}
} // namespace math_gui::plugins
