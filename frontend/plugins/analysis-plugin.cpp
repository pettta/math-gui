#include "analysis-plugin.h"

#include "imgui.h"
#include "implot.h"

#include <cmath>

namespace math_gui::plugins
{
void RenderAnalysisWindow(ImGuiRenderer::FrameState& state)
{
    ImGui::Begin("Analysis Window", &state.show_analysis_window);


    if (ImGui::CollapsingHeader("Stein Fourier Analysis")) {
        if (ImGui::TreeNode("Weierstrass Function")) {
            constexpr double kPi = 3.14159265358979323846;

            static float b_param = 0.5f;   // 0 < b < 1
            static int   a_param = 13;     // integer > 1; default ab = 6.5 > 5.712

            ImGui::TextUnformatted("W(x) = sum_{n=1}^{inf} b^n cos(a^n x)");
            ImGui::SliderFloat("b (0 < b < 1)", &b_param, 0.01f, 0.99f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
            ImGui::SliderInt("a (integer > 1)", &a_param, 2, 20, "%d", ImGuiSliderFlags_AlwaysClamp);

            const double ab = static_cast<double>(b_param) * static_cast<double>(a_param);
            const double threshold = 1.0 + 3.0 * kPi / 2.0;            // ~= 5.712389
            const bool condition_met = ab > threshold;
            const ImVec4 ok_color(0.2f, 0.8f, 0.3f, 1.0f);
            const ImVec4 bad_color(0.9f, 0.2f, 0.2f, 1.0f);
            ImGui::TextColored(condition_met ? ok_color : bad_color,
                "ab = %.3f  %s  1 + 3pi/2 = %.3f   ->   %s",
                ab, condition_met ? ">" : "<=", threshold,
                condition_met ? "condition satisfied" : "condition NOT satisfied");

            ImGui::Separator();

            constexpr int kSamples = 10000;
            static double wx_x[kSamples];
            static double wx_y[kSamples];

            const auto weierstrass = [](double x, double b, int a) {
                double sum = 0.0;
                double bn = b;                            // b^1
                double an = static_cast<double>(a);       // a^1
                for (int n = 1; n <= 100; ++n)            // hard cap on number of terms
                {
                    if (!std::isfinite(an))               // a^n overflowed -> cos(a^n x) is NaN
                    {
                        break;
                    }
                    sum += bn * std::cos(an * x);
                    if (bn < 1e-7)                        // remaining tail negligible (0 < b < 1)
                    {
                        break;
                    }
                    bn *= b;
                    an *= static_cast<double>(a);
                }
                return sum;
            };

            const double x_lo = -kPi;
            const double x_hi = kPi;
            for (int i = 0; i < kSamples; ++i)
            {
                const double x = x_lo + (x_hi - x_lo) * static_cast<double>(i) / static_cast<double>(kSamples - 1);
                wx_x[i] = x;
                wx_y[i] = weierstrass(x, static_cast<double>(b_param), a_param);
            }

            ImPlot::SetNextAxesLimits(-kPi, kPi, -2.0, 2.0, ImPlotCond_Once);
            if (ImPlot::BeginPlot("Weierstrass W(x)", ImVec2(-1, 400.0f)))
            {
                ImPlot::SetupAxes("x", "W(x)");
                ImPlot::PlotLine("W(x)", wx_x, wx_y, kSamples);
                ImPlot::EndPlot();
            }

            ImGui::TreePop();
            ImGui::Spacing();
        }
    }


    if (ImGui::Button("Close Me"))
        state.show_analysis_window = false;
    ImGui::End();
}
} // namespace math_gui::plugins
