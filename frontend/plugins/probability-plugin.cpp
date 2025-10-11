#include "probability-plugin.h"

#include "imgui.h"
#include "implot.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace math_gui::plugins
{
namespace
{
enum class DistributionType
{
    Continuous,
    Discrete
};

using DistributionVariant = std::variant<
    boost::math::normal_distribution<float>,
    boost::math::lognormal_distribution<float>>;

struct DistributionDefinition
{
    std::array<float, 2> domain;
    std::array<float, 2> renderDomain;
    DistributionType type;
    std::vector<float> parameters;
    std::vector<std::array<float, 2>> parameterDomains;
    std::function<DistributionVariant(const std::vector<float>&)> factory;
};

const std::array<const char*, 2> kDistributionNames = {
    "normal distribution",
    "lognormal distribution"
};

const std::unordered_map<std::string, DistributionDefinition> kDistributions = {
    {
        "normal distribution",
        DistributionDefinition{
            {-std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()},
            {-10.0f, 10.0f},
            DistributionType::Continuous,
            {0.0f, 1.0f},
            {{-std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()}, {1e-6f, std::numeric_limits<float>::infinity()}},
            [](const std::vector<float>& params) {
                const float mean = params.size() > 0 ? params[0] : 0.0f;
                const float stddev = params.size() > 1 ? std::max(params[1], 1e-6f) : 1.0f;
                return DistributionVariant{boost::math::normal_distribution<float>(mean, stddev)};
            }
        }
    },
    {
        "lognormal distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 10.0f},
            DistributionType::Continuous,
            {0.0f, 0.25f},
            {{-std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()}, {1e-6f, std::numeric_limits<float>::infinity()}},
            [](const std::vector<float>& params) {
                const float location = params.size() > 0 ? params[0] : 0.0f;
                const float scale = params.size() > 1 ? std::max(params[1], 1e-6f) : 0.25f;
                return DistributionVariant{boost::math::lognormal_distribution<float>(location, scale)};
            }
        }
    }
};

} // namespace

void RenderProbabilityWindow(ImGuiRenderer::FrameState& state)
{
    if (!state.show_probability_window)
    {
        return;
    }

    ImGui::Begin("Probability Window", &state.show_probability_window);

    static int item_selected_idx = 0;
    const char* combo_preview_value = kDistributionNames[item_selected_idx];

    if (ImGui::BeginCombo("Continuous Distributions", combo_preview_value))
    {
        static ImGuiTextFilter filter;
        if (ImGui::IsWindowAppearing())
        {
            ImGui::SetKeyboardFocusHere();
            filter.Clear();
        }
        ImGui::SetNextItemShortcut(ImGuiMod_Ctrl | ImGuiKey_F);
        filter.Draw("##Filter", -FLT_MIN);

        for (int n = 0; n < static_cast<int>(kDistributionNames.size()); ++n)
        {
            const bool is_selected = (item_selected_idx == n);
            if (filter.PassFilter(kDistributionNames[n]))
            {
                if (ImGui::Selectable(kDistributionNames[n], is_selected))
                {
                    item_selected_idx = n;
                }
            }
        }
        ImGui::EndCombo();
    }

    static bool should_render = false;
    ImGui::Checkbox("Render PDF/PMF", &should_render);

    if (should_render)
    {
        const auto* distribution_name = kDistributionNames[item_selected_idx];
        const auto& definition = kDistributions.at(distribution_name);
        const auto distribution = definition.factory(definition.parameters);

        static float xs1[1001];
        static float ys1[1001];
        constexpr int sample_count = static_cast<int>(std::size(xs1));

        ImPlot::SetNextAxesLimits(definition.renderDomain[0], definition.renderDomain[1], 0.0f, 2.0f, ImPlotCond_Once);
        if (ImPlot::BeginPlot("Line Plots"))
        {
            const ImPlotRect plot_limits = ImPlot::GetPlotLimits();
            float view_start = static_cast<float>(plot_limits.X.Min);
            float view_end = static_cast<float>(plot_limits.X.Max);
            if (!std::isfinite(view_start) || !std::isfinite(view_end))
            {
                view_start = definition.renderDomain[0];
                view_end = definition.renderDomain[1];
            }

            float sample_start = std::max(view_start, definition.domain[0]);
            float sample_end = std::min(view_end, definition.domain[1]);
            if (!std::isfinite(sample_start))
            {
                sample_start = view_start;
            }
            if (!std::isfinite(sample_end))
            {
                sample_end = view_end;
            }
            if (sample_end <= sample_start)
            {
                sample_end = sample_start + 1.0f;
            }

            const float sample_range = sample_end - sample_start;
            const float step = (sample_count > 1) ? sample_range / static_cast<float>(sample_count - 1) : 0.0f;

            for (int i = 0; i < sample_count; ++i)
            {
                const float x = sample_start + step * static_cast<float>(i);
                xs1[i] = x;
                ys1[i] = definition.type == DistributionType::Continuous
                             ? std::visit([x](const auto& dist) { return static_cast<float>(boost::math::pdf(dist, x)); }, distribution)
                             : 0.0f;
            }

            definition.type == DistributionType::Continuous
                ? ImPlot::PlotLine("PDF", xs1, ys1, sample_count)
                : ImPlot::PlotStairs("PMF", xs1, ys1, sample_count);
            ImPlot::EndPlot();
        }
    }

    ImGui::Text("Hello from probability window!");
    if (ImGui::Button("Close Me"))
    {
        state.show_probability_window = false;
    }
    ImGui::End();
}

} // namespace math_gui::plugins
