#include "probability-plugin.h"

#include "imgui.h"
#include "implot.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <string>
#include <vector>



// === Boost Distributions === // 
// Continuous distributions 
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/exponential.hpp>

// Discrete distributions 
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>


// ===  Custom Distributions === // 
// Continuous distributions 
// Discrete distributions 
#include "../utils/distributions/beta_binomial.cpp"



namespace math_gui::plugins
{
namespace
{
enum class DistributionType
{
    Continuous,
    Discrete
};

using PdfFunction = std::function<float(float)>;
using PdfFactory = std::function<PdfFunction(const std::vector<float>&)>;
using CdfFunction = std::function<float(float)>;
using CdfFactory = std::function<CdfFunction(const std::vector<float>&)>;

struct DistributionDefinition
{
    std::array<float, 2> domain;
    std::array<float, 2> renderDomain;
    DistributionType type;
    std::vector<float> parameters;
    std::vector<std::array<float, 2>> parameterDomains;
    std::vector<std::string> parameterNames;
    std::vector<std::string> parameterRanges;
    std::vector<std::string> parameterDescriptions;
    PdfFactory makePdf;
    std::vector<bool> parameterIsIntegral;
    CdfFactory makeCdf;
};

struct DistributionEntry
{
    std::string label;
    DistributionDefinition definition;
};

const std::vector<DistributionEntry> kDistributions = {
    {
        "normal distribution",
        DistributionDefinition{
            {-std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()},
            {-10.0f, 10.0f},
            DistributionType::Continuous,
            {0.0f, 1.0f},
            {{-std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()}, {1e-6f, std::numeric_limits<float>::infinity()}},
            {(const char*)u8"μ", "σ"},
            {(const char*)u8"(-∞, ∞)", "σ   > 0"},
            {"Mean", "Standard deviation"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float mean = params.size() > 0 ? params[0] : 0.0f;
                const float stddev = params.size() > 1 ? std::max(params[1], 1e-6f) : 1.0f;
                const boost::math::normal_distribution<float> distribution(mean, stddev);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, x));
                };
            },
            {false, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float mean = params.size() > 0 ? params[0] : 0.0f;
                const float stddev = params.size() > 1 ? std::max(params[1], 1e-6f) : 1.0f;
                const boost::math::normal_distribution<float> distribution(mean, stddev);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, x));
                };
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
            {"μ", "σ"},
            {"(-∞, ∞)", "σ > 0"},
            {"Location", "Scale"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float location = params.size() > 0 ? params[0] : 0.0f;
                const float scale = params.size() > 1 ? std::max(params[1], 1e-6f) : 0.25f;
                const boost::math::lognormal_distribution<float> distribution(location, scale);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, x));
                };
            },
            {false, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float location = params.size() > 0 ? params[0] : 0.0f;
                const float scale = params.size() > 1 ? std::max(params[1], 1e-6f) : 0.25f;
                const boost::math::lognormal_distribution<float> distribution(location, scale);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, x));
                };
            }
        }
    },
    {
        "beta distribution",
        DistributionDefinition{
            {0.0f, 1.0f},
            {0.0f, 1.0f},
            DistributionType::Continuous,
            {2.0f, 5.0f},
            {{1e-6f, std::numeric_limits<float>::infinity()}, {1e-6f, std::numeric_limits<float>::infinity()}},
            {"α", "β"},
            {"α > 0", "β > 0"},
            {"First shape", "Second shape"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float alpha = params.size() > 0 ? std::max(params[0], 1e-6f) : 2.0f;
                const float beta = params.size() > 1 ? std::max(params[1], 1e-6f) : 5.0f;
                const boost::math::beta_distribution<float> distribution(alpha, beta);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, x));
                };
            },
            {false, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float alpha = params.size() > 0 ? std::max(params[0], 1e-6f) : 2.0f;
                const float beta = params.size() > 1 ? std::max(params[1], 1e-6f) : 5.0f;
                const boost::math::beta_distribution<float> distribution(alpha, beta);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, x));
                };
            }
        }
    },
    {
        "exponential distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 10.0f},
            DistributionType::Continuous,
            {1.0f},
            {{1e-6f, std::numeric_limits<float>::infinity()}},
            {"λ"},
            {"λ > 0"},
            {"Rate parameter"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float lambda = params.size() > 0 ? std::max(params[0], 1e-6f) : 1.0f;
                const boost::math::exponential_distribution<float> distribution(lambda);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, x));
                };
            },
            {false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float lambda = params.size() > 0 ? std::max(params[0], 1e-6f) : 1.0f;
                const boost::math::exponential_distribution<float> distribution(lambda);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, x));
                };
            }
        }
    },
    {
        "binomial distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 100.0f},
            DistributionType::Discrete,
            {10.0f, 0.5f},
            {{1e-6f, std::numeric_limits<float>::infinity()}, {0.0f, 1.0f}},
            {"n", "p"},
            {"n ≥ 1", "0 ≤ p ≤ 1"},
            {"Number of trials", "Probability of success"},
            [](const std::vector<float>& params) -> PdfFunction {
                const int trials = params.size() > 0 ? static_cast<int>(std::max(params[0], 1e-6f)) : 10;
                const float p = params.size() > 1 ? std::clamp(params[1], 0.0f, 1.0f) : 0.5f;
                const boost::math::binomial_distribution<float> distribution(trials, p);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, std::round(x)));
                };
            },
            {true, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const int trials = params.size() > 0 ? static_cast<int>(std::max(params[0], 1e-6f)) : 10;
                const float p = params.size() > 1 ? std::clamp(params[1], 0.0f, 1.0f) : 0.5f;
                const boost::math::binomial_distribution<float> distribution(trials, p);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, std::round(x)));
                };
            }
        }
    },
    {
        "beta binomial distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 100.0f},
            DistributionType::Discrete,
            {10.0f, 2.0f, 5.0f},
            {{1.0f, 200.0f}, {1e-5f, std::numeric_limits<float>::infinity()}, {1e-5f, std::numeric_limits<float>::infinity()}},
            {"n", "α", "β"},
            {"n ≥ 1", "α > 0", "β > 0"},
            {"Trials", "Beta prior α", "Beta prior β"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float trials_value = params.size() > 0 ? std::max(params[0], 1.0f) : 10.0f;
                const int trials = static_cast<int>(std::round(trials_value));
                const float alpha = params.size() > 1 ? std::max(params[1], 1e-5f) : 2.0f;
                const float beta = params.size() > 2 ? std::max(params[2], 1e-5f) : 5.0f;
                const boost::math::beta_binomial_distribution<float> distribution(trials, alpha, beta);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, std::round(std::clamp(x, 0.0f, static_cast<float>(distribution.trials())))));
                };
            },
            {true, false, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float trials_value = params.size() > 0 ? std::max(params[0], 1.0f) : 10.0f;
                const int trials = static_cast<int>(std::round(trials_value));
                const float alpha = params.size() > 1 ? std::max(params[1], 1e-5f) : 2.0f;
                const float beta = params.size() > 2 ? std::max(params[2], 1e-5f) : 5.0f;
                const boost::math::beta_binomial_distribution<float> distribution(trials, alpha, beta);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, std::round(std::clamp(x, 0.0f, static_cast<float>(distribution.trials())))));
                };
            }
        }
    },
    {
        "negative binomial distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 100.0f},
            DistributionType::Discrete,
            {10.0f, 0.5f},
            {{1e-6f, std::numeric_limits<float>::infinity()}, {0.0f, 1.0f}},
            {"r", "p"},
            {"r > 0", "0 ≤ p ≤ 1"},
            {"Failures until stop", "Success probability"},
            [](const std::vector<float>& params) -> PdfFunction {
                const int k = params.size() > 0 ? static_cast<int>(std::max(params[0], 1e-6f)) : 10;
                const float p = params.size() > 1 ? std::clamp(params[1], 0.0f, 1.0f) : 0.5f;
                const boost::math::negative_binomial_distribution<float> distribution(k, p);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, std::round(x)));
                };
            },
            {true, false},
            [](const std::vector<float>& params) -> CdfFunction {
                const int k = params.size() > 0 ? static_cast<int>(std::max(params[0], 1e-6f)) : 10;
                const float p = params.size() > 1 ? std::clamp(params[1], 0.0f, 1.0f) : 0.5f;
                const boost::math::negative_binomial_distribution<float> distribution(k, p);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, std::round(x)));
                };
            }
        }
    },
    {
        "poisson distribution",
        DistributionDefinition{
            {0.0f, std::numeric_limits<float>::infinity()},
            {0.0f, 100.0f},
            DistributionType::Discrete,
            {10.0f},
            {{1e-6f, std::numeric_limits<float>::infinity()}},
            {"λ"},
            {"λ > 0"},
            {"Mean events"},
            [](const std::vector<float>& params) -> PdfFunction {
                const float mean = params.size() > 0 ? std::max(params[0], 1e-6f) : 10.0f;
                const boost::math::poisson_distribution<float> distribution(mean);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::pdf(distribution, std::round(x)));
                };
            },
            {false},
            [](const std::vector<float>& params) -> CdfFunction {
                const float mean = params.size() > 0 ? std::max(params[0], 1e-6f) : 10.0f;
                const boost::math::poisson_distribution<float> distribution(mean);
                return [distribution](float x) -> float {
                    return static_cast<float>(boost::math::cdf(distribution, std::round(x)));
                };
            }
        }
    },
};

std::vector<std::vector<float>> gDistributionParameters = [] {
    std::vector<std::vector<float>> values;
    values.reserve(kDistributions.size());
    for (const auto& entry : kDistributions)
    {
        values.push_back(entry.definition.parameters);
    }
    return values;
}();

} // namespace

void RenderProbabilityWindow(ImGuiRenderer::FrameState& state)
{
    if (!state.show_probability_window)
    {
        return;
    }

    ImGui::Begin("Probability Window", &state.show_probability_window);

    static int item_selected_idx = 0;
    item_selected_idx = std::clamp(item_selected_idx, 0, static_cast<int>(kDistributions.size()) - 1);
    const auto& current_entry = kDistributions[item_selected_idx];

    const ImVec4 continuous_color(0.2f, 0.8f, 0.3f, 1.0f);
    const ImVec4 discrete_color(0.9f, 0.2f, 0.2f, 1.0f);
    const ImVec4 preview_color = current_entry.definition.type == DistributionType::Continuous ? continuous_color : discrete_color;
    ImGui::PushStyleColor(ImGuiCol_Text, preview_color);
    const bool opened = ImGui::BeginCombo("Distributions", current_entry.label.c_str());
    ImGui::PopStyleColor();

    if (opened)
    {
        static ImGuiTextFilter filter;
        if (ImGui::IsWindowAppearing())
        {
            ImGui::SetKeyboardFocusHere();
            filter.Clear();
        }
        ImGui::SetNextItemShortcut(ImGuiMod_Ctrl | ImGuiKey_F);
        filter.Draw("##Filter", -FLT_MIN);

        for (int n = 0; n < static_cast<int>(kDistributions.size()); ++n)
        {
            const bool is_selected = (item_selected_idx == n);
            const auto& entry = kDistributions[n];
            const ImVec4 entry_color = entry.definition.type == DistributionType::Continuous ? continuous_color : discrete_color;
            ImGui::PushStyleColor(ImGuiCol_Text, entry_color);
            const bool passes_filter = filter.PassFilter(entry.label.c_str());
            if (passes_filter)
            {
                if (ImGui::Selectable(entry.label.c_str(), is_selected))
                {
                    item_selected_idx = n;
                }
            }
            ImGui::PopStyleColor();
        }
        ImGui::EndCombo();
    }

    static bool should_render = false;
    ImGui::Checkbox("Render PDF/PMF", &should_render);

    const auto& definition = kDistributions[item_selected_idx].definition;
    auto& parameter_values = gDistributionParameters[item_selected_idx];

    if (ImGui::BeginTable("ParameterTable", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_BordersOuter))
    {
        ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Range", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Description", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableHeadersRow();

        for (std::size_t idx = 0; idx < parameter_values.size(); ++idx)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            const std::string& name = idx < definition.parameterNames.size() ? definition.parameterNames[idx] : ("Param " + std::to_string(idx + 1));
            ImGui::TextUnformatted(name.c_str());

            ImGui::TableSetColumnIndex(1);
            const std::string& range = idx < definition.parameterRanges.size() ? definition.parameterRanges[idx] : "";
            ImGui::TextUnformatted(range.c_str());

            ImGui::TableSetColumnIndex(2);
            const std::string& description = idx < definition.parameterDescriptions.size() ? definition.parameterDescriptions[idx] : "";
            ImGui::TextWrapped("%s", description.c_str());

            ImGui::TableSetColumnIndex(3);
            const auto& bounds = definition.parameterDomains[idx];
            const bool is_integral = idx < definition.parameterIsIntegral.size() ? definition.parameterIsIntegral[idx] : false;
            const float min_bound = bounds[0];
            const float max_bound = bounds[1];

            ImGui::PushID(static_cast<int>(idx));
            float step = std::max((max_bound - min_bound) * 0.005f, 0.001f);
            if (is_integral)
            {
                step = 1.0f;
            }
            if (!std::isfinite(step) || step <= 0.0f)
            {
                step = is_integral ? 1.0f : 0.1f;
            }

            if (ImGui::DragFloat("##ParamValue", &parameter_values[idx], step, min_bound, max_bound, "%.4f", ImGuiSliderFlags_AlwaysClamp))
            {
                parameter_values[idx] = std::clamp(parameter_values[idx], min_bound, max_bound);
                if (is_integral)
                {
                    parameter_values[idx] = std::round(parameter_values[idx]);
                }
            }
            else
            {
                parameter_values[idx] = std::clamp(parameter_values[idx], min_bound, max_bound);
                if (is_integral)
                {
                    parameter_values[idx] = std::round(parameter_values[idx]);
                }
            }
            if (is_integral)
            {
                ImGui::SameLine();
                ImGui::Text("%d", static_cast<int>(parameter_values[idx]));
            }
            ImGui::PopID();
        }

        ImGui::EndTable();
    }

    if (should_render)
    {
        const auto pdf = definition.makePdf(parameter_values);
        const auto cdf = definition.makeCdf(parameter_values);

        static float xs1[1001];
        static float ys1[1001];
        static float ys2[1001];
        constexpr int sample_count = static_cast<int>(std::size(xs1));

    ImPlot::SetNextAxesLimits(definition.renderDomain[0], definition.renderDomain[1], 0.0f, 2.0f, ImPlotCond_Once);
    ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.2f, 0.6f, 0.9f, 1.0f));
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

            auto adjust_sample_to_domain = [](float value, float lower, float upper) {
                if (std::isfinite(lower) && value <= lower)
                {
                    const float toward = std::isfinite(upper) ? upper : std::numeric_limits<float>::infinity();
                    value = std::nextafter(lower, toward);
                }
                if (std::isfinite(upper) && value >= upper)
                {
                    const float toward = std::isfinite(lower) ? lower : -std::numeric_limits<float>::infinity();
                    value = std::nextafter(upper, toward);
                }
                return value;
            };

            for (int i = 0; i < sample_count; ++i)
            {
                float x = sample_start + step * static_cast<float>(i);
                x = adjust_sample_to_domain(x, definition.domain[0], definition.domain[1]);
                xs1[i] = x;
                if (definition.type == DistributionType::Continuous)
                {
                    try
                    {
                        ys1[i] = pdf(x);
                        if (!std::isfinite(ys1[i]))
                        {
                            ys1[i] = 0.0f;
                        }
                        ys2[i] = cdf(x);
                        if (!std::isfinite(ys2[i]))
                        {
                            ys2[i] = 0.0f;
                        }
                    }
                    catch (const std::exception&)
                    {
                        ys1[i] = 0.0f;
                        ys2[i] = 0.0f;
                    }
                }
                else if (definition.type == DistributionType::Discrete)
                {
                    try
                    {
                        const float x_rounded = std::round(x);
                        ys1[i] = pdf(x_rounded);
                        if (!std::isfinite(ys1[i]))
                        {
                            ys1[i] = 0.0f;
                        }
                        ys2[i] = cdf(x_rounded);
                        if (!std::isfinite(ys2[i]))
                        {
                            ys2[i] = 0.0f;
                        }
                    }
                    catch (const std::exception&)
                    {
                        ys1[i] = 0.0f;
                        ys2[i] = 0.0f;
                    }
                }
                else
                {
                    ys1[i] = 0.0f;
                    ys2[i] = 0.0f;
                }
            }

            definition.type == DistributionType::Continuous
                ? ImPlot::PlotLine("PDF", xs1, ys1, sample_count)
                : ImPlot::PlotStairs("PMF", xs1, ys1, sample_count);

            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.9f, 0.4f, 0.2f, 1.0f));
            definition.type == DistributionType::Continuous
                ? ImPlot::PlotLine("CDF", xs1, ys2, sample_count)
                : ImPlot::PlotStairs("CDF", xs1, ys2, sample_count);
            ImPlot::PopStyleColor();
            ImPlot::EndPlot();
        }
        ImPlot::PopStyleColor();
    }

    if (ImGui::Button("Close Me"))
    {
        state.show_probability_window = false;
    }
    ImGui::End();
}

} // namespace math_gui::plugins
