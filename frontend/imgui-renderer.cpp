#include "imgui-renderer.h"

#include <stdexcept>
#include <algorithm>
#include <array>
#include <functional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

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
    DistributionType type;
    std::vector<float> parameters;
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
            {-10.0f, 10.0f},
            DistributionType::Continuous,
            {0.0f, 1.0f},
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
            {0.0f, 10.0f},
            DistributionType::Continuous,
            {0.0f, 0.25f},
            [](const std::vector<float>& params) {
                const float location = params.size() > 0 ? params[0] : 0.0f;
                const float scale = params.size() > 1 ? std::max(params[1], 1e-6f) : 0.25f;
                return DistributionVariant{boost::math::lognormal_distribution<float>(location, scale)};
            }
        }
    }
};
} // namespace

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

        ImGui::Begin("Hello, world!");                                // Create a window called "Hello, world!" and append into it.

        ImGui::Text("This is some useful text.");                     // Display some text (you can use a format strings too)
        ImGui::Checkbox("Demo Window", &state.show_demo_window);      // Edit bools storing our window open/close state
        ImGui::Checkbox("Linear Algebra Window", &state.show_linear_algebra_window);
        ImGui::Checkbox("Probability Window", &state.show_probability_window);
        ImGui::ColorEdit3("clear color", state.clear_color);          // Edit 3 floats representing a color

        ImGuiIO& io = ImGui::GetIO(); 
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
        ImGui::End();
    }

    if (state.show_linear_algebra_window)
    {
        ImGui::Begin("Linear Algebra Window", &state.show_linear_algebra_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
        ImGui::Text("Hello from linear algebra window!");
        if (ImGui::Button("Close Me"))
            state.show_linear_algebra_window = false;
        ImGui::End();
    }

    if (state.show_probability_window)
    {
        
        ImGui::Begin("Probability Window", &state.show_probability_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
        
        static int item_selected_idx = 0; // Here we store our selection data as an index.
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

            for (int n = 0; n < static_cast<int>(kDistributionNames.size()); n++)
            {
                const bool is_selected = (item_selected_idx == n);
                if (filter.PassFilter(kDistributionNames[n]))
                    if (ImGui::Selectable(kDistributionNames[n], is_selected))
                        item_selected_idx = n;
            }
            ImGui::EndCombo();
        }

        static bool should_render = false;
        ImGui::Checkbox("Render PDF/PMF", &should_render);

        if (should_render) {
            const auto* distribution_name = kDistributionNames[item_selected_idx];
            const auto& definition = kDistributions.at(distribution_name);
            const auto distribution = definition.factory(definition.parameters);

            static float xs1[1001];
            static float ys1[1001];
            constexpr int sample_count = static_cast<int>(std::size(xs1));
            const float domain_start = definition.domain[0];
            const float domain_end = definition.domain[1];
            const float domain_range = domain_end - domain_start;
            const float domain_step = sample_count > 1 ? domain_range / static_cast<float>(sample_count - 1) : 0.0f;

            for (int i = 0; i < sample_count; ++i) {
            const float x = domain_start + domain_step * static_cast<float>(i);
            xs1[i] = x;
            ys1[i] = definition.type == DistributionType::Continuous
                     ? std::visit([x](const auto& dist) { return static_cast<float>(boost::math::pdf(dist, x)); }, distribution)
                     : 0.0f;
            }

            constexpr float y_min = 0.0f;
            constexpr float y_max = 2.0f;
            ImPlot::SetNextAxesLimits(domain_start, domain_end, y_min, y_max, ImPlotCond_Once);
            if (ImPlot::BeginPlot("Line Plots")) {
            definition.type == DistributionType::Continuous
                ? ImPlot::PlotLine("PDF", xs1, ys1, sample_count)
                : ImPlot::PlotStairs("PMF", xs1, ys1, sample_count);
            ImPlot::EndPlot();
            }
        }
        
        
        ImGui::Text("Hello from probability window!");
        if (ImGui::Button("Close Me"))
            state.show_probability_window = false;
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
    ImPlot::DestroyContext();
    ImPlot3D::DestroyContext();
    initialized_ = false;
}
