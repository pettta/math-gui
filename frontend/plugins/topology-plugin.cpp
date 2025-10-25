#include "topology-plugin.h"

#include "imgui.h"
#include "implot.h"

#include <cmath>
#include <set>
#include <optional>
#include <utility>
#include <vector>


namespace
{

constexpr int kMinCoordinate = 1;
constexpr int kMaxCoordinate = 9;
constexpr int kGridLength = kMaxCoordinate - kMinCoordinate + 1;
constexpr int kPointCount = kGridLength * kGridLength;

using PointPath = std::vector<ImPlotPoint>;
using GridCoord = std::pair<int, int>;
using VisitedSet = std::set<GridCoord>;

bool IsWithinGrid(int x, int y)
{
    return x >= kMinCoordinate && x <= kMaxCoordinate && y >= kMinCoordinate && y <= kMaxCoordinate;
}

std::optional<ImPlotPoint> AdvanceInCoordinate(const ImPlotPoint& point, const std::string& order_type, const VisitedSet& visited)
{
    const int rounded_x = static_cast<int>(std::round(point.x));
    const int rounded_y = static_cast<int>(std::round(point.y));

    if (!IsWithinGrid(rounded_x, rounded_y))
    {
        return std::nullopt;
    }

    if (visited.size() >= static_cast<std::size_t>(kPointCount))
    {
        return std::nullopt;
    }

    const auto is_unvisited = [&visited](int x, int y) {
        return visited.find(GridCoord{x, y}) == visited.end();
    };

    if (order_type == "dictionary")
    {
        for (int x = rounded_x; x <= kMaxCoordinate; ++x)
        {
            const int y_start = (x == rounded_x) ? rounded_y + 1 : kMinCoordinate;
            for (int y = y_start; y <= kMaxCoordinate; ++y)
            {
                if (!IsWithinGrid(x, y))
                {
                    continue;
                }

                if (is_unvisited(x, y))
                {
                    return ImPlotPoint(static_cast<double>(x), static_cast<double>(y));
                }
            }
        }
        return std::nullopt;
    }

    // (x0, y0) < (x1, y1) IF (x0-y0) < (x1-y1) OR ((x0-y0) == (x1-y1) AND x0 < x1)
    if (order_type == "difference")
    {
        const int current_diff = rounded_x - rounded_y;

        for (int y = rounded_y + 1; y <= kMaxCoordinate; ++y)
        {
            const int x_candidate = y + current_diff;
            if (!IsWithinGrid(x_candidate, y))
            {
                continue;
            }

            if (is_unvisited(x_candidate, y))
            {
                return ImPlotPoint(static_cast<double>(x_candidate), static_cast<double>(y));
            }
        }

        const int max_diff = kMaxCoordinate - kMinCoordinate;
        for (int diff = current_diff + 1; diff <= max_diff; ++diff)
        {
            for (int y = kMinCoordinate; y <= kMaxCoordinate; ++y)
            {
                const int x_candidate = y + diff;
                if (!IsWithinGrid(x_candidate, y))
                {
                    continue;
                }

                if (is_unvisited(x_candidate, y))
                {
                    return ImPlotPoint(static_cast<double>(x_candidate), static_cast<double>(y));
                }
            }
        }

        return std::nullopt;
    }

    // (x0, y0) < (x1, y1) IF (x0+y0) < (x1+y1) OR ((x0+y0) == (x1+y1) AND x0 < x1)
    if (order_type == "sum")
    {
        const int current_sum = rounded_x + rounded_y; 
        for (int y = rounded_y + 1; y <= kMaxCoordinate; ++y)
        {
            const int x_candidate = current_sum - y;
            if (!IsWithinGrid(x_candidate, y))
            {
                continue;
            }

            if (is_unvisited(x_candidate, y))
            {
                return ImPlotPoint(static_cast<double>(x_candidate), static_cast<double>(y));
            }
        }
        const int max_sum = 2 * kMaxCoordinate;
        for (int sum = current_sum + 1; sum <= max_sum; ++sum)
        {
            for (int y = kMinCoordinate; y <= kMaxCoordinate; ++y)
            {
                const int x_candidate = sum - y;
                if (!IsWithinGrid(x_candidate, y))
                {
                    continue;
                }

                if (is_unvisited(x_candidate, y))
                {
                    return ImPlotPoint(static_cast<double>(x_candidate), static_cast<double>(y));
                }
            }
        }
        return std::nullopt;
    }

    return std::nullopt;
}

PointPath NextPath(const ImPlotPoint& point, const std::string& order_type)
{
    PointPath path;
    VisitedSet visited;
    visited.emplace(static_cast<int>(std::round(point.x)), static_cast<int>(std::round(point.y)));

    std::optional<ImPlotPoint> current = AdvanceInCoordinate(point, order_type, visited);
    while (current)
    {
        path.push_back(*current);
        visited.emplace(static_cast<int>(std::round(current->x)), static_cast<int>(std::round(current->y)));
        current = AdvanceInCoordinate(*current, order_type, visited);
    }
    return path;
}


PointPath GetPathForActiveRelation(
    const ImPlotPoint& point,
    bool dictionary_active,
    bool difference_active,
    bool sum_active)
{
    if (dictionary_active)
    {
        return NextPath(point, "dictionary");
    }
    if (difference_active)
    {
        return NextPath(point, "difference");
    }
    if (sum_active)
    {
        return NextPath(point, "sum");
    }
    return {};
}

} // namespace


namespace math_gui::plugins
{
void RenderTopologyWindow(ImGuiRenderer::FrameState& state)
{
    ImGui::Begin("Topology Window", &state.show_topology_window);


    if (ImGui::CollapsingHeader("Munkres Practice Problem Visualizations")) {
        if (ImGui::TreeNode("1.3.12 Plots")){
            static float xs1[kPointCount];
            static float ys1[kPointCount];
            static bool points_initialized = false;

            if (!points_initialized)
            {
                int index = 0;
                for (int x = kMinCoordinate; x <= kMaxCoordinate; ++x)
                {
                    for (int y = kMinCoordinate; y <= kMaxCoordinate; ++y)
                    {
                        xs1[index] = static_cast<float>(x);
                        ys1[index] = static_cast<float>(y);
                        ++index;
                    }
                }
                points_initialized = true;
            }

            static bool dictionary_order_active = false;
            static bool difference_order_active = false;
            static bool sum_order_active = false;
            static std::optional<ImPlotPoint> selected_point;
            static PointPath relation_path;

            auto update_relation_path = [&]() {
                if (selected_point)
                {
                    relation_path = GetPathForActiveRelation(
                        *selected_point,
                        dictionary_order_active,
                        difference_order_active,
                        sum_order_active);
                }
                else
                {
                    relation_path.clear();
                }
            };

            const ImVec4 active_color(0.3f, 0.7f, 0.3f, 1.0f);
            const ImVec4 active_hover_color(0.34f, 0.78f, 0.34f, 1.0f);
            const ImVec4 active_pressed_color(0.26f, 0.62f, 0.26f, 1.0f);

            auto draw_order_button = [&](const char* label, bool is_active) {
                if (is_active)
                {
                    ImGui::PushStyleColor(ImGuiCol_Button, active_color);
                    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, active_hover_color);
                    ImGui::PushStyleColor(ImGuiCol_ButtonActive, active_pressed_color);
                }

                const bool clicked = ImGui::Button(label);

                if (is_active)
                {
                    ImGui::PopStyleColor(3);
                }

                return clicked;
            };

            if (draw_order_button("Dictionary Order", dictionary_order_active)) {
                dictionary_order_active = true;
                difference_order_active = false;
                sum_order_active = false;
                update_relation_path();
            }
            ImGui::SameLine();
            if (draw_order_button("Difference Order Relation", difference_order_active)) {
                dictionary_order_active = false;
                difference_order_active = true;
                sum_order_active = false;
                update_relation_path();
            }
            ImGui::SameLine();
            if (draw_order_button("Sum Order Relation", sum_order_active)) {
                dictionary_order_active = false;
                difference_order_active = false;
                sum_order_active = true;
                update_relation_path();
            }
            ImGui::Separator();

            if (!dictionary_order_active && !difference_order_active && !sum_order_active)
            {
                ImGui::TextUnformatted("Select an order relation to enable navigation arrows.");
            }

            const ImVec2 plot_size(680.0f, 680.0f);
            ImPlot::SetNextAxesLimits(
                static_cast<double>(kMinCoordinate) - 1.0,
                static_cast<double>(kMaxCoordinate) + 1.0,
                static_cast<double>(kMinCoordinate) - 1.0,
                static_cast<double>(kMaxCoordinate) + 1.0,
                ImPlotCond_Once);

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(18.0f, 18.0f));
            if (ImPlot::BeginPlot("Scatter Plot", plot_size)) {
                ImPlot::SetupAxes("x", "y");

                ImPlot::PushStyleVar(ImPlotStyleVar_FillAlpha, 0.25f);
                ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 6, ImPlot::GetColormapColor(1), IMPLOT_AUTO, ImPlot::GetColormapColor(1));
                ImPlot::PlotScatter("Grid Points", xs1, ys1, kPointCount);
                ImPlot::PopStyleVar();

                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
                {
                    const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
                    const int rounded_x = static_cast<int>(std::round(mouse.x));
                    const int rounded_y = static_cast<int>(std::round(mouse.y));
                    if (IsWithinGrid(rounded_x, rounded_y))
                    {
                        selected_point = ImPlotPoint(static_cast<double>(rounded_x), static_cast<double>(rounded_y));
                        update_relation_path();
                    }
                }

                if (selected_point)
                {
                    const double selected_x = selected_point->x;
                    const double selected_y = selected_point->y;
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 8.0f, ImVec4(1.0f, 0.9f, 0.6f, 0.35f), 2.0f, ImVec4(1.0f, 0.8f, 0.0f, 1.0f));
                    ImPlot::PlotScatter("##SelectedPoint", &selected_x, &selected_y, 1);

                    if (!relation_path.empty())
                    {
                        const ImVec4 line_color(1.0f, 1.0f, 1.0f, 1.0f);

                        ImDrawList* draw_list = ImPlot::GetPlotDrawList();
                        if (draw_list != nullptr)
                        {
                            ImPlot::PushPlotClipRect();
                        }

                        for (std::size_t idx = 0; idx < relation_path.size(); ++idx)
                        {
                            const ImPlotPoint& from = (idx == 0) ? *selected_point : relation_path[idx - 1];
                            const ImPlotPoint& to = relation_path[idx];

                            ImGui::PushID(static_cast<int>(idx));

                            const double seg_x[2] = {from.x, to.x};
                            const double seg_y[2] = {from.y, to.y};

                            ImPlot::SetNextLineStyle(line_color, 1.5f);
                            ImPlot::PlotLine("##Segment", seg_x, seg_y, 2);

                            if (draw_list != nullptr)
                            {
                                const ImVec2 from_pixels = ImPlot::PlotToPixels(from.x, from.y);
                                const ImVec2 to_pixels = ImPlot::PlotToPixels(to.x, to.y);
                                ImVec2 direction(to_pixels.x - from_pixels.x, to_pixels.y - from_pixels.y);
                                const float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
                                if (length > 1e-3f)
                                {
                                    direction.x /= length;
                                    direction.y /= length;
                                    const float arrow_size = 6.0f;
                                    const ImVec2 perpendicular(-direction.y, direction.x);
                                    const ImVec2 tip = to_pixels;
                                    const ImVec2 left = ImVec2(
                                        tip.x - direction.x * arrow_size + perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size + perpendicular.y * (arrow_size * 0.5f));
                                    const ImVec2 right = ImVec2(
                                        tip.x - direction.x * arrow_size - perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size - perpendicular.y * (arrow_size * 0.5f));
                                    const ImU32 arrow_color = IM_COL32(255, 255, 255, 255);
                                    draw_list->AddTriangleFilled(tip, left, right, arrow_color);
                                }
                            }

                            ImGui::PopID();
                        }

                        if (draw_list != nullptr)
                        {
                            ImPlot::PopPlotClipRect();
                        }

                    }
                    else
                    {
                        relation_path.clear();
                    }
                }

                ImPlot::EndPlot();
            }
            ImPlot::PopStyleVar();
            ImGui::TreePop();
            ImGui::Spacing();
       
        }

        if (ImGui::TreeNode("1.7.12 Plots")) {
            // Future topology visualizations can be added here
            ImGui::TextUnformatted("Visualizations for equivalence classes in countably infinite sets.");
            
            ImGui::TreePop();
            ImGui::Spacing();
        }
    }



    if (ImGui::Button("Close Me"))
        state.show_topology_window = false;
    ImGui::End();
}
} // namespace math_gui::plugins 