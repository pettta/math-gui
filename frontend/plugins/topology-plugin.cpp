#include "topology-plugin.h"

#include "imgui.h"
#include "implot.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <limits>
#include <set>
#include <optional>
#include <utility>
#include <vector>
#include <numeric>


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

        if (ImGui::TreeNode("1.7.1 Plots")) {
            constexpr int kPairGridMax = 10;
            constexpr int kMinPairSum = 2;
            constexpr int kMaxPairSum = kPairGridMax * 2;
            constexpr int kClassCount = kMaxPairSum - kMinPairSum + 1;

            struct ClassVisualization
            {
                std::vector<double> pair_x;
                std::vector<double> pair_y;
                std::vector<double> z_x;
                std::vector<double> z_y;
                std::vector<double> q_x;
                std::vector<double> q_y;
            };

            static bool equivalence_data_initialized = false;
            static std::vector<double> pair_grid_x;
            static std::vector<double> pair_grid_y;
            static std::array<ClassVisualization, kClassCount> equivalence_classes{};
            static std::array<int, kClassCount> class_sizes{};
            static std::array<ImVec4, kClassCount> class_colors{};
            static double q_min_value = std::numeric_limits<double>::infinity();
            static double q_max_value = -std::numeric_limits<double>::infinity();

            if (!equivalence_data_initialized)
            {
                pair_grid_x.reserve(kPairGridMax * kPairGridMax);
                pair_grid_y.reserve(kPairGridMax * kPairGridMax);

                for (int x = 1; x <= kPairGridMax; ++x)
                {
                    for (int y = 1; y <= kPairGridMax; ++y)
                    {
                        pair_grid_x.push_back(static_cast<double>(x));
                        pair_grid_y.push_back(static_cast<double>(y));

                        const int sum = x + y;
                        const int idx = sum - kMinPairSum;
                        ClassVisualization& storage = equivalence_classes.at(idx);
                        storage.pair_x.push_back(static_cast<double>(x));
                        storage.pair_y.push_back(static_cast<double>(y));

                        double q_value = 0.0;
                        if (y % 2 == 0)
                        {
                            q_value = static_cast<double>(y - 1) / static_cast<double>(x);
                        }
                        else
                        {
                            q_value = -static_cast<double>(y) / static_cast<double>(x);
                        }

                        storage.q_x.push_back(q_value);
                        storage.q_y.push_back(0.0);
                        q_min_value = std::min(q_min_value, q_value);
                        q_max_value = std::max(q_max_value, q_value);
                    }
                }

                const int colormap_size = ImPlot::GetColormapSize();
                for (int class_idx = 0; class_idx < kClassCount; ++class_idx)
                {
                    ClassVisualization& storage = equivalence_classes[class_idx];
                    class_sizes[class_idx] = static_cast<int>(storage.pair_x.size());

                    const double column = static_cast<double>(class_idx + 1);
                    storage.z_x.resize(storage.pair_x.size(), column);
                    storage.z_y.resize(storage.pair_x.size());
                    for (int member_idx = 0; member_idx < class_sizes[class_idx]; ++member_idx)
                    {
                        storage.z_y[member_idx] = static_cast<double>(member_idx + 1);
                    }

                    const ImVec4 color = ImPlot::GetColormapColor(class_idx % colormap_size);
                    class_colors[class_idx] = color;

                    storage.q_y.assign(storage.q_x.size(), 0.0);
                }

                equivalence_data_initialized = true;
            }

            const ImVec2 subplots_size(-1, 340.0f);
            if (ImPlot::BeginSubplots("##1_7_12_subplots", 1, 3, subplots_size))
            {
                if (ImPlot::BeginPlot("Q Representation", ImVec2(-1, -1)))
                {
                    const double q_min = std::isfinite(q_min_value) ? q_min_value : -1.0;
                    const double q_max = std::isfinite(q_max_value) ? q_max_value : 1.0;
                    const double x_min = (q_min == q_max) ? q_min - 1.0 : q_min - 0.5;
                    const double x_max = (q_min == q_max) ? q_max + 1.0 : q_max + 0.5;

                    ImPlot::SetupAxes("q", "Stack Height");
                    ImPlot::SetupAxesLimits(x_min, x_max, -0.5, 0.5, ImPlotCond_Once);

                    ImDrawList* draw_list = ImPlot::GetPlotDrawList();
                    if (draw_list != nullptr)
                    {
                        ImPlot::PushPlotClipRect();
                    }

                    char label_buffer[32];
                    for (int class_idx = 0; class_idx < kClassCount; ++class_idx)
                    {
                        const auto& storage = equivalence_classes[class_idx];
                        if (storage.q_x.empty())
                        {
                            continue;
                        }

                        const ImVec4 color = class_colors[class_idx];
                        const ImVec4 fill_color(color.x, color.y, color.z, 0.25f);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 7.0f, fill_color, 1.5f, color);
                        std::snprintf(label_buffer, sizeof(label_buffer), "q-class %d", class_idx + 1);
                        ImPlot::PlotScatter(label_buffer, storage.q_x.data(), storage.q_y.data(), class_sizes[class_idx]);

                        if (draw_list != nullptr && storage.q_x.size() >= 2)
                        {
                            ImGui::PushID(class_idx);
                            for (int member_idx = 1; member_idx < class_sizes[class_idx]; ++member_idx)
                            {
                                ImGui::PushID(member_idx);
                                const ImVec2 start_pixels = ImPlot::PlotToPixels(
                                    storage.q_x[member_idx - 1],
                                    storage.q_y[member_idx - 1]);
                                const ImVec2 end_pixels = ImPlot::PlotToPixels(
                                    storage.q_x[member_idx],
                                    storage.q_y[member_idx]);

                                ImVec2 delta(end_pixels.x - start_pixels.x, end_pixels.y - start_pixels.y);
                                const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
                                if (length > 1e-3f)
                                {
                                    ImVec2 normal(-delta.y, delta.x);
                                    normal.x /= length;
                                    normal.y /= length;
                                    const float curvature = std::min(length * 0.3f, 40.0f);
                                    const ImVec2 control1 = ImVec2(
                                        start_pixels.x + delta.x * 0.33f + normal.x * curvature,
                                        start_pixels.y + delta.y * 0.33f + normal.y * curvature);
                                    const ImVec2 control2 = ImVec2(
                                        start_pixels.x + delta.x * 0.66f + normal.x * curvature,
                                        start_pixels.y + delta.y * 0.66f + normal.y * curvature);

                                    const ImU32 color_u32 = ImGui::GetColorU32(color);
                                    draw_list->AddBezierCubic(start_pixels, control1, control2, end_pixels, color_u32, 1.5f);

                                    ImVec2 tangent(end_pixels.x - control2.x, end_pixels.y - control2.y);
                                    const float tangent_length = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y);
                                    if (tangent_length > 1e-3f)
                                    {
                                        tangent.x /= tangent_length;
                                        tangent.y /= tangent_length;
                                        const float arrow_size = 8.0f;
                                        const ImVec2 perpendicular(-tangent.y, tangent.x);
                                        const ImVec2 tip = end_pixels;
                                        const ImVec2 left = ImVec2(
                                            tip.x - tangent.x * arrow_size + perpendicular.x * (arrow_size * 0.5f),
                                            tip.y - tangent.y * arrow_size + perpendicular.y * (arrow_size * 0.5f));
                                        const ImVec2 right = ImVec2(
                                            tip.x - tangent.x * arrow_size - perpendicular.x * (arrow_size * 0.5f),
                                            tip.y - tangent.y * arrow_size - perpendicular.y * (arrow_size * 0.5f));
                                        draw_list->AddTriangleFilled(tip, left, right, color_u32);
                                    }
                                }
                                ImGui::PopID();
                            }
                            ImGui::PopID();
                        }
                    }

                    if (draw_list != nullptr)
                    {
                        ImPlot::PopPlotClipRect();
                    }

                    ImPlot::EndPlot();
                }

                if (ImPlot::BeginPlot("Z^{+} \u00D7 Z^{+}", ImVec2(-1, -1)))
                {
                    ImPlot::SetupAxes("x", "y");
                    ImPlot::SetupAxesLimits(0.5, static_cast<double>(kPairGridMax) + 0.5, 0.5, static_cast<double>(kPairGridMax) + 0.5, ImPlotCond_Once);

                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 3.0f, ImVec4(0.7f, 0.7f, 0.9f, 0.9f));
                    ImPlot::PlotScatter("Points", pair_grid_x.data(), pair_grid_y.data(), static_cast<int>(pair_grid_x.size()));

                    ImDrawList* draw_list = ImPlot::GetPlotDrawList();
                    if (draw_list != nullptr)
                    {
                        ImPlot::PushPlotClipRect();
                    }

                    char label_buffer[32];
                    for (int class_idx = 0; class_idx < kClassCount; ++class_idx)
                    {
                        const auto& storage = equivalence_classes[class_idx];
                        if (storage.pair_x.empty())
                        {
                            continue;
                        }

                        const ImVec4 color = class_colors[class_idx];
                        const ImVec4 fill_color(color.x, color.y, color.z, 0.08f);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 10.0f, fill_color, 2.0f, color);
                        std::snprintf(label_buffer, sizeof(label_buffer), "Class %d", class_idx + 1);
                        ImPlot::PlotScatter(label_buffer, storage.pair_x.data(), storage.pair_y.data(), class_sizes[class_idx]);

                        if (draw_list != nullptr && class_sizes[class_idx] >= 2)
                        {
                            ImGui::PushID(class_idx);
                            std::vector<int> indices(class_sizes[class_idx]);
                            std::iota(indices.begin(), indices.end(), 0);
                            std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
                                if (storage.pair_x[lhs] == storage.pair_x[rhs])
                                {
                                    return storage.pair_y[lhs] < storage.pair_y[rhs];
                                }
                                return storage.pair_x[lhs] > storage.pair_x[rhs];
                            });

                            for (std::size_t member_idx = 1; member_idx < indices.size(); ++member_idx)
                            {
                                ImGui::PushID(static_cast<int>(member_idx));
                                const int from_idx = indices[member_idx - 1];
                                const int to_idx = indices[member_idx];

                                const ImVec2 start_pixels = ImPlot::PlotToPixels(
                                    storage.pair_x[from_idx],
                                    storage.pair_y[from_idx]);
                                const ImVec2 end_pixels = ImPlot::PlotToPixels(
                                    storage.pair_x[to_idx],
                                    storage.pair_y[to_idx]);

                                draw_list->AddLine(start_pixels, end_pixels, ImGui::GetColorU32(color), 1.5f);

                                ImVec2 direction(end_pixels.x - start_pixels.x, end_pixels.y - start_pixels.y);
                                const float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
                                if (length > 1e-3f)
                                {
                                    direction.x /= length;
                                    direction.y /= length;
                                    const float arrow_size = 8.0f;
                                    const ImVec2 perpendicular(-direction.y, direction.x);
                                    const ImVec2 tip = end_pixels;
                                    const ImVec2 left = ImVec2(
                                        tip.x - direction.x * arrow_size + perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size + perpendicular.y * (arrow_size * 0.5f));
                                    const ImVec2 right = ImVec2(
                                        tip.x - direction.x * arrow_size - perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size - perpendicular.y * (arrow_size * 0.5f));
                                    draw_list->AddTriangleFilled(tip, left, right, ImGui::GetColorU32(color));
                                }

                                ImGui::PopID();
                            }
                            ImGui::PopID();
                        }
                    }

                    if (draw_list != nullptr)
                    {
                        ImPlot::PopPlotClipRect();
                    }

                    ImPlot::EndPlot();
                }

                if (ImPlot::BeginPlot("Z^{+} / ~"))
                {
                    ImPlot::SetupAxes("n", "Index");
                    ImPlot::SetupAxesLimits(0.5, static_cast<double>(kClassCount) + 0.5, 0.5, static_cast<double>(kPairGridMax) + 0.5, ImPlotCond_Once);

                    ImDrawList* draw_list = ImPlot::GetPlotDrawList();
                    if (draw_list != nullptr)
                    {
                        ImPlot::PushPlotClipRect();
                    }

                    char label_buffer[32];
                    for (int class_idx = 0; class_idx < kClassCount; ++class_idx)
                    {
                        const auto& storage = equivalence_classes[class_idx];
                        if (storage.z_x.empty())
                        {
                            continue;
                        }

                        const ImVec4 color = class_colors[class_idx];
                        const ImVec4 fill_color(color.x, color.y, color.z, 0.25f);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 9.0f, fill_color, 1.5f, color);
                        std::snprintf(label_buffer, sizeof(label_buffer), "[%d]", class_idx + 1);
                        ImPlot::PlotScatter(label_buffer, storage.z_x.data(), storage.z_y.data(), class_sizes[class_idx]);

                        if (draw_list != nullptr && class_sizes[class_idx] >= 2)
                        {
                            ImGui::PushID(class_idx);
                            for (int member_idx = 1; member_idx < class_sizes[class_idx]; ++member_idx)
                            {
                                ImGui::PushID(member_idx);

                                const ImVec2 start_pixels = ImPlot::PlotToPixels(
                                    storage.z_x[member_idx - 1],
                                    storage.z_y[member_idx - 1]);
                                const ImVec2 end_pixels = ImPlot::PlotToPixels(
                                    storage.z_x[member_idx],
                                    storage.z_y[member_idx]);

                                draw_list->AddLine(start_pixels, end_pixels, ImGui::GetColorU32(color), 1.5f);

                                ImVec2 direction(end_pixels.x - start_pixels.x, end_pixels.y - start_pixels.y);
                                const float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
                                if (length > 1e-3f)
                                {
                                    direction.x /= length;
                                    direction.y /= length;
                                    const float arrow_size = 8.0f;
                                    const ImVec2 perpendicular(-direction.y, direction.x);
                                    const ImVec2 tip = end_pixels;
                                    const ImVec2 left = ImVec2(
                                        tip.x - direction.x * arrow_size + perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size + perpendicular.y * (arrow_size * 0.5f));
                                    const ImVec2 right = ImVec2(
                                        tip.x - direction.x * arrow_size - perpendicular.x * (arrow_size * 0.5f),
                                        tip.y - direction.y * arrow_size - perpendicular.y * (arrow_size * 0.5f));
                                    draw_list->AddTriangleFilled(tip, left, right, ImGui::GetColorU32(color));
                                }

                                ImGui::PopID();
                            }
                            ImGui::PopID();
                        }
                    }

                    if (draw_list != nullptr)
                    {
                        ImPlot::PopPlotClipRect();
                    }

                    ImPlot::EndPlot();
                }

                ImPlot::EndSubplots();
            }

            ImGui::TreePop();
            ImGui::Spacing();
        }
    }



    if (ImGui::Button("Close Me"))
        state.show_topology_window = false;
    ImGui::End();
}
} // namespace math_gui::plugins 