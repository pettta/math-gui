#pragma once

#include <functional>
#include <string>

#include "imgui.h"
#include "implot.h"
#include "implot3d.h" 

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

class ImGuiBackend {
public:
    virtual ~ImGuiBackend() = default;

    virtual void initializeBackend() = 0;
    virtual void newFrame() = 0;
    virtual void renderDrawData(ImDrawData* drawData) = 0;
    virtual void shutdownBackend() = 0;
};

class ImGuiRenderer {
public:
    explicit ImGuiRenderer(ImGuiBackend& backend);
    ~ImGuiRenderer();

    struct FrameState {
        bool show_demo_window = true;
        bool show_linear_algebra_window = false; 
        bool show_probability_window = false; 
        bool show_topology_window = false;
        bool show_analysis_window = false;
        float clear_color[4] = {0.45f, 0.55f, 0.60f, 1.00f};

        // GIF recorder: desired state toggled by the header button / Cmd+Shift+G.
        // The macOS engine reacts to transitions and writes gif_status back for display.
        bool gif_recording = false;
        std::string gif_status;
    };

    void initialize();
    ImGuiIO& io();
    FrameState& frameState();
    const FrameState& frameState() const;

    void beginFrame();
    void businessLogic(FrameState& state);
    void endFrame();
    void renderFrame(const std::function<void(FrameState&)>& drawUI);

    void shutdown();

private:
    ImGuiBackend& backend_;
    bool initialized_{false};
    FrameState initialState_{};
};
