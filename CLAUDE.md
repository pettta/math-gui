# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Math GUI is a cross-platform C++17 desktop app for interactively visualizing math
(probability distributions, topology, linear algebra). The UI is built with Dear ImGui +
ImPlot/ImPlot3D on top of SDL2. It renders through **Metal on macOS** and **Vulkan on
Linux/Windows** — these are two separate executables sharing the same frontend code.

## Build & Run

The CMake project name is `vulkan_guide` (legacy), but the macOS executable is `metal_engine`
and the Linux/Windows executable is `vulkan_engine`.

```bash
# One-shot: configure, build (Release), and launch the right engine for this platform
python3 local_setup.py

# Manual equivalent
cmake -B build .
cmake --build build --config Release
cd bin && ./metal_engine        # macOS;  ./vulkan_engine on Linux

# First clone only — Eigen is a git submodule (Boost is fetched automatically by CMake)
git submodule update --init --recursive
```

- Binaries are emitted to `bin/`. **Run from inside `bin/`** — the font loader resolves
  `../frontend/utils/fonts/...` relative to the working directory.
- There is no test suite, linter config, or CI in this repo. CMake passes `-Wall -Wformat`.
- Windows requires editing the `PATH_PREFIX` in `local_setup.py` to match the local MSVC install
  (see comments there). Deps: MSVC `cl`, CMake>=3.21, Nmake, Windows 10/11 SDK, VulkanSDK>=1.4.321.
- macOS deps (via brew): `clang++`, `cmake>=3.21`, `make`, and SDL2 (`find_package(SDL2)`).

## Architecture

`CMakeLists.txt` branches on `APPLE`: the macOS path (ObjC++/Metal) and the else path
(Vulkan + shader compilation) are almost entirely independent. The shared, platform-agnostic
code lives in `frontend/`.

**Backend abstraction (the key seam).** `frontend/imgui-renderer.h` defines two classes:
- `ImGuiBackend` — pure virtual interface: `initializeBackend / newFrame / renderDrawData / shutdownBackend`.
- `ImGuiRenderer` — owns the ImGui/ImPlot/ImPlot3D contexts, font loading, and the per-frame
  loop. It holds an `ImGuiBackend&` and is completely platform-unaware.

Concrete backends implement `ImGuiBackend`:
- `frontend/metal_imgui_backend.mm` (`MetalImguiBackend`)
- `frontend/vulkan_imgui_backend.cpp` (`VulkanImguiBackend`, uses Vulkan dynamic rendering)

**Entry points wire backend → renderer → main loop:**
- macOS: `metal-engine/metal_engine.mm` (`main`) sets up SDL2+Metal, constructs
  `MetalImguiBackend`, hands it to `ImGuiRenderer`, runs the loop, calling
  `setFrameResources(...)` each frame with the current command buffer/encoder.
- Vulkan: `vulkan-engine/main.cpp` → `VulkanEngine` (`vk_engine.cpp`). `init_imgui()` builds the
  `VulkanImguiBackend` + `ImGuiRenderer`; the draw loop calls `setFrameResources(...)` and
  `updateSwapchainInfo(...)` on resize. The Vulkan engine (in `shared/` + `vulkan-engine/`) is a
  full vkguide-style renderer; ImGui is drawn into the swapchain image on top.

**UI logic.** `ImGuiRenderer::businessLogic(FrameState&)` in `imgui-renderer.cpp` is where the
actual app windows live. `FrameState` holds the toggle bools (demo / linear algebra / probability
/ topology windows) and clear color. Each enabled window dispatches to a plugin.

**Plugins** (`frontend/plugins/`) — each exposes a single `RenderXxxWindow(FrameState&)` in
namespace `math_gui::plugins`:
- `probability-plugin.cpp` (the main feature, ~1200 lines). Distributions are **data-driven**:
  `kDistributions` is a `std::vector<DistributionEntry>` where each `DistributionDefinition`
  carries its domain, parameters (with names/ranges/descriptions/integral flags), and factory
  lambdas (`PdfFactory`/`CdfFactory`/`PpfFactory`/`StatisticFactory`) that capture the current
  parameter vector. **To add a distribution, append an entry to `kDistributions`** — the UI
  (sliders, plots, statistics table) is generated from the definition. Backed by
  `boost::math::distributions`.
  - `frontend/utils/distributions/beta_binomial.cpp` is a vendored/extended Boost-style
    distribution header used by the plugin.
- `topology-plugin.cpp`, plus a stubbed Linear Algebra window inline in `imgui-renderer.cpp`.

## Third-party

`third_party/` vendors imgui, implot, implot3d, volk, fmt, etc. (checked in, not submodules).
Eigen is a git submodule. Boost (headers only, math + random) is downloaded at configure time via
`cmake/BoostExternal.cmake` (`ExternalProject_Add boost_ep`, exposed as `Boost::headers`).

On Vulkan builds, GLSL shaders in `shaders/*.{vert,frag,comp}` are compiled to `.spv` by
`glslangValidator` via the `Shaders` custom target.
