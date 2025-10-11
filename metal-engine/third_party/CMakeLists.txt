set(IMGUI_DIR "${PROJECT_SOURCE_DIR}/third_party/imgui")

set(IMGUI_SOURCES
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_impl_sdl2.cpp
    ${IMGUI_DIR}/imgui_impl_metal.mm
)

add_library(imgui STATIC ${IMGUI_SOURCES})

set_target_properties(imgui PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(imgui
    PUBLIC
        ${IMGUI_DIR}
)

target_include_directories(imgui
    PRIVATE
        ${SDL2_INCLUDE_DIRS}
)

target_compile_options(imgui PRIVATE
    -Wall
    -Wformat
    $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-weak>
    $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>
)
