find_package(Vulkan REQUIRED)

add_library(vkbootstrap STATIC)
add_library(glm INTERFACE)
add_library(vma INTERFACE)

add_library(stb_image INTERFACE)

add_subdirectory(fastgltf)
add_subdirectory(fmt EXCLUDE_FROM_ALL)
add_subdirectory(SDL EXCLUDE_FROM_ALL)

target_sources(vkbootstrap PRIVATE 
  vkbootstrap/VkBootstrap.h
  vkbootstrap/VkBootstrap.cpp
  )

target_include_directories(vkbootstrap PUBLIC vkbootstrap)
target_link_libraries(vkbootstrap PUBLIC Vulkan::Vulkan $<$<BOOL:UNIX>:${CMAKE_DL_LIBS}>)
set_property(TARGET vkbootstrap PROPERTY CXX_STANDARD 20)

#both vma and glm and header only libs so we only need the include path
target_include_directories(vma INTERFACE vma)
target_include_directories(glm INTERFACE glm)

#add_library(sdl2 INTERFACE)
#target_include_directories(sdl2 INTERFACE $ENV{VULKAN_SDK}/Include/SDL2 )
#target_link_directories(sdl2 INTERFACE $ENV{VULKAN_SDK}/Lib )
#target_link_libraries(sdl2 INTERFACE SDL2 SDL2main)

add_library(imgui STATIC)
target_include_directories(imgui PUBLIC imgui)
target_sources(imgui PRIVATE 
  imgui/imgui.h
  imgui/imgui.cpp

  imgui/imgui_demo.cpp
  imgui/imgui_draw.cpp
  imgui/imgui_widgets.cpp
  imgui/imgui_tables.cpp
  imgui/imgui_impl_vulkan.cpp
  imgui/imgui_impl_sdl2.cpp
  )
set_property(TARGET imgui PROPERTY CXX_STANDARD 20)
target_link_libraries(imgui PUBLIC Vulkan::Vulkan SDL2::SDL2)

add_library(implot STATIC)
target_include_directories(implot PUBLIC implot)
target_sources(implot PRIVATE
  implot/implot.cpp
  implot/implot_items.cpp
  implot/implot_demo.cpp
)
target_link_libraries(implot PUBLIC imgui)

add_library(implot3d STATIC)
target_include_directories(implot3d PUBLIC implot3d)
target_sources(implot3d PRIVATE
  implot3d/implot3d.cpp
  implot3d/implot3d_items.cpp
  implot3d/implot3d_demo.cpp
  implot3d/implot3d_meshes.cpp
)
target_link_libraries(implot3d PUBLIC implot)

add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE eigen)

target_include_directories(stb_image INTERFACE stb_image)
