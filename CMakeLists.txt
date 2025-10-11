# CMakeList.txt : CMake project for vulkan_guide, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.20)

project("vulkan_guide")

if(APPLE)
  enable_language(OBJCXX)

  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  find_package(SDL2 REQUIRED)

  set(_OUTPUT_DIR "${PROJECT_SOURCE_DIR}/bin")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${_OUTPUT_DIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${_OUTPUT_DIR}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${_OUTPUT_DIR}")

  if(CMAKE_CONFIGURATION_TYPES)
    foreach(_CONFIG ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER "${_CONFIG}" _CONFIG_UPPER)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_CONFIG_UPPER} "${_OUTPUT_DIR}")
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_CONFIG_UPPER} "${_OUTPUT_DIR}")
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_CONFIG_UPPER} "${_OUTPUT_DIR}")
    endforeach()
  endif()

  add_subdirectory(metal-engine/third_party)

  add_executable(metal_engine
    metal-engine/metal_engine.mm
    frontend/imgui-renderer.cpp
    frontend/metal_imgui_backend.mm
  )

  target_link_libraries(metal_engine PRIVATE
    imgui
  )

  if(TARGET SDL2::SDL2)
    target_link_libraries(metal_engine PRIVATE SDL2::SDL2)
  else()
    target_include_directories(metal_engine PRIVATE ${SDL2_INCLUDE_DIRS})
    target_link_libraries(metal_engine PRIVATE ${SDL2_LIBRARIES})
  endif()

  target_link_libraries(metal_engine PRIVATE
    "-framework Metal"
    "-framework MetalKit"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-framework QuartzCore"
  )

  target_compile_options(metal_engine PRIVATE
    -Wall
    -Wformat
    $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-weak>
    $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>
  )

else()

  find_package(Vulkan REQUIRED)

  add_subdirectory(third_party)

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

  add_subdirectory(shared)

  add_subdirectory(vulkan-engine)
  find_program(GLSL_VALIDATOR glslangValidator HINTS /usr/bin /usr/local/bin $ENV{VULKAN_SDK}/Bin/ $ENV{VULKAN_SDK}/Bin32/)

  file(GLOB_RECURSE GLSL_SOURCE_FILES
    "${PROJECT_SOURCE_DIR}/shaders/*.frag"
    "${PROJECT_SOURCE_DIR}/shaders/*.vert"
    "${PROJECT_SOURCE_DIR}/shaders/*.comp"
  )

  foreach(GLSL ${GLSL_SOURCE_FILES})
    message(STATUS "BUILDING SHADER")
    get_filename_component(FILE_NAME ${GLSL} NAME)
    set(SPIRV "${PROJECT_SOURCE_DIR}/shaders/${FILE_NAME}.spv")
    message(STATUS ${GLSL})
    message(STATUS COMMAND ${GLSL_VALIDATOR} -V ${GLSL} -o ${SPIRV})
    add_custom_command(
    OUTPUT ${SPIRV}
    COMMAND ${GLSL_VALIDATOR} -V ${GLSL} -o ${SPIRV}
    DEPENDS ${GLSL})
    list(APPEND SPIRV_BINARY_FILES ${SPIRV})
  endforeach(GLSL)

  add_custom_target(
    Shaders 
    DEPENDS ${SPIRV_BINARY_FILES}
  )

endif()
