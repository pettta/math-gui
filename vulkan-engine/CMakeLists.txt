
# Add source to this project's executable.
add_executable (vulkan_engine 
  main.cpp
  vk_engine.cpp
  vk_engine.h
)

set_property(TARGET vulkan_engine PROPERTY CXX_STANDARD 20)
set_property(TARGET vulkan_engine PROPERTY VS_DEBUGGER_WORKING_DIRECTORY $<TARGET_FILE_DIR:vulkan_engine>)

target_include_directories(vulkan_engine PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


target_link_libraries(vulkan_engine vkguide_shared vkbootstrap imgui)

add_custom_command(TARGET vulkan_engine POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:vulkan_engine> $<TARGET_FILE_DIR:vulkan_engine>
  COMMAND_EXPAND_LISTS
  )
