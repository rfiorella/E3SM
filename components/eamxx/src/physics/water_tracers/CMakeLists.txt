# Water tracer support library
# Note: SCREAM_TRACE_WATER and SCREAM_NUM_WATER_TRACERS are configured
# in src/CMakeLists.txt and propagated via add_compile_definitions()

add_library(water_tracers eamxx_water_tracers.cpp)

# Compile definitions now set at physics level (physics/CMakeLists.txt)
# to ensure they propagate to all consumers including atm_process

target_link_libraries(water_tracers physics_share scream_share)

if (TARGET eamxx_physics)
  # add to eamxx-physics
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()

# add tests
add_subdirectory(tests)
