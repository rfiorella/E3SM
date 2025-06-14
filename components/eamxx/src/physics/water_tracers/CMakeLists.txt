add_library(water_tracers eamxx_water_tracers.cpp)
target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER)
target_link_libraries(water_tracers physics_share scream_share)

if (TARGET eamxx_physics)
  # add to eamxx-physics
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()

# add tests
#add_subdirectory(tests)
