# Water tracer support for EAMxx

# Set default tracer count (bulk water only)
if(NOT DEFINED SCREAM_NUM_TRACERS)
  set(SCREAM_NUM_TRACERS 1 CACHE STRING "Number of water tracer slots including bulk")
endif()

# Initialize tracer name and longname lists with bulk water
set(TRACER_NAME_LIST "\"bulk_H2O\"")
set(TRACER_LONGNAME_LIST "\"Bulk Water (Total)\"")

# Function to add a water tracer
# This is a simplified version - the full implementation is in cmake/add_water_tracer.cmake
function(add_water_tracer)
  set(options "")
  set(oneValueArgs NAME LONGNAME KIND RATIO_STANDARD PARTITION_GROUP_ID)
  set(multiValueArgs "")
  cmake_parse_arguments(TRACER "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Validate required arguments
  if(NOT TRACER_NAME OR NOT TRACER_LONGNAME)
    message(FATAL_ERROR "add_water_tracer requires NAME and LONGNAME")
  endif()

  # Default KIND to evaporation for Group 1 test tracers
  if(NOT TRACER_KIND)
    set(TRACER_KIND "evaporation")
  endif()

  # Increment SCREAM_NUM_TRACERS
  math(EXPR SCREAM_NUM_TRACERS "${SCREAM_NUM_TRACERS} + 1")
  set(SCREAM_NUM_TRACERS ${SCREAM_NUM_TRACERS} PARENT_SCOPE)

  # Append to tracer lists (propagate to parent scope)
  set(TRACER_NAME_LIST "${TRACER_NAME_LIST}, \"${TRACER_NAME}\"" PARENT_SCOPE)
  set(TRACER_LONGNAME_LIST "${TRACER_LONGNAME_LIST}, \"${TRACER_LONGNAME}\"" PARENT_SCOPE)

  message(STATUS "Registered water tracer: ${TRACER_NAME} (KIND=${TRACER_KIND})")
endfunction()

# Generate water_tracer_config.hpp from template
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/water_tracer_config.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/water_tracer_config.hpp
  @ONLY
)

# Add to include path
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(water_tracers eamxx_water_tracers.cpp)
target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER)
target_link_libraries(water_tracers physics_share scream_share)

if (TARGET eamxx_physics)
  # add to eamxx-physics
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()

# add tests
add_subdirectory(tests)

message(STATUS "Water tracer configuration: SCREAM_NUM_TRACERS=${SCREAM_NUM_TRACERS}")
