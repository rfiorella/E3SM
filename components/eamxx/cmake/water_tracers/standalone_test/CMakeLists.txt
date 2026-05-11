# Standalone test for water tracer registry generation
# This minimal CMakeLists.txt tests the registry mechanism without full EAMxx dependencies

cmake_minimum_required(VERSION 3.18)
project(WaterTracerRegistryTest CXX)

set(CMAKE_CXX_STANDARD 17)

# Load the registry mechanism from parent CMakeLists.txt
# For standalone testing, we'll inline the relevant parts here

# Global tracer list
set(SCREAM_WATER_TRACER_LIST "" CACHE INTERNAL "List of configured tracers")

# Function to register a water tracer
function(add_water_tracer)
  cmake_parse_arguments(TRACER "TAG" "NAME;ISOTOPOLOGUE" "" ${ARGN})

  if (NOT DEFINED TRACER_NAME)
    message(FATAL_ERROR "add_water_tracer: NAME is required")
  endif()
  if (NOT DEFINED TRACER_ISOTOPOLOGUE)
    message(FATAL_ERROR "add_water_tracer: ISOTOPOLOGUE is required")
  endif()

  # Set TAG flag (default false)
  if (TRACER_TAG)
    set(is_tag "true")
  else()
    set(is_tag "false")
  endif()

  # Append to global list
  set(entry "${TRACER_NAME}|${TRACER_ISOTOPOLOGUE}|${is_tag}")
  set(SCREAM_WATER_TRACER_LIST "${SCREAM_WATER_TRACER_LIST};${entry}"
      CACHE INTERNAL "List of configured tracers")
endfunction()

# Load tracer configuration
if (NOT SCREAM_WATER_TRACERS_FILE)
  message(FATAL_ERROR "Must specify -DSCREAM_WATER_TRACERS_FILE=<path>")
endif()

if (NOT EXISTS "${SCREAM_WATER_TRACERS_FILE}")
  message(FATAL_ERROR "File not found: ${SCREAM_WATER_TRACERS_FILE}")
endif()

include("${SCREAM_WATER_TRACERS_FILE}")

# Remove empty first element if present (artifact of list initialization)
if (SCREAM_WATER_TRACER_LIST)
  list(REMOVE_AT SCREAM_WATER_TRACER_LIST 0)
endif()

# Validation
list(LENGTH SCREAM_WATER_TRACER_LIST tracer_count)
if (tracer_count EQUAL 0)
  message(FATAL_ERROR "Tracer list is empty")
endif()

message(STATUS "Raw tracer list: ${SCREAM_WATER_TRACER_LIST}")

# Validate first tracer is H2O
list(GET SCREAM_WATER_TRACER_LIST 0 first_tracer)
string(REPLACE "|" ";" first_parts "${first_tracer}")
list(GET first_parts 1 first_iso)
if (NOT first_iso STREQUAL "H2O" AND NOT first_iso STREQUAL "H216O")
  message(FATAL_ERROR "First tracer must be H2O or H216O, found: ${first_iso}")
endif()

# Validate known isotopologues and check duplicates
set(known_isotopologues "H2O;H216O;HDO;H218O;H217O;HTO")
set(seen_names "")
foreach(entry ${SCREAM_WATER_TRACER_LIST})
  string(REPLACE "|" ";" parts "${entry}")
  list(GET parts 0 name)
  list(GET parts 1 iso)

  if (NOT iso IN_LIST known_isotopologues)
    message(FATAL_ERROR "Unknown isotopologue: ${iso}")
  endif()

  if (name IN_LIST seen_names)
    message(FATAL_ERROR "Duplicate tracer name: ${name}")
  endif()
  list(APPEND seen_names "${name}")
endforeach()

set(WTRC_MAX_CNST ${tracer_count})
message(STATUS "Tracers configured: ${tracer_count}")

# Generate registry
set(SCREAM_WATER_TRACER_REGISTRY_ENTRIES "")
foreach(entry ${SCREAM_WATER_TRACER_LIST})
  string(REPLACE "|" ";" parts "${entry}")
  list(GET parts 0 name)
  list(GET parts 1 iso)
  list(GET parts 2 is_tag)
  set(SCREAM_WATER_TRACER_REGISTRY_ENTRIES
      "${SCREAM_WATER_TRACER_REGISTRY_ENTRIES}  {\"${name}\", catalog_index(\"${iso}\"), ${is_tag}},\n")
endforeach()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/../water_tracer_registry.gen.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/water_tracer_registry.gen.hpp"
  @ONLY
)

# Print generated file for verification
file(READ "${CMAKE_CURRENT_BINARY_DIR}/water_tracer_registry.gen.hpp" GENERATED_CONTENT)
message(STATUS "Generated registry:\n${GENERATED_CONTENT}")
