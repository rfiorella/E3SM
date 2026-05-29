# add_water_tracer.cmake
# CMake function for registering water tracers in EAMxx
#
# Usage:
#   add_water_tracer(
#     NAME <tracer_name>
#     LONGNAME <descriptive_name>
#     [KIND <bulk|evaporation|isotope|...>]
#     [RATIO_STANDARD <standard_name>]
#     [PARTITION_GROUP_ID <group_id>]
#   )
#
# Effect:
#   - Registers tracer metadata in scream::water_tracers::registry
#   - Increments SCREAM_NUM_TRACERS at configure time
#   - Tracer applies to ALL water reservoirs automatically
#   - Units are implicitly kg/kg (water mass mixing ratio)

function(add_water_tracer)
  set(options "")
  set(oneValueArgs NAME LONGNAME KIND RATIO_STANDARD PARTITION_GROUP_ID)
  set(multiValueArgs "")
  cmake_parse_arguments(TRACER "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Validate required arguments
  if(NOT TRACER_NAME)
    message(FATAL_ERROR "add_water_tracer requires NAME argument")
  endif()

  if(NOT TRACER_LONGNAME)
    message(FATAL_ERROR "add_water_tracer requires LONGNAME argument")
  endif()

  # Default KIND to evaporation for Group 1 test tracers
  if(NOT TRACER_KIND)
    set(TRACER_KIND "evaporation")
  endif()

  # Validate KIND
  set(VALID_KINDS bulk evaporation isotope tagged_partition tagged_diagnostic auxiliary)
  if(NOT TRACER_KIND IN_LIST VALID_KINDS)
    message(WARNING "Unknown tracer KIND: ${TRACER_KIND}. Valid values: ${VALID_KINDS}")
  endif()

  # Increment SCREAM_NUM_TRACERS
  if(NOT DEFINED SCREAM_NUM_TRACERS)
    set(SCREAM_NUM_TRACERS 1 CACHE STRING "Number of water tracer slots including bulk")
  endif()

  math(EXPR SCREAM_NUM_TRACERS_NEW "${SCREAM_NUM_TRACERS} + 1")
  set(SCREAM_NUM_TRACERS ${SCREAM_NUM_TRACERS_NEW} CACHE STRING "Number of water tracer slots including bulk" FORCE)

  # Append to global tracer lists
  get_property(TRACER_NAMES GLOBAL PROPERTY SCREAM_TRACER_NAMES)
  get_property(TRACER_LONGNAMES GLOBAL PROPERTY SCREAM_TRACER_LONGNAMES)
  get_property(TRACER_KINDS GLOBAL PROPERTY SCREAM_TRACER_KINDS)

  list(APPEND TRACER_NAMES "${TRACER_NAME}")
  list(APPEND TRACER_LONGNAMES "${TRACER_LONGNAME}")
  list(APPEND TRACER_KINDS "${TRACER_KIND}")

  set_property(GLOBAL PROPERTY SCREAM_TRACER_NAMES "${TRACER_NAMES}")
  set_property(GLOBAL PROPERTY SCREAM_TRACER_LONGNAMES "${TRACER_LONGNAMES}")
  set_property(GLOBAL PROPERTY SCREAM_TRACER_KINDS "${TRACER_KINDS}")

  # Store optional metadata
  if(TRACER_RATIO_STANDARD)
    set_property(GLOBAL APPEND PROPERTY SCREAM_TRACER_RATIO_STANDARDS "${TRACER_RATIO_STANDARD}")
  endif()

  if(TRACER_PARTITION_GROUP_ID)
    set_property(GLOBAL APPEND PROPERTY SCREAM_TRACER_PARTITION_GROUPS "${TRACER_PARTITION_GROUP_ID}")
  endif()

  message(STATUS "Registered water tracer: ${TRACER_NAME} (${TRACER_LONGNAME}) KIND=${TRACER_KIND}")
  message(STATUS "  SCREAM_NUM_TRACERS is now ${SCREAM_NUM_TRACERS_NEW}")
endfunction()

# Initialize tracer lists if not already set
get_property(TRACER_NAMES_INIT GLOBAL PROPERTY SCREAM_TRACER_NAMES SET)
if(NOT TRACER_NAMES_INIT)
  # Bulk water (slot 0) is always present
  set_property(GLOBAL PROPERTY SCREAM_TRACER_NAMES "bulk_H2O")
  set_property(GLOBAL PROPERTY SCREAM_TRACER_LONGNAMES "Bulk Water (Total)")
  set_property(GLOBAL PROPERTY SCREAM_TRACER_KINDS "bulk")

  # Ensure SCREAM_NUM_TRACERS is set
  if(NOT DEFINED SCREAM_NUM_TRACERS)
    set(SCREAM_NUM_TRACERS 1 CACHE STRING "Number of water tracer slots including bulk")
  endif()

  message(STATUS "Initialized water tracer system: SCREAM_NUM_TRACERS=${SCREAM_NUM_TRACERS}")
endif()
