if (NOT SCREAM_LIB_ONLY)
  # Create test_support lib. Subfolders will add sources as we scan them
  add_library(eamxx_test_support)
  target_link_libraries (eamxx_test_support PUBLIC Catch2::Catch2)
endif()

# Note: Water tracer configuration (SCREAM_TRACE_WATER, SCREAM_NUM_WATER_TRACERS,
# WTRC_MAX_CNST) is set at components/eamxx/CMakeLists.txt level to ensure it
# propagates to all subdirectories including tests

add_subdirectory(share)

add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(control)
if (PROJECT_NAME STREQUAL "E3SM")
  add_subdirectory(mct_coupling)
endif()

option (EAMXX_ENABLE_PYSCREAM "Whether to enable interfaces to call eamxx from python" OFF)
if (EAMXX_ENABLE_PYSCREAM)
  add_subdirectory(python)
endif()
