# Water tracer tests for EAMxx
#
# This directory contains unit tests and validation tests for water tracer functionality

# Test 1: Header compilation test (standalone, no EAMxx build required)
add_test(
  NAME water_tracer_headers_compile
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test_build_water_tracers.sh
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

# Future tests will be added here:
# - CMake add_water_tracer function test
# - Prototype qv_extension_test
# - Tracer ratio validation test (PR 5)
