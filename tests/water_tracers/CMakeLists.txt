# Water tracer tests for EAMxx
#
# This directory contains unit tests and validation tests for water tracer functionality

# Test 1: Header compilation test (standalone, no EAMxx build required)
add_test(
  NAME water_tracer_headers_compile
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test_build_water_tracers.sh
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

# Unit tests for qv tracer access in each process (PR 2d)
# These tests validate that each process correctly accesses qv slot-0 from tracer-aware fields

add_executable(test_qv_rrtmgp test_qv_rrtmgp.cpp)
target_link_libraries(test_qv_rrtmgp scream_share Catch2::Catch2 Kokkos::kokkos)
add_test(NAME test_qv_rrtmgp COMMAND test_qv_rrtmgp)

add_executable(test_qv_zm test_qv_zm.cpp)
target_link_libraries(test_qv_zm scream_share Catch2::Catch2 Kokkos::kokkos)
add_test(NAME test_qv_zm COMMAND test_qv_zm)

add_executable(test_qv_homme test_qv_homme.cpp)
target_link_libraries(test_qv_homme scream_share Catch2::Catch2 Kokkos::kokkos)
add_test(NAME test_qv_homme COMMAND test_qv_homme)

add_executable(test_qv_surface test_qv_surface.cpp)
target_link_libraries(test_qv_surface scream_share Catch2::Catch2 Kokkos::kokkos)
add_test(NAME test_qv_surface COMMAND test_qv_surface)

# Future tests will be added here:
# - CMake add_water_tracer function test
# - Prototype qv_extension_test
# - Tracer ratio validation test (PR 5)
