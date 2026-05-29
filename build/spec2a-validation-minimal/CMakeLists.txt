cmake_minimum_required(VERSION 3.18)
project(spec2a_validation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define SCREAM_NUM_TRACERS
add_compile_definitions(SCREAM_NUM_TRACERS=1)

# Include directories for EAMxx headers
include_directories(
  ${CMAKE_SOURCE_DIR}/../../src/share/field
  ${CMAKE_SOURCE_DIR}/../../src/share/grid
  ${CMAKE_SOURCE_DIR}/../../src/share/util
)

# Test 1: Compile field_tag.hpp with TRACER
add_executable(test_field_tag_compiles test_field_tag.cpp)

# Test 2: Test explicit indexing pattern
add_executable(test_explicit_pattern test_explicit.cpp)
target_compile_definitions(test_explicit_pattern PRIVATE SCREAM_TRACER_ACCESS=EXPLICIT)

# Test 3: Test subview pattern
add_executable(test_subview_pattern test_subview.cpp)
target_compile_definitions(test_subview_pattern PRIVATE SCREAM_TRACER_ACCESS=SUBVIEW)
