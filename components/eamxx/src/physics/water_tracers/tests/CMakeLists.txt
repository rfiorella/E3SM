include(ScreamUtils)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Water tracer subview accessor tests
# These test the fundamental subview extraction infrastructure
# Tests run in all configurations (OFF and ON)
CreateUnitTest(water_tracer_subview_accessor "water_tracer_subview_accessor.cpp"
  THREADS 1
  LABELS "wtrc;water_tracers;unit"
)

# Water mass passive copy test (N=2)
# Only build when SCREAM_TRACE_WATER=ON and WTRC_MAX_CNST >= 2
if (SCREAM_TRACE_WATER AND WTRC_MAX_CNST GREATER_EQUAL 2)
  CreateUnitTest(water_mass_passive_copy_n2 "water_mass_passive_copy_n2.cpp"
    LIBS p3
    THREADS 1
    LABELS "wtrc;water_tracers;p3;physics"
  )
endif()

# Water tracer registry test (N=4)
# Only build when SCREAM_TRACE_WATER=ON and WTRC_MAX_CNST == 4
if (SCREAM_TRACE_WATER AND WTRC_MAX_CNST EQUAL 4)
  CreateUnitTest(water_tracer_registry_n4 "water_tracer_registry_n4_test.cpp"
    LIBS water_tracers ekat
    THREADS 1
    LABELS "wtrc;water_tracers;registry;unit"
  )
endif()

# Water tracer registry BFB test (bulk-only vs default)
# Only build when SCREAM_TRACE_WATER=ON and WTRC_MAX_CNST == 1
if (SCREAM_TRACE_WATER AND WTRC_MAX_CNST EQUAL 1)
  CreateUnitTest(water_tracer_registry_bulk_vs_default_bfb "water_tracer_registry_bfb_test.cpp"
    LIBS water_tracers ekat
    THREADS 1
    LABELS "wtrc;water_tracers;registry;bfb;unit"
  )
endif()

# Water isotope fractionation tests - temporarily disabled
# The water_isotopes.hpp header requires additional dependencies
# that are not yet configured. These tests will be re-enabled
# when full isotope support is implemented.

# set(WISO_TESTS_SRCS
#     water_isotopes_alphaeq_vi.cpp)

# # All tests should understand the same baseline args
# if (SCREAM_ENABLE_BASELINE_TESTS)
#   if (SCREAM_ONLY_GENERATE_BASELINES)
#     set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
#     # We don't want to do thread spreads when generating. That
#     # could cause race conditions in the file system.
#     set(WISO_THREADS "${SCREAM_TEST_MAX_THREADS}")
#   else()
#     set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
#     set(WISO_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
#   endif()
# else()
#   set(BASELINE_FILE_ARG "-n") # no baselines
#   set(WISO_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
# endif()

# CreateUnitTest(water_tracer_tests "${WISO_TESTS_SRCS}"
#   LIBS water_tracers
#   THREADS ${WISO_THREADS}
#   EXE_ARGS "--args ${BASELINE_FILE_ARG}"
#   LABELS "wtrc;wiso;physics"
#   )