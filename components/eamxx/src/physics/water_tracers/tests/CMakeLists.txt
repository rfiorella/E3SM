include(ScreamUtils)

include_directories(${CMAKE_CURRENT_SOURCE_DATA}/..)

set(WISO_TESTS_SRCS
    water_isotopes_alphaeq_vi.cpp)


# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(WISO_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(WISO_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(WISO_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(water_tracer_tests "${WISO_TESTS_SRCS}"
  LIBS water_tracers
  THREADS ${WISO_THREADS}
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  LABELS "wtrc;wiso;physics"
  )