include (ScreamUtils)

# water_tracers + water_tracer_transfers with the P3 replay: tracers with a
# uniform ratio R to bulk water must keep it through P3's clips, process
# rates, sedimentation and homogeneous freezing, in every phase and in the
# surface precipitation (see check_water_tracer_transfers.py). Two runs:
#  - p3 only:  [water_tracers, p3, water_tracer_transfers], no SHOC replay
#  - mac:      [water_tracers, mac], mac = [shoc, cld_fraction, p3,
#              water_tracer_transfers] subcycled, both replays
CreateADUnitTestExec(shoc_cld_p3_water_tracer_transfers LIBS shoc cld_fraction p3 water_tracers)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

set (RUN_T0 2021-10-12-45000)
set (NUM_STEPS 4)
# p3 only: P3 runs once per step, with the dt it gets inside a mac group
set (P3_TIME_STEP 300)
# mac: SHOC, P3 and the replay run in a subcycled group, as in CIME runs (dt = 300 s)
set (ATM_TIME_STEP 1800)
set (MAC_SUBCYCLES 6)

foreach (f input_p3.yaml output_p3.yaml input.yaml output.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f})
endforeach()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/tracers.yaml COPYONLY)

# Tolerance of the roundoff checks, relative to the total water at a point
if (SCREAM_DOUBLE_PRECISION)
  set (WTT_RTOL 1e-12)
else()
  set (WTT_RTOL 1e-5)
endif()

set (WTT_CHECK ${CMAKE_CURRENT_SOURCE_DIR}/../shoc_water_tracer_transfers/check_water_tracer_transfers.py)
set (WTT_HEADERS
  --phases-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp
  --p3-header ${SCREAM_SRC_DIR}/physics/p3/p3_water_transfers.hpp)

# p3 only
CreateUnitTestFromExec(p3_water_tracer_transfers shoc_cld_p3_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_p3.yaml"
  LABELS p3 water_tracers physics
  FIXTURES_SETUP p3_water_tracer_transfers)
add_test (NAME p3_water_tracer_transfers_check
  COMMAND ${WTT_CHECK}
    -f p3_water_tracer_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    ${WTT_HEADERS}
    --tracer wtrc_t=1 wtrc_half=0.5
    --min-p3-transfers 8
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (p3_water_tracer_transfers_check PROPERTIES
  LABELS "p3;water_tracers;physics"
  FIXTURES_REQUIRED p3_water_tracer_transfers)

# mac
CreateUnitTestFromExec(shoc_cld_p3_water_tracer_transfers shoc_cld_p3_water_tracer_transfers
  EXE_ARGS "--args -ifile=input.yaml"
  LABELS shoc cld p3 water_tracers physics
  FIXTURES_SETUP shoc_cld_p3_water_tracer_transfers)
add_test (NAME shoc_cld_p3_water_tracer_transfers_check
  COMMAND ${WTT_CHECK}
    -f shoc_cld_p3_water_tracer_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    ${WTT_HEADERS}
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --min-p3-transfers 8
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_cld_p3_water_tracer_transfers_check PROPERTIES
  LABELS "shoc;cld;p3;water_tracers;physics"
  FIXTURES_REQUIRED shoc_cld_p3_water_tracer_transfers)
