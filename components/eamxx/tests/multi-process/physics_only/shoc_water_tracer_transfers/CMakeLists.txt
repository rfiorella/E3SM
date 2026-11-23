include (ScreamUtils)

# water_tracers + SHOC + water_tracer_transfers (SHOC replay only): tracers with
# a uniform ratio R to bulk water must keep it through SHOC's transport and
# condensation, in every phase (see check_water_tracer_transfers.py).
CreateADUnitTestExec(shoc_water_tracer_transfers LIBS shoc water_tracers)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

set (RUN_T0 2021-10-12-45000)
set (NUM_STEPS 4)
set (ATM_TIME_STEP 1800)
# SHOC and the replay run in a subcycled group, as in CIME runs (SHOC dt = 300 s)
set (MAC_SUBCYCLES 6)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/tracers.yaml COPYONLY)

CreateUnitTestFromExec(shoc_water_tracer_transfers shoc_water_tracer_transfers
  EXE_ARGS "--args -ifile=input.yaml"
  LABELS shoc water_tracers physics
  FIXTURES_SETUP shoc_water_tracer_transfers)

# Tolerance of the roundoff checks, relative to the total water at a point
if (SCREAM_DOUBLE_PRECISION)
  set (WTT_RTOL 1e-12)
else()
  set (WTT_RTOL 1e-5)
endif()
add_test (NAME shoc_water_tracer_transfers_check
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_tracer_transfers.py
    -f shoc_water_tracer_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --phases-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --exact-qr-qi
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_water_tracer_transfers_check PROPERTIES
  LABELS "shoc;water_tracers;physics"
  FIXTURES_REQUIRED shoc_water_tracer_transfers)
