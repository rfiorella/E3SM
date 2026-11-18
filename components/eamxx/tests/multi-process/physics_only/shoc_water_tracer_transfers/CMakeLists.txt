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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/tracers.yaml COPYONLY)

# Two runs of the same tree:
#  - no evaporation: surf_evap = 0 over an ocean that evaporates (the land part
#    dews by the same amount), so the tracers get no surface flux
#  - evap: surf_evap > 0 and half ocean, half land and sea ice. The bulk water
#    and the tracers then gain water through the surface, which SHOC applies in
#    the same solve. ocn_ssq is chosen so that columns span all four cases of
#    surface_tracer_flux: ocean evaporating or dewing, land evaporating or dewing.
# check_water_tracer_transfers.py must pass on both, with the same tolerance.
foreach (POSTFIX no_evap evap)
  if (POSTFIX STREQUAL no_evap)
    set (SURF_EVAP 0.0)
    set (OCNFRAC   1.0)
    set (OCN_SSQ   0.02)
    set (OUT_PREFIX shoc_water_tracer_transfers_output)
  else()
    set (SURF_EVAP 1.5e-6)
    set (OCNFRAC   0.5)
    set (OCN_SSQ   0.01)
    set (OUT_PREFIX shoc_water_tracer_transfers_evap_output)
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
endforeach()

CreateUnitTestFromExec(shoc_water_tracer_transfers shoc_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_no_evap.yaml"
  LABELS shoc water_tracers physics
  FIXTURES_SETUP shoc_water_tracer_transfers)
CreateUnitTestFromExec(shoc_water_tracer_transfers_evap shoc_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_evap.yaml"
  LABELS shoc water_tracers physics
  FIXTURES_SETUP shoc_water_tracer_transfers_evap)

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

# With evaporation on: the tracers' surface flux must close against bulk water
add_test (NAME shoc_water_tracer_transfers_evap_check
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_tracer_transfers.py
    -f shoc_water_tracer_transfers_evap_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --phases-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --exact-qr-qi
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_water_tracer_transfers_evap_check PROPERTIES
  LABELS "shoc;water_tracers;physics"
  FIXTURES_REQUIRED shoc_water_tracer_transfers_evap)
