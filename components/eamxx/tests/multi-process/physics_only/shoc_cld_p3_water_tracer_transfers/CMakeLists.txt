include (ScreamUtils)

# water_tracers + water_tracer_transfers with the P3 replay: tracers with a
# uniform ratio R to bulk water must keep it through P3's clips, process
# rates, sedimentation and homogeneous freezing, in every phase and in the
# surface precipitation (see check_water_tracer_transfers.py). Two runs:
#  - p3 only:  [water_tracers, p3, water_tracer_transfers], no SHOC replay
#  - mac:      [water_tracers, mac], mac = [shoc, cld_fraction, p3,
#              water_tracer_transfers] subcycled, both replays
# and the driver link (AtmosphereDriver::setup_water_tracer_transfer_links):
#  - link:     the mac run with no water_tracers_file for the transfers; the
#              link must copy it. Same checks, same numbers.
#  - neg_*:    mis-configured trees that the link must reject during init
#              (incl. water_tracers after shoc, and SHOC's flux fixer on: SHOC
#              applies the tracer surface fluxes).
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

set (WTT_NEG_CASES order group nowtrc subcycle wtrc_after_shoc fixer)
# These two need a process that updates bulk water, which no linked physics does:
# wtt_link_modifier_test.cpp adds a test-only one (wtt_dummy_modifier) to the driver test
set (WTT_NEG_MOD_CASES modifier_qc modifier_group)
set (WTT_NEG_inputs)
foreach (c ${WTT_NEG_CASES} ${WTT_NEG_MOD_CASES})
  list (APPEND WTT_NEG_inputs input_neg_${c}.yaml)
endforeach()

# The AD test driver, with the dummy modifier process registered
CreateUnitTestExec(wtt_link_modifier_test
  SOURCES wtt_link_modifier_test.cpp
  LIBS scream_control eamxx_io eamxx_diagnostics shoc cld_fraction p3 water_tracers)
target_include_directories(wtt_link_modifier_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

foreach (f input_p3.yaml output_p3.yaml ${WTT_NEG_inputs})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f})
endforeach()

# The mac and link runs come in two variants (tracer surface fluxes, applied by SHOC):
#  - no evaporation: surf_evap = 0 over an ocean that evaporates (the land part
#    dews by the same amount), so the tracers get no surface flux
#  - evap: surf_evap > 0, half ocean and half land and sea ice. ocn_ssq is chosen
#    so that columns span all four cases of surface_tracer_flux: ocean evaporating
#    or dewing, land evaporating or dewing. The bulk water and the tracers then
#    gain water through the surface, in the same SHOC solve.
# (The p3-only run has no SHOC, so it has no surface fluxes to vary.)
foreach (v no_evap evap)
  if (v STREQUAL no_evap)
    set (SURF_EVAP 0.0)
    set (OCNFRAC   1.0)
    set (OCN_SSQ   0.02)
    set (SUFFIX "")
  else()
    set (SURF_EVAP 1.5e-6)
    set (OCNFRAC   0.5)
    set (OCN_SSQ   0.01)
    set (SUFFIX "_evap")
  endif()
  foreach (base input output input_link output_link)
    set (OUT_YAML output${SUFFIX}.yaml)
    set (OUT_PREFIX shoc_cld_p3_water_tracer_transfers${SUFFIX}_output)
    if (base MATCHES "link")
      set (OUT_YAML output_link${SUFFIX}.yaml)
      set (OUT_PREFIX shoc_cld_p3_wtt_link${SUFFIX}_output)
    endif()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${base}.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/${base}${SUFFIX}.yaml)
  endforeach()
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

# mac with evaporation on: R = 1 tracers must still close against bulk water
# (and against the surface precipitation), and R = 0.5 tracers must stay half
CreateUnitTestFromExec(shoc_cld_p3_water_tracer_transfers_evap shoc_cld_p3_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_evap.yaml"
  LABELS shoc cld p3 water_tracers physics
  FIXTURES_SETUP shoc_cld_p3_water_tracer_transfers_evap)
add_test (NAME shoc_cld_p3_water_tracer_transfers_evap_check
  COMMAND ${WTT_CHECK}
    -f shoc_cld_p3_water_tracer_transfers_evap_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    ${WTT_HEADERS}
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --min-p3-transfers 8
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_cld_p3_water_tracer_transfers_evap_check PROPERTIES
  LABELS "shoc;cld;p3;water_tracers;physics"
  FIXTURES_REQUIRED shoc_cld_p3_water_tracer_transfers_evap)

# link: nothing set by hand
CreateUnitTestFromExec(shoc_cld_p3_wtt_link shoc_cld_p3_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_link.yaml"
  LABELS shoc cld p3 water_tracers physics
  FIXTURES_SETUP shoc_cld_p3_wtt_link)
add_test (NAME shoc_cld_p3_wtt_link_check
  COMMAND ${WTT_CHECK}
    -f shoc_cld_p3_wtt_link_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    ${WTT_HEADERS}
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --min-p3-transfers 8
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_cld_p3_wtt_link_check PROPERTIES
  LABELS "shoc;cld;p3;water_tracers;physics"
  FIXTURES_REQUIRED shoc_cld_p3_wtt_link)

# link with evaporation on
CreateUnitTestFromExec(shoc_cld_p3_wtt_link_evap shoc_cld_p3_water_tracer_transfers
  EXE_ARGS "--args -ifile=input_link_evap.yaml"
  LABELS shoc cld p3 water_tracers physics
  FIXTURES_SETUP shoc_cld_p3_wtt_link_evap)
add_test (NAME shoc_cld_p3_wtt_link_evap_check
  COMMAND ${WTT_CHECK}
    -f shoc_cld_p3_wtt_link_evap_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    ${WTT_HEADERS}
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --min-p3-transfers 8
    --rtol ${WTT_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_cld_p3_wtt_link_evap_check PROPERTIES
  LABELS "shoc;cld;p3;water_tracers;physics"
  FIXTURES_REQUIRED shoc_cld_p3_wtt_link_evap)

# negative tests: the run must stop in the link, with the given message. With
# PASS_REGULAR_EXPRESSION ctest ignores the exit code, so the regex is the check.
# The error printer wraps lines, so the regexes accept a newline between words
set (WTT_NEG_REGEX_order    "process 'p3' must come BEFORE[ \n]+'water_tracer_transfers'")
set (WTT_NEG_REGEX_group    "must be in the same atm process[ \n]+group")
set (WTT_NEG_REGEX_nowtrc   "requires process 'water_tracers',[ \n]+which is not in the atm process tree")
set (WTT_NEG_REGEX_subcycle "must not be sub-cycled[ \n]+individually")
set (WTT_NEG_REGEX_wtrc_after_shoc "Process 'water_tracers' must come BEFORE[ \n]+'shoc'")
set (WTT_NEG_REGEX_fixer  "'shoc::check_flux_state_consistency'[ \n]+is[ \n]+true,[ \n]+but[ \n]+SHOC[ \n]+applies[ \n]+the[ \n]+water[ \n]+tracer[ \n]+surface[ \n]+fluxes")
set (WTT_NEG_REGEX_modifier_qc    "'wtt_dummy_modifier' in group 'mac' updates bulk water field[ \n]+'qc' between SHOC/P3")
set (WTT_NEG_REGEX_modifier_group "'wtt_dummy_modifier' in group 'mac' updates the group[ \n]+'tracers' between SHOC/P3")
foreach (c ${WTT_NEG_CASES} ${WTT_NEG_MOD_CASES})
  set (neg_exe shoc_cld_p3_water_tracer_transfers)
  if (c IN_LIST WTT_NEG_MOD_CASES)
    set (neg_exe wtt_link_modifier_test)
  endif()
  CreateUnitTestFromExec(wtt_link_neg_${c} ${neg_exe}
    EXE_ARGS "--args -ifile=input_neg_${c}.yaml"
    LABELS water_tracers physics)
  set_tests_properties (wtt_link_neg_${c} PROPERTIES
    PASS_REGULAR_EXPRESSION "${WTT_NEG_REGEX_${c}}")
endforeach()
