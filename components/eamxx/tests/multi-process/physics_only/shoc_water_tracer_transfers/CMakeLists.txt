include (ScreamUtils)

# water_tracers + SHOC + water_tracer_transfers (SHOC replay only): tracers with
# a uniform ratio R to bulk water must keep it through SHOC's transport and
# condensation, in every phase (see check_water_tracer_transfers.py), and the
# rebuilt stage temperatures (water_transfer_stage_T; for SHOC a bookkeeping
# check only, and P3's stages must hold the fill value).
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

# The generic-tracer runs below: tracers.yaml, no alpha override, no kinetic
# isotope fractionation
set (WTT_TRACERS_FILE tracers.yaml)
set (WTT_ALPHA_OVERRIDE 0.0)
set (WTT_EXTRA_IC "")
set (WTT_KINETIC false)
set (WTT_KINETIC_ICE_S empirical)
set (WTT_KINETIC_OUTPUT "")

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

# Tolerance of the roundoff checks, relative to the total water at a point,
# and of the stage temperature checks, relative to T_mid
if (SCREAM_DOUBLE_PRECISION)
  set (WTT_RTOL 1e-12)
  set (WTT_T_RTOL 1e-13)
else()
  set (WTT_RTOL 1e-5)
  set (WTT_T_RTOL 1e-5)
endif()
add_test (NAME shoc_water_tracer_transfers_check
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_tracer_transfers.py
    -f shoc_water_tracer_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    --phases-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp
    --tracer wtrc_t=1 wtrc_half=0.5
    --exact-qr-qi
    --rtol ${WTT_RTOL}
    --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
    --T-rtol ${WTT_T_RTOL}
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
    --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
    --T-rtol ${WTT_T_RTOL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (shoc_water_tracer_transfers_evap_check PROPERTIES
  LABELS "shoc;water_tracers;physics"
  FIXTURES_REQUIRED shoc_water_tracer_transfers_evap)

# Fractionating isotopes (HDO, H218O, with the equilibrium alphas), next to a
# generic R = 1 tracer, in the evap configuration (surface fluxes on), see
# check_water_isotope_transfers.py:
#  - iso:          the 4-step evap run above, with tracers_iso.yaml. The
#                  tracers must be conserved (column
#                  + surface flux) and stay >= 0 (to roundoff), condensing
#                  cloud must be enriched relative to the vapor, and where SHOC
#                  made new cloud from vapor R_qc/R_qv must be the Rayleigh
#                  result for alpha_l at the stage temperature.
#  - iso_newcloud: one SHOC call (dt = 300 s) from a state with no cloud
#                  liquid, so that all of its condensation is new cloud: the
#                  last check at every condensing point.
# and with kinetic_fractionation on (the replay's kinetic factors):
#  - iso_kinetic:  the iso run: the same checks (net condensation stays at
#                  equilibrium, so the Rayleigh check of new cloud still holds)
#  - iso_kinetic_onecall: one SHOC call (dt = 300 s) from the initial state,
#                  where every phase has R = 1 and the IC's cloud exists. Where
#                  SHOC's net transfer evaporates cloud, the cloud left must
#                  have R_qc = f^(1/akel - 1) (f the fraction of the cloud
#                  left), with akel computed in the check from alpha_l and the
#                  relative humidity over liquid of SHOC's vapor after
#                  transport, and (in double precision) differ from the
#                  equilibrium result at a point at least.
if (EAMXX_ENABLE_WATER_ISOTOPES)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_iso.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_iso.yaml COPYONLY)
  # Tolerance of R_qc/R_qv against the Rayleigh result, and the new cloud
  # points each run must have (the 4-step run has few: cloud that SHOC
  # transports is rarely exactly absent where it condenses)
  # The one-call kinetic run starts from the IC's cloud, which leaves traces
  # (< 1e-8 of the condensation) at its new-cloud points, with R = 1: a looser
  # tolerance there. Its check of SHOC's net evaporation has its own, and must
  # resolve the kinetic factor at a point at least (not in single precision:
  # SHOC evaporates cloud close to saturation, where akel is within 1e-4 of
  # alpha_l).
  if (SCREAM_DOUBLE_PRECISION)
    set (WTT_ALPHA_RTOL 1e-10)
    set (WTT_ALPHA_RTOL_ONECALL 1e-8)
    set (WTT_EVAP_RTOL 1e-12)
    set (WTT_MIN_EVAP_KINETIC 1)
  else()
    set (WTT_ALPHA_RTOL 1e-3)
    set (WTT_ALPHA_RTOL_ONECALL 1e-3)
    set (WTT_EVAP_RTOL 1e-4)
    set (WTT_MIN_EVAP_KINETIC 0)
  endif()
  set (WTT_MIN_NEW_CLOUD_iso 1)
  set (WTT_MIN_NEW_CLOUD_iso_newcloud 100)
  set (WTT_MIN_NEW_CLOUD_iso_kinetic 1)
  set (WTT_MIN_NEW_CLOUD_iso_kinetic_onecall 0)
  set (WTT_TRACERS_FILE tracers_iso.yaml)
  set (WTT_ALPHA_OVERRIDE 0.0)
  set (SURF_EVAP 1.5e-6)
  set (OCNFRAC   0.5)
  set (OCN_SSQ   0.01)
  foreach (POSTFIX iso iso_newcloud iso_kinetic iso_kinetic_onecall)
    set (OUT_PREFIX shoc_water_${POSTFIX}_transfers_output)
    set (NUM_STEPS 4)
    set (ATM_TIME_STEP 1800)
    set (MAC_SUBCYCLES 6)
    set (WTT_EXTRA_IC "")
    set (WTT_KINETIC false)
    set (WTT_KINETIC_OUTPUT "")
    set (kinetic_args)
    set (alpha_rtol ${WTT_ALPHA_RTOL})
    if (POSTFIX MATCHES "newcloud|onecall")
      set (NUM_STEPS 1)
      set (ATM_TIME_STEP 300)
      set (MAC_SUBCYCLES 1)
    endif()
    if (POSTFIX STREQUAL iso_newcloud)
      set (WTT_EXTRA_IC "qc: 0.0")
    endif()
    if (POSTFIX MATCHES "kinetic")
      set (WTT_KINETIC true)
      # (the saturation inputs, which only the kinetic runs request)
      set (WTT_KINETIC_OUTPUT "- p_dry_mid\n  - pseudo_density_dry")
      set (kinetic_args --kinetic)
      if (POSTFIX MATCHES "onecall")
        list (APPEND kinetic_args --shoc-evap-rayleigh --evap-rtol ${WTT_EVAP_RTOL}
                                  --min-evap-kinetic ${WTT_MIN_EVAP_KINETIC})
        set (alpha_rtol ${WTT_ALPHA_RTOL_ONECALL})
      endif()
    endif()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_iso.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
    CreateUnitTestFromExec(shoc_water_${POSTFIX}_transfers shoc_water_tracer_transfers
      EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
      LABELS shoc water_tracers water_isotopes physics
      FIXTURES_SETUP shoc_water_${POSTFIX}_transfers)
    add_test (NAME shoc_water_${POSTFIX}_transfers_check
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_isotope_transfers.py
        -f shoc_water_${POSTFIX}_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
        --tracer wtrc_hdo=HDO wtrc_o18=H218O --reference wtrc_t
        --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
        --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
        --dt ${ATM_TIME_STEP}
        --rtol-budget ${WTT_RTOL} --neg-rtol ${WTT_RTOL}
        --alpha-check --alpha-rtol ${alpha_rtol}
        --min-new-cloud ${WTT_MIN_NEW_CLOUD_${POSTFIX}}
        ${kinetic_args}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties (shoc_water_${POSTFIX}_transfers_check PROPERTIES
      LABELS "shoc;water_tracers;water_isotopes;physics"
      FIXTURES_REQUIRED shoc_water_${POSTFIX}_transfers)
  endforeach()
endif()
