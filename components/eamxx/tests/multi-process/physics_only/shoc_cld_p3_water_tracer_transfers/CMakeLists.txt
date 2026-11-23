include (ScreamUtils)

# water_tracers + water_tracer_transfers with the P3 replay: tracers with a
# uniform ratio R to bulk water must keep it through P3's clips, process
# rates, sedimentation and homogeneous freezing, in every phase and in the
# surface precipitation (see check_water_tracer_transfers.py). The rebuilt
# stage temperatures (water_transfer_stage_T) are checked too. Two runs:
#  - p3 only:  [water_tracers, p3, water_tracer_transfers], no SHOC replay.
#              Nothing but P3 changes T_mid, so the rebuilt P3 entry
#              temperature must be the previous step's T_mid.
#  - mac:      [water_tracers, mac], mac = [shoc, cld_fraction, p3,
#              water_tracer_transfers] subcycled, both replays
# and the driver link (AtmosphereDriver::setup_water_tracer_transfer_links):
#  - link:     the mac run with no water_tracers_file for the transfers; the
#              link must copy it. Same checks, same numbers.
#  - neg_*:    mis-configured trees that the link must reject during init
#              (incl. water_tracers after shoc, SHOC's flux fixer on: SHOC
#              applies the tracer surface fluxes, and a process between SHOC
#              and P3 that updates T_mid).
#  - link_iso, link_iso_manual, neg_iso_*: the isotope options (only in builds
#              with EAMXX_ENABLE_WATER_ISOTOPES, at the end of this file)
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
# These need a process that updates bulk water or T_mid, which no linked physics does:
# wtt_link_modifier_test.cpp adds a test-only one (wtt_dummy_modifier) to the driver test
set (WTT_NEG_MOD_CASES modifier_qc modifier_group modifier_T)
set (WTT_NEG_inputs)
foreach (c ${WTT_NEG_CASES} ${WTT_NEG_MOD_CASES})
  list (APPEND WTT_NEG_inputs input_neg_${c}.yaml)
endforeach()

# The AD test driver, with the dummy modifier process registered
CreateUnitTestExec(wtt_link_modifier_test
  SOURCES wtt_link_modifier_test.cpp
  LIBS scream_control eamxx_io eamxx_diagnostics shoc cld_fraction p3 water_tracers)
target_include_directories(wtt_link_modifier_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

# The generic-tracer runs: tracers.yaml, no alpha override, no kinetic isotope
# fractionation (the isotope runs, at the end, set their own)
set (WTT_TRACERS_FILE tracers.yaml)
set (WTT_ALPHA_OVERRIDE 0.0)
set (WTT_KINETIC false)
set (WTT_KINETIC_ICE_S empirical)
set (OUT_YAML output_p3.yaml)
set (OUT_PREFIX p3_water_tracer_transfers_output)
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

# Tolerance of the roundoff checks, relative to the total water at a point,
# and of the stage temperature checks, relative to T_mid
if (SCREAM_DOUBLE_PRECISION)
  set (WTT_RTOL 1e-12)
  set (WTT_T_RTOL 1e-13)
else()
  set (WTT_RTOL 1e-5)
  set (WTT_T_RTOL 1e-5)
endif()
# dt of P3's calls in the mac runs
math (EXPR MAC_DT "${ATM_TIME_STEP}/${MAC_SUBCYCLES}")

set (WTT_CHECK ${CMAKE_CURRENT_SOURCE_DIR}/../shoc_water_tracer_transfers/check_water_tracer_transfers.py)
set (WTT_HEADERS
  --phases-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp
  --p3-header ${SCREAM_SRC_DIR}/physics/p3/p3_water_transfers.hpp
  --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
  --T-rtol ${WTT_T_RTOL})

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
    --dt ${P3_TIME_STEP}
    --T-from-previous-slice
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
    --dt ${MAC_DT}
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
    --dt ${MAC_DT}
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
    --dt ${MAC_DT}
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
    --dt ${MAC_DT}
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
set (WTT_NEG_REGEX_modifier_T     "'wtt_dummy_modifier'[ \n]+in[ \n]+group[ \n]+'mac'[ \n]+updates[ \n]+'T_mid'[ \n]+between[ \n]+SHOC/P3")
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

# Isotope runs (tracers_iso.yaml: HDO and H218O, next to a generic R = 1
# tracer), in builds that accept isotope entries:
#  - iso:           mac, equilibrium alphas, evap configuration (surface
#                   fluxes on). The tracers must be conserved (column +
#                   precipitation + surface flux) and stay >= 0 (to
#                   roundoff), condensing cloud must be enriched relative to
#                   the vapor, and ice that P3 mostly deposited from vapor
#                   must have R_qi/R_qv close to alpha_i (see
#                   check_water_isotope_transfers.py).
#  - iso_alpha1p5,
#    iso_alpha0p7:  the same with alpha_override = 1.5 and 0.7, far from any
#                   physical alpha: conservation and positivity only.
#  - iso_alpha1:    mac with alpha_override = 1, no ocean (ocnfrac = 0, so the
#                   isotopes' bulk-formula ocean flux, which has its own alpha,
#                   is 0 and their surface flux is the land part, as for
#                   wtrc_t): every isotope field must equal wtrc_t's bit for
#                   bit, so the fractionation code is exact at alpha = 1.
#  - p3_iso_alpha1: the same in the p3-only run (no SHOC, no surface flux).
#  - iso_kinetic,
#    iso_kinetic_model: the iso run with kinetic_fractionation on, with the
#                   empirical supersaturation over ice and with P3's RH over
#                   ice. Conservation and positivity; ice that P3 mostly
#                   deposited from vapor must have R_qi/R_qv close to the
#                   effective factor (akci below 253.15 K); and against the iso
#                   run (kinetic off, same bulk water bit for bit): where the
#                   rain evaporated below cloud into subsaturated air, the rain
#                   left is more enriched and has a lower deuterium excess, and
#                   cold ice deposited from vapor is less enriched relative to
#                   the vapor.
if (EAMXX_ENABLE_WATER_ISOTOPES)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_iso.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_iso.yaml COPYONLY)
  set (WTT_TRACERS_FILE tracers_iso.yaml)
  set (WTT_ISO_CHECK ${CMAKE_CURRENT_SOURCE_DIR}/../shoc_water_tracer_transfers/check_water_isotope_transfers.py)

  foreach (v iso iso_alpha1p5 iso_alpha0p7 iso_alpha1 iso_kinetic iso_kinetic_model)
    set (SURF_EVAP 1.5e-6)
    set (OCNFRAC   0.5)
    set (OCN_SSQ   0.01)
    set (WTT_KINETIC false)
    set (WTT_KINETIC_ICE_S empirical)
    if (v MATCHES "kinetic")
      set (WTT_ALPHA_OVERRIDE 0.0)
      set (WTT_KINETIC true)
      if (v MATCHES "model")
        set (WTT_KINETIC_ICE_S model)
      endif()
    elseif (v STREQUAL iso)
      set (WTT_ALPHA_OVERRIDE 0.0)
    elseif (v STREQUAL iso_alpha1p5)
      set (WTT_ALPHA_OVERRIDE 1.5)
    elseif (v STREQUAL iso_alpha0p7)
      set (WTT_ALPHA_OVERRIDE 0.7)
    else()
      set (WTT_ALPHA_OVERRIDE 1.0)
      set (OCNFRAC 0.0)
    endif()
    set (OUT_YAML output_${v}.yaml)
    set (OUT_PREFIX shoc_cld_p3_water_${v}_transfers_output)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml ${CMAKE_CURRENT_BINARY_DIR}/input_${v}.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_iso.yaml ${CMAKE_CURRENT_BINARY_DIR}/output_${v}.yaml)
    CreateUnitTestFromExec(shoc_cld_p3_water_${v}_transfers shoc_cld_p3_water_tracer_transfers
      EXE_ARGS "--args -ifile=input_${v}.yaml"
      LABELS shoc cld p3 water_tracers water_isotopes physics
      FIXTURES_SETUP shoc_cld_p3_water_${v}_transfers)

    set (out_nc shoc_cld_p3_water_${v}_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
    set (fixtures shoc_cld_p3_water_${v}_transfers)
    if (v STREQUAL iso_alpha1)
      # Bit for bit: every field of each isotope against wtrc_t's
      set (pairs)
      foreach (iso wtrc_hdo wtrc_o18)
        foreach (f qv qc qr qi precip_liq_surf_mass precip_ice_surf_mass residual_shoc residual_p3 qv_surf_flux)
          list (APPEND pairs "${iso}_${f}=wtrc_t_${f}")
        endforeach()
      endforeach()
      add_test (NAME shoc_cld_p3_water_${v}_transfers_bfb
        COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files -s ${out_nc} -c ${pairs}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      # and the uniform-ratio checks of the generic tracers
      add_test (NAME shoc_cld_p3_water_${v}_transfers_check
        COMMAND ${WTT_CHECK} -f ${out_nc}
          ${WTT_HEADERS}
          --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
          --tracer wtrc_t=1 wtrc_hdo=1 wtrc_o18=1
          --min-p3-transfers 8
          --rtol ${WTT_RTOL}
          --dt ${MAC_DT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      set_tests_properties (shoc_cld_p3_water_${v}_transfers_bfb PROPERTIES
        LABELS "shoc;cld;p3;water_tracers;water_isotopes;physics"
        FIXTURES_REQUIRED shoc_cld_p3_water_${v}_transfers)
    else()
      set (alpha_args)
      if (v STREQUAL iso)
        set (alpha_args --alpha-check)
      elseif (v MATCHES "kinetic")
        # against the iso run, which must have run first
        set (alpha_args --alpha-check --kinetic --kinetic-ice-S ${WTT_KINETIC_ICE_S}
          --equilibrium-file shoc_cld_p3_water_iso_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
        list (APPEND fixtures shoc_cld_p3_water_iso_transfers)
      endif()
      add_test (NAME shoc_cld_p3_water_${v}_transfers_check
        COMMAND ${WTT_ISO_CHECK} -f ${out_nc}
          --tracer wtrc_hdo=HDO wtrc_o18=H218O --reference wtrc_t
          --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
          --p3-header ${SCREAM_SRC_DIR}/physics/p3/p3_water_transfers.hpp
          --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
          --dt ${ATM_TIME_STEP} --dt-p3 ${MAC_DT}
          --rtol-budget ${WTT_RTOL} --neg-rtol ${WTT_RTOL}
          ${alpha_args}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()
    set_tests_properties (shoc_cld_p3_water_${v}_transfers_check PROPERTIES
      LABELS "shoc;cld;p3;water_tracers;water_isotopes;physics"
      FIXTURES_REQUIRED "${fixtures}")
  endforeach()

  set (WTT_ALPHA_OVERRIDE 1.0)
  set (OUT_YAML output_p3_iso_alpha1.yaml)
  set (OUT_PREFIX p3_water_iso_alpha1_transfers_output)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_p3.yaml ${CMAKE_CURRENT_BINARY_DIR}/input_p3_iso_alpha1.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_p3_iso.yaml ${CMAKE_CURRENT_BINARY_DIR}/output_p3_iso_alpha1.yaml)
  CreateUnitTestFromExec(p3_water_iso_alpha1_transfers shoc_cld_p3_water_tracer_transfers
    EXE_ARGS "--args -ifile=input_p3_iso_alpha1.yaml"
    LABELS p3 water_tracers water_isotopes physics
    FIXTURES_SETUP p3_water_iso_alpha1_transfers)
  set (out_nc p3_water_iso_alpha1_transfers_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
  set (pairs)
  foreach (iso wtrc_hdo wtrc_o18)
    foreach (f qv qc qr qi precip_liq_surf_mass precip_ice_surf_mass residual_p3)
      list (APPEND pairs "${iso}_${f}=wtrc_t_${f}")
    endforeach()
  endforeach()
  add_test (NAME p3_water_iso_alpha1_transfers_bfb
    COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files -s ${out_nc} -c ${pairs}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_test (NAME p3_water_iso_alpha1_transfers_check
    COMMAND ${WTT_CHECK} -f ${out_nc}
      ${WTT_HEADERS}
      --tracer wtrc_t=1 wtrc_hdo=1 wtrc_o18=1
      --min-p3-transfers 8
      --rtol ${WTT_RTOL}
      --dt ${P3_TIME_STEP}
      --T-from-previous-slice
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  foreach (t bfb check)
    set_tests_properties (p3_water_iso_alpha1_transfers_${t} PROPERTIES
      LABELS "p3;water_tracers;water_isotopes;physics"
      FIXTURES_REQUIRED p3_water_iso_alpha1_transfers)
  endforeach()

  # Driver link of the isotope options. AtmosphereDriver::setup_water_tracer_transfer_links
  # copies the six options of WaterIsotopes::configure from water_tracers to
  # water_tracer_transfers, unless the latter sets one itself (then they must agree).
  #  - link_iso:        mac run with the isotope tracers and kinetic_fractionation on,
  #                     the options set NON-DEFAULT on water_tracers only
  #  - link_iso_manual: the same, the options set by hand on both processes
  # The two outputs must be bit for bit equal (all variables, tolerance 0): the
  # replay's alphas depend on liquid_vapor_formulation, ice_vapor_formulation (isocam3
  # differs for HDO only) and diffusivity_formulation (kinetic is on), so a copy
  # missing in the driver would change the results. The other three options
  # (standard_ratio_formulation, ocean_enrichment_formulation,
  # kinetic_fractionation_formulation) are read by water_tracer_transfers but do
  # not enter its results (they act on the tracers' initial and ocean ratios and
  # on the ocean flux, in water_tracers): a missing copy of those is NOT visible
  # in the outputs, and only the negative tests below cover them (a mismatch
  # must stop the run). The checks here are conservation and positivity; the
  # alpha check is left out, since it assumes the default formulations.
  set (ISO_OPT_LIST
    "liquid_vapor_formulation: majoube_1971"
    "ice_vapor_formulation: isocam3"
    "diffusivity_formulation: cappa_2003"
    "kinetic_fractionation_formulation: none")
  function (wtt_iso_opts out indent)
    set (txt)
    foreach (o ${ARGN})
      string (APPEND txt "${indent}${o}\n")
    endforeach()
    string (REGEX REPLACE "\n$" "" txt "${txt}")
    set (${out} "${txt}" PARENT_SCOPE)
  endfunction()

  set (SURF_EVAP 1.5e-6)
  set (OCNFRAC   0.5)
  set (OCN_SSQ   0.01)
  foreach (v link_iso link_iso_manual)
    wtt_iso_opts (WTRC_ISO_OPTS "    " ${ISO_OPT_LIST})
    if (v STREQUAL link_iso)
      set (XFER_ISO_OPTS "")
    else()
      wtt_iso_opts (XFER_ISO_OPTS "      " ${ISO_OPT_LIST})
    endif()
    set (OUT_YAML output_${v}.yaml)
    set (OUT_PREFIX shoc_cld_p3_wtt_${v}_output)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_link_iso.yaml ${CMAKE_CURRENT_BINARY_DIR}/input_${v}.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_iso.yaml ${CMAKE_CURRENT_BINARY_DIR}/output_${v}.yaml)
    CreateUnitTestFromExec(shoc_cld_p3_wtt_${v} shoc_cld_p3_water_tracer_transfers
      EXE_ARGS "--args -ifile=input_${v}.yaml"
      LABELS shoc cld p3 water_tracers water_isotopes physics
      FIXTURES_SETUP shoc_cld_p3_wtt_${v})
  endforeach()
  set (link_iso_nc        shoc_cld_p3_wtt_link_iso_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
  set (link_iso_manual_nc shoc_cld_p3_wtt_link_iso_manual_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
  add_test (NAME shoc_cld_p3_wtt_link_iso_check
    COMMAND ${WTT_ISO_CHECK} -f ${link_iso_nc}
      --tracer wtrc_hdo=HDO wtrc_o18=H218O --reference wtrc_t
      --stages-header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
      --p3-header ${SCREAM_SRC_DIR}/physics/p3/p3_water_transfers.hpp
      --stage-T-header ${SCREAM_SRC_DIR}/physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp
      --dt ${ATM_TIME_STEP} --dt-p3 ${MAC_DT}
      --rtol-budget ${WTT_RTOL} --neg-rtol ${WTT_RTOL}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  # every output field (compare-nc-files without -c would also compare 'time', which it cannot)
  set (iso_link_fields)
  file (STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/output_iso.yaml iso_out_lines REGEX "^  - ")
  foreach (l ${iso_out_lines})
    string (REGEX REPLACE "^  - +" "" f "${l}")
    list (APPEND iso_link_fields "${f}=${f}")
  endforeach()
  add_test (NAME shoc_cld_p3_wtt_link_iso_bfb
    COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files -s ${link_iso_nc} -t ${link_iso_manual_nc} -tol 0
      -c ${iso_link_fields}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (shoc_cld_p3_wtt_link_iso_check PROPERTIES
    LABELS "shoc;cld;p3;water_tracers;water_isotopes;physics"
    FIXTURES_REQUIRED shoc_cld_p3_wtt_link_iso)
  set_tests_properties (shoc_cld_p3_wtt_link_iso_bfb PROPERTIES
    LABELS "shoc;cld;p3;water_tracers;water_isotopes;physics"
    FIXTURES_REQUIRED "shoc_cld_p3_wtt_link_iso;shoc_cld_p3_wtt_link_iso_manual")

  # negative tests: an isotope option set on both processes, to different values,
  # must stop the run in the link. One per kind: an option that changes the
  # replay's results (ice_vapor_formulation) and one that does not
  # (ocean_enrichment_formulation).
  foreach (c ice_vapor ocean_enrichment)
    if (c STREQUAL ice_vapor)
      set (WTRC_ISO_OPTS "    ice_vapor_formulation: isocam3")
      set (XFER_ISO_OPTS "      ice_vapor_formulation: merlivat_nief_1967")
      set (WTT_NEG_REGEX_iso_${c} "'water_tracer_transfers::ice_vapor_formulation'[ \n]+\\(merlivat_nief_1967\\)[ \n]+differs[ \n]+from[ \n]+'water_tracers::ice_vapor_formulation'[ \n]+\\(isocam3\\)")
    else()
      set (WTRC_ISO_OPTS "    ocean_enrichment_formulation: LGM")
      set (XFER_ISO_OPTS "      ocean_enrichment_formulation: modern")
      set (WTT_NEG_REGEX_iso_${c} "'water_tracer_transfers::ocean_enrichment_formulation'[ \n]+\\(modern\\)[ \n]+differs[ \n]+from[ \n]+'water_tracers::ocean_enrichment_formulation'[ \n]+\\(LGM\\)")
    endif()
    set (OUT_YAML output_link_iso.yaml)
    set (OUT_PREFIX shoc_cld_p3_wtt_link_iso_neg_output)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_link_iso.yaml ${CMAKE_CURRENT_BINARY_DIR}/input_neg_iso_${c}.yaml)
    CreateUnitTestFromExec(wtt_link_neg_iso_${c} shoc_cld_p3_water_tracer_transfers
      EXE_ARGS "--args -ifile=input_neg_iso_${c}.yaml"
      LABELS water_tracers water_isotopes physics)
    set_tests_properties (wtt_link_neg_iso_${c} PROPERTIES
      PASS_REGULAR_EXPRESSION "${WTT_NEG_REGEX_iso_${c}}")
  endforeach()
endif()
