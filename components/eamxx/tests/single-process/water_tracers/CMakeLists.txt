include (ScreamUtils)

# Standalone AD runs of the water_tracers process on a mesh-free point grid.
# All inputs are constant ICs (no input files). The process only computes the
# tracer surface fluxes (SHOC applies them), so tracers and qv stay at their
# initial values.

set (NUM_STEPS 2)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)
# output.yaml skips the t0 output, so the output file is stamped one step later
set (FIRST_OUTPUT_T 2021-10-12-46800)

CreateADUnitTestExec(water_tracers_standalone LIBS water_tracers)

# Generic tracers only.
set (POSTFIX generic)
set (TRACERS_FILE tracers_generic.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_generic.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/tracers_generic.yaml COPYONLY)
CreateUnitTestFromExec(water_tracers_standalone water_tracers_standalone
  EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
  LABELS water_tracers physics
  FIXTURES_SETUP water_tracers_standalone_generic)

# wtrc_bf has init: ratio and r_init = r_ocn = 1, and nothing changes it after
# initialization, so every phase must equal its bulk phase.
add_test (NAME water_tracers_standalone_phase_ic_check
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s water_tracers_standalone_output_generic.INSTANT.nsteps_x1.np1.${FIRST_OUTPUT_T}.nc
          -c "wtrc_bf_qv=qv" "wtrc_bf_qc=qc" "wtrc_bf_qr=qr" "wtrc_bf_qi=qi"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (water_tracers_standalone_phase_ic_check PROPERTIES
          LABELS "water_tracers;physics"
          FIXTURES_REQUIRED water_tracers_standalone_generic)

# The tracer surface fluxes (check_surface_fluxes.py): an R = 1 proportional tracer
# gets exactly surf_evap, and R_ocn != R_land tracers a mix of the two parts.
add_test (NAME water_tracers_standalone_surface_flux_check
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_surface_fluxes.py
          -f water_tracers_standalone_output_generic.INSTANT.nsteps_x1.np1.${FIRST_OUTPUT_T}.nc
          --surf-evap 2.0e-5 --ocnfrac 0.5 --ssq 0.02 --re 0.0015 --ustar 0.3
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (water_tracers_standalone_surface_flux_check PROPERTIES
          LABELS "water_tracers;physics"
          FIXTURES_REQUIRED water_tracers_standalone_generic)

if (EAMXX_ENABLE_WATER_ISOTOPES)
  # Same tracers, plus isotope entries in the same file/process.
  set (POSTFIX isotopes)
  set (TRACERS_FILE tracers_with_isotopes.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_with_isotopes.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_with_isotopes.yaml COPYONLY)
  CreateUnitTestFromExec(water_tracers_standalone_isotopes water_tracers_standalone
    EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
    LABELS water_tracers water_isotopes physics)
else()
  # Negative test: a build without isotope support must reject isotope
  # entries with a clear build-time-aware error, not silently ignore them
  # or require a different process name.
  set (POSTFIX isotopes_rejected)
  set (TRACERS_FILE tracers_with_isotopes.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_with_isotopes.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_with_isotopes.yaml COPYONLY)
  CreateUnitTestFromExec(water_tracers_standalone_isotopes_rejected water_tracers_standalone
    EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
    LABELS water_tracers physics)
  set_tests_properties(water_tracers_standalone_isotopes_rejected PROPERTIES
    WILL_FAIL TRUE)
endif()
