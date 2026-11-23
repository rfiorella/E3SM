include (ScreamUtils)

set (TEST_BASE_NAME shoc_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS shoc physics
  LIBS shoc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR NUM_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS shoc physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

# Check tendency calculation
foreach (NRANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
  set (fname shoc_standalone_output.INSTANT.nsteps_x1.np${NRANKS}.${RUN_T0}.nc)
  set (tname shoc_tend_check_np${NRANKS})
  add_test (NAME ${tname}
    COMMAND ${script} -f ${fname}
    -v T_mid qv tke horiz_winds sgs_buoy_flux eddy_diff_mom qc cldfrac_liq
    -t shoc_T_mid_tend shoc_qv_tend shoc_tke_tend shoc_horiz_winds_tend shoc_sgs_buoy_flux_tend shoc_eddy_diff_mom_tend shoc_qc_tend shoc_cldfrac_liq_tend
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (${tname} PROPERTIES
    LABELS "shoc;physics"
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${NRANKS}_omp1)
endforeach()

# Check that the content of sliced vars U/V and corresponding surf_mom_flux
# components in all output NC files matches the content of corresponding components
# of the vector quantities
foreach (RANK RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (nc_file shoc_standalone_output.INSTANT.nsteps_x1.np${RANK}.${RUN_T0}.nc)
  add_test (NAME check_U_V_slices_np${RANK}
            COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
            -s ${nc_file} -c "horiz_winds(:,:,1,:)=U(:,:,:)" "horiz_winds(:,:,2,:)=V(:,:,:)"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (check_U_V_slices_np${RANK} PROPERTIES
            FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${RANK}_omp1)
  add_test (NAME check_surf_mom_flux_slices_np${RANK}
            COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
            -s ${nc_file} -c "surf_mom_flux(:,:,1)=surf_mom_flux_U(:,:)"
                             "surf_mom_flux(:,:,2)=surf_mom_flux_V(:,:)"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (check_surf_mom_flux_slices_np${RANK} PROPERTIES
      LABELS "shoc;physics"
      FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${RANK}_omp1)
endforeach()

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME})
endif()

################################
#        MUST FAIL tests       #
################################

# These tests ensure that the compare-nc-files script doesn't always pass
set (nc_file shoc_standalone_output.INSTANT.nsteps_x1.np${TEST_RANK_START}.${RUN_T0}.nc)

# Legitimate diff
add_test (NAME check_U_V_slices_fail_diff
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "horiz_winds(:,:,1,:)=T_mid(:,:,:)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_diff PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_START}_omp1)

# Wrong layout
add_test (NAME check_U_V_slices_fail_layout
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "horiz_winds(:,:,1,:)=U(:,1,:)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_layout PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_START}_omp1)

# Missing variable(s)
add_test (NAME check_U_V_slices_fail_missing
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "T=T"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_missing PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_START}_omp1)

# In builds with water tracers SHOC exports bulk water after turbulent transport
# and after condensation (shoc_water_stages), which must reproduce the qv and qc
# SHOC returns (see check_water_stages.py).
if (EAMXX_ENABLE_WATER_TRACERS)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_water_transfers.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_water_transfers.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_water_transfers.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_water_transfers.yaml)
  CreateUnitTestFromExec(shoc_standalone_water_transfers shoc_standalone
    EXE_ARGS "--args -ifile=input_water_transfers.yaml"
    LABELS shoc physics
    FIXTURES_SETUP shoc_standalone_water_transfers)

  # output_water_transfers.yaml skips t0, so its file is stamped one step later
  math (EXPR WT_FIRST_OUTPUT_SECONDS "45000 + ${ATM_TIME_STEP}")
  add_test (NAME shoc_water_stages_bookkeeping_check
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_stages.py
      -f shoc_standalone_water_transfers_output.INSTANT.nsteps_x1.np1.2021-10-12-${WT_FIRST_OUTPUT_SECONDS}.nc
      --header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (shoc_water_stages_bookkeeping_check PROPERTIES
    LABELS "shoc;physics"
    FIXTURES_REQUIRED shoc_standalone_water_transfers)
endif()

# The water tracers are diffused in the same SHOC solve as bulk water, so a
# tracer with R = 1 must come out of the first SHOC call with exactly the
# transported state that shoc_water_stages records.
if (EAMXX_ENABLE_WATER_TRACERS)
  set (WT_TIME_STEP 300)
  set (WT_NUM_STEPS 2)
  CreateADUnitTestExec(shoc_water_tracers_standalone LIBS shoc water_tracers)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_water_transport.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_water_transport.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_water_transport.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_water_transport.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_water_transport.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_water_transport.yaml COPYONLY)
  CreateUnitTestFromExec(shoc_water_tracers_transport shoc_water_tracers_standalone
    EXE_ARGS "--args -ifile=input_water_transport.yaml"
    LABELS shoc water_tracers physics
    FIXTURES_SETUP shoc_water_tracers_transport)

  add_test (NAME shoc_water_stages_transport_check
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_water_stages.py
      -f shoc_water_tracers_transport_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
      --header ${SCREAM_SRC_DIR}/physics/shoc/shoc_water_transfers.hpp
      --tracer wtrc_t
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (shoc_water_stages_transport_check PROPERTIES
    LABELS "shoc;water_tracers;physics"
    FIXTURES_REQUIRED shoc_water_tracers_transport)
endif()
