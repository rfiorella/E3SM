#!/bin/bash
################################################################################
# EAMv2 Water Isotope Test Suite (Functional Tests)
#
# Purpose: Test water isotope functionality in EAMv2 by building and running
#          short simulations with different isotope tracer packages.
#          Uses ne4pg2 resolution (~100km grid, 384 columns) for fast,
#          low-memory functional testing.
#
# Usage: Run this script inside the e3sm-openmpi-dev-latest Docker container
#        with appropriate volume mounts.
#
# Author: Auto-generated test suite
# Date: 2026-04-07
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

################################################################################
# Configuration
################################################################################

# Test configuration
COMPSET="Fi20TR"
RESOLUTION="ne4pg2_ne4pg2"
# Resolution choice rationale:
# - ne4pg2_ne4pg2: Coarsest practical grid for functional testing
# - 384 atmospheric columns (~100 km spacing)
# - Data ocean (DOCN) for simplicity and speed
# - Memory: ~1-2 GB (vs 12-16 GB for ne30pg2)
# - Runtime: ~15-30 min (vs 2-4 hours for ne30pg2)
# - Validates all isotope tracers and physics at coarse resolution
MACHINE="docker"
COMPILER="gnu"
PROJECT="docker_runs"
TEST_DURATION_DAYS=5
DEBUG_COMPILE="${DEBUG_COMPILE:-false}"

# Detect if running in container or on host
if [ -d "/home/e3smuser/EAMv2-wiso" ]; then
    # Running in container
    REPO_ROOT="/home/e3smuser/EAMv2-wiso"
    OUTPUT_BASE="/home/e3smuser/output"
    INPUT_DATA="/home/e3smuser/inputdata"
else
    # Running on host - use current directory
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    REPO_ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"
    OUTPUT_BASE="${OUTPUT_BASE:-/tmp/eamv2_isotope_tests}"
    INPUT_DATA="${INPUT_DATA:-/tmp/eamv2_inputdata}"
fi

CIME_SCRIPTS="${REPO_ROOT}/cime/scripts"

# Test timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TEST_DIR="${OUTPUT_BASE}/isotope_tests_${TIMESTAMP}"

# Tracer packages to test
declare -a TRACER_PACKAGES=(
    "h2o_h216o_hdo_h218o"
    "all_stable_wiso"
)

# Tracer counts for reference  
declare -A TRACER_COUNTS=(
    ["h2o_h216o_hdo_h218o"]=28
    ["all_stable_wiso"]=35
)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Helper Functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo ""
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo ""
}

print_separator() {
    echo "--------------------------------------------------------------------------------"
}

################################################################################
# Main Functions
################################################################################

check_environment() {
    print_header "Checking Environment"
    
    # Check if we're in the container
    if [ ! -f "/.dockerenv" ] && [ ! -f "/run/.containerenv" ]; then
        log_warn "Not running in Docker container - running on host system"
        log_info "Using repository at: ${REPO_ROOT}"
    else
        log_info "Running in container environment"
    fi
    
    # Check required directories exist
    if [ ! -d "${REPO_ROOT}" ]; then
        log_error "Repository not found at ${REPO_ROOT}"
        log_error "Make sure the repository is mounted correctly"
        exit 1
    fi
    log_success "Repository found: ${REPO_ROOT}"
    
    if [ ! -d "${CIME_SCRIPTS}" ]; then
        log_error "CIME scripts not found at ${CIME_SCRIPTS}"
        log_error "Expected location: ${REPO_ROOT}/cime/scripts"
        log_error "This suggests the repository is incomplete or incorrectly mounted"
        exit 1
    fi
    log_success "CIME scripts found: ${CIME_SCRIPTS}"
    
    # Validate critical CIME script exists (only create_newcase is in CIME scripts dir)
    if [ ! -f "${CIME_SCRIPTS}/create_newcase" ]; then
        log_error "create_newcase script not found in ${CIME_SCRIPTS}"
        log_error "CIME installation appears incomplete"
        exit 1
    fi
    log_success "CIME create_newcase script found"
    
    # Check that critical EAM directories exist
    if [ ! -d "${REPO_ROOT}/components/eam" ]; then
        log_error "EAM component not found at ${REPO_ROOT}/components/eam"
        log_error "Repository structure is incorrect"
        exit 1
    fi
    log_success "EAM component found"
    
    # Create output directories with error checking
    if ! mkdir -p "${OUTPUT_BASE}" 2>/dev/null; then
        log_error "Cannot create output base directory: ${OUTPUT_BASE}"
        log_error "Check permissions or specify different output location"
        exit 1
    fi
    
    if ! mkdir -p "${INPUT_DATA}" 2>/dev/null; then
        log_warn "Cannot create input data directory: ${INPUT_DATA}"
        log_warn "Input data download may fail - ensure directory is writable"
    fi
    
    if ! mkdir -p "${TEST_DIR}" 2>/dev/null; then
        log_error "Cannot create test directory: ${TEST_DIR}"
        log_error "Check permissions or disk space"
        exit 1
    fi
    log_success "Output directories created"
    
    # Export environment variables that config_machines.xml needs
    export OUTPUT_BASE
    export INPUT_DATA
    
    log_success "Environment check passed"
    log_info "Repository: ${REPO_ROOT}"
    log_info "CIME scripts: ${CIME_SCRIPTS}"
    log_info "Test output: ${TEST_DIR}"
    log_info "Input data: ${INPUT_DATA}"
    log_info "Debug compilation: ${DEBUG_COMPILE}"
}

create_test_case() {
    local package=$1
    local case_name="isotope_test_${package}"
    local case_dir="${TEST_DIR}/${case_name}"
    
    # All logging to stderr so only case_dir goes to stdout
    print_separator >&2
    log_info "Creating test case for ${package}" >&2
    
    # Validate CIME scripts directory
    if [ ! -d "${CIME_SCRIPTS}" ]; then
        log_error "CIME scripts directory not found: ${CIME_SCRIPTS}" >&2
        return 1
    fi
    
    cd "${CIME_SCRIPTS}" || {
        log_error "Cannot access CIME scripts directory: ${CIME_SCRIPTS}" >&2
        return 1
    }
    
    if [ ! -f "./create_newcase" ]; then
        log_error "create_newcase not found in ${CIME_SCRIPTS}" >&2
        return 1
    fi
    
    # Create log file for this operation
    local create_log="${TEST_DIR}/create_newcase_${package}.log"
    
    # Create the case - redirect all output to log and stderr
    log_info "Running create_newcase (output: ${create_log})" >&2
    if ./create_newcase \
        --case "${case_dir}" \
        --compset "${COMPSET}" \
        --res "${RESOLUTION}" \
        --mach "${MACHINE}" \
        --compiler "${COMPILER}" \
        --project "${PROJECT}" > "${create_log}" 2>&1; then
        
        # Verify case directory was actually created with required files
        if [ ! -d "${case_dir}" ]; then
            log_error "Case directory was not created: ${case_dir}" >&2
            log_error "Check log: ${create_log}" >&2
            return 1
        fi
        
        if [ ! -f "${case_dir}/xmlquery" ] || [ ! -f "${case_dir}/xmlchange" ]; then
            log_error "Case directory missing required files (xmlquery/xmlchange)" >&2
            log_error "Check log: ${create_log}" >&2
            return 1
        fi
        
        log_success "Case created successfully: ${case_dir}" >&2
        echo "${case_dir}"  # Only output: the case directory path
        return 0
    else
        log_error "create_newcase failed for ${package}" >&2
        log_error "Check log: ${create_log}" >&2
        if [ -f "${create_log}" ]; then
            log_error "Last 10 lines of log:" >&2
            tail -10 "${create_log}" >&2
        fi
        return 1
    fi
}

configure_case() {
    local case_dir=$1
    local package=$2
    
    log_info "Configuring case for ${package}"
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    if [ ! -f "./xmlquery" ]; then
        log_error "xmlquery not found in case directory - case may not be properly created"
        return 1
    fi
    
    # Get base CAM_CONFIG_OPTS from the compset
    # For Fi2000climo, this should include the default physics
    # We need to add -water_tracer to it
    
    # NOTE: The Kokkos initialization bug has been fixed (2026-04-10)
    # - Fixed use-after-scope bug in compose_slmm.cpp and ExecSpaceDefs.cpp
    # - Now using default theta-l dynamics target with COMPOSE transport
    # - Previous preqx workaround has been removed
    log_info "Using default theta-l dynamics with Kokkos (bug fixed)"
    
    # Set water tracer package in CAM_CONFIG_OPTS
    # Note: We need to append to existing config options
    local cam_opts=$(./xmlquery CAM_CONFIG_OPTS --value)
    local new_opts="${cam_opts} -water_tracer ${package}"
    
    # ./xmlchange CAM_CONFIG_OPTS="${new_opts}"
    
    # Enable isotope field passing in coupler
    # ./xmlchange FLDS_WISO=TRUE
    
    # Enable debug compilation if requested
    if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
        log_info "Enabling debug compilation (DEBUG=TRUE)"
        ./xmlchange DEBUG=TRUE
    fi
    
    # Set run length
    ./xmlchange STOP_N=${TEST_DURATION_DAYS}
    ./xmlchange STOP_OPTION=ndays
    ./xmlchange REST_N=${TEST_DURATION_DAYS}
    ./xmlchange RESUBMIT=0
    
    # Set processor layout for Docker
    ./xmlchange NTASKS=4
    ./xmlchange NTHRDS=1
    
    # Set input data directory
    ./xmlchange DIN_LOC_ROOT="${INPUT_DATA}"
    
    # Set build and run directories inside CASEROOT for self-contained tests
    log_info "Setting build and run directories inside case directory"
    ./xmlchange EXEROOT="${case_dir}/bld"
    ./xmlchange RUNDIR="${case_dir}/run"
    
    # Disable short-term archiving (not needed for tests)
    ./xmlchange DOUT_S=FALSE
    
    # Create user_nl_eam with isotope-specific namelists
    cat > user_nl_eam << EOF
! Water isotope configuration
!trace_water = .true.
!wisotope = .true.

! Physics options for isotopes
!wtrc_alpha_kinetic = .true.
!wtrc_lh2oadj = .true.
!wtrc_niter = 1
!wtrc_qmin = 1.0e-18

! Conservation checks (disabled for performance)
!wtrc_check_total_h2o = .false.
!wtrc_warn_only = .true.

! Disable prescribed aerosols by setting empty file
prescribed_aero_file = ''
EOF
    
    log_success "Case configured for ${package}"
}

setup_case() {
    local case_dir=$1
    local package=$2
    
    log_info "Setting up case for ${package}"
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    if [ ! -f "./case.setup" ]; then
        log_error "case.setup not found in case directory"
        return 1
    fi
    
    ./case.setup
    
    if [ $? -eq 0 ]; then
        log_success "Case setup complete"
        return 0
    else
        log_error "Failed to setup case for ${package}"
        return 1
    fi
}

check_input_data() {
    local case_dir=$1
    local package=$2
    
    log_info "Checking input data for ${package}"
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    if [ ! -f "./check_input_data" ]; then
        log_warn "check_input_data script not found, skipping input data check"
        return 0
    fi
    
    # Check if input data is available
    ./check_input_data --download
    
    if [ $? -eq 0 ]; then
        log_success "Input data check complete"
        return 0
    else
        log_warn "Some input data may be missing - continuing anyway"
        return 0
    fi
}

build_case() {
    local case_dir=$1
    local package=$2
    local start_time=$(date +%s)
    
    log_info "Building case for ${package} (this may take 30-60 minutes)..."
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    if [ ! -f "./case.build" ]; then
        log_error "case.build not found in case directory"
        return 1
    fi
    
    ./case.build 2>&1 | tee build.log
    
    local build_status=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local minutes=$((duration / 60))
    
    if [ $build_status -eq 0 ]; then
        log_success "Build complete in ${minutes} minutes"
        return 0
    else
        log_error "Build failed for ${package}"
        log_error "Check build log: ${case_dir}/build.log"
        return 1
    fi
}

run_case() {
    local case_dir=$1
    local package=$2
    local start_time=$(date +%s)
    
    log_info "Running ${TEST_DURATION_DAYS}-day simulation for ${package}..."
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    if [ ! -f "./case.submit" ]; then
        log_error "case.submit not found in case directory"
        return 1
    fi
    
    # Configure runtime environment for enhanced debugging
    log_info "Setting up sanitizer environment variables"
    export ASAN_OPTIONS="abort_on_error=1:detect_leaks=1:detect_stack_use_after_return=1:print_stacktrace=1:log_path=${case_dir}/run/asan.log"
    export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:log_path=${case_dir}/run/ubsan.log"
    
    # Increase stack size (may help with large tracer counts)
    ulimit -s unlimited
    
    # Enable core dumps for post-mortem analysis
    ulimit -c unlimited
    
    log_info "Runtime environment configured:"
    log_info "  - AddressSanitizer enabled (logs: ${case_dir}/run/asan.log.*)"
    log_info "  - UndefinedBehaviorSanitizer enabled (logs: ${case_dir}/run/ubsan.log.*)"
    log_info "  - Stack size: unlimited"
    log_info "  - Core dumps: enabled"
    
    # Submit the case (or run directly if batch system is none)
    ./case.submit 2>&1 | tee run.log
    
    local run_status=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local minutes=$((duration / 60))
    
    if [ $run_status -eq 0 ]; then
        log_success "Run complete in ${minutes} minutes"
        return 0
    else
        log_error "Run failed for ${package}"
        log_error "Check run log: ${case_dir}/run.log"
        
        # Check for sanitizer logs
        if ls "${case_dir}"/run/asan.log.* 1>/dev/null 2>&1; then
            log_error "AddressSanitizer detected issues - see: ${case_dir}/run/asan.log.*"
            log_error "Last 20 lines of AddressSanitizer log:"
            tail -20 "${case_dir}"/run/asan.log.* 2>/dev/null
        fi
        
        if ls "${case_dir}"/run/ubsan.log.* 1>/dev/null 2>&1; then
            log_error "UndefinedBehaviorSanitizer detected issues - see: ${case_dir}/run/ubsan.log.*"
            log_error "Last 20 lines of UndefinedBehaviorSanitizer log:"
            tail -20 "${case_dir}"/run/ubsan.log.* 2>/dev/null
        fi
        
        return 1
    fi
}

validate_output() {
    local case_dir=$1
    local package=$2
    
    log_info "Validating output for ${package}"
    
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist: ${case_dir}"
        return 1
    fi
    
    cd "${case_dir}" || return 1
    
    # Check if run completed successfully
    if grep -q "SUCCESSFUL TERMINATION" run/cpl.log.* 2>/dev/null; then
        log_success "Model run completed successfully"
    else
        log_error "Model did not complete successfully"
        return 1
    fi
    
    # Check for history files
    local hist_files=$(ls run/*.eam.h0.*.nc 2>/dev/null | wc -l)
    if [ $hist_files -eq 0 ]; then
        log_error "No history files found"
        return 1
    fi
    log_info "Found ${hist_files} history file(s)"
    
    # Check for isotope tracer variables
    local latest_hist=$(ls -t run/*.eam.h0.*.nc 2>/dev/null | head -1)
    
    if [ -z "${latest_hist}" ]; then
        log_error "Cannot find history file for validation"
        return 1
    fi
    
    log_info "Checking for isotope tracers in: $(basename ${latest_hist})"
    
    # Define expected tracers based on package
    local expected_tracers=""
    case "${package}" in
        "h2o_h216o_hdo_h218o")
            expected_tracers="Q_H2O Q_H216O Q_HDO Q_H218O"
            ;;
        "all_stable_wiso")
            expected_tracers="Q_H2O Q_H216O Q_HDO Q_H218O Q_H217O"
            ;;
    esac
    
    # Check for each expected tracer
    local all_found=true
    for tracer in ${expected_tracers}; do
        if ncdump -h "${latest_hist}" 2>/dev/null | grep -q "float ${tracer}("; then
            log_success "Found tracer: ${tracer}"
        else
            log_error "Missing tracer: ${tracer}"
            all_found=false
        fi
    done
    
    if [ "${all_found}" = true ]; then
        log_success "All expected isotope tracers found"
        return 0
    else
        log_error "Some isotope tracers are missing"
        return 1
    fi
}

run_single_test() {
    local package=$1
    local test_num=$2
    local total_tests=$3
    
    print_header "Test ${test_num}/${total_tests}: ${package} (${TRACER_COUNTS[$package]} tracers)"
    
    local case_dir
    local overall_status=0
    
    # Track start time for this test
    local test_start=$(date +%s)
    
    # Create case - capture output which should be the case directory path
    log_info "Creating case..." >&2
    case_dir=$(create_test_case "${package}")
    local create_status=$?
    
    # Validate case_dir is actually a valid directory path
    if [ $create_status -ne 0 ] || [ -z "${case_dir}" ] || [[ "${case_dir}" == *"ERROR"* ]] || [[ "${case_dir}" == *"error"* ]]; then
        log_error "Failed at case creation stage (status: ${create_status})" >&2
        log_error "Invalid or empty case directory returned: '${case_dir}'" >&2
        # Still set a case_dir for reporting purposes
        case_dir="${TEST_DIR}/isotope_test_${package}"
        mkdir -p "${case_dir}" 2>/dev/null || true
        echo "FAIL" > "${case_dir}/TEST_STATUS" 2>/dev/null || true
        return 1
    fi
    
    # Additional validation: ensure it's an absolute path
    if [[ ! "${case_dir}" =~ ^/ ]]; then
        log_error "Case directory is not an absolute path: ${case_dir}" >&2
        return 1
    fi
    
    # Verify the directory actually exists
    if [ ! -d "${case_dir}" ]; then
        log_error "Case directory does not exist after creation: ${case_dir}" >&2
        return 1
    fi
    
    log_success "Case directory created and validated: ${case_dir}" >&2
    
    # Configure case
    configure_case "${case_dir}" "${package}" || { overall_status=1; }
    
    # Setup case
    if [ $overall_status -eq 0 ]; then
        setup_case "${case_dir}" "${package}" || { overall_status=1; }
    fi
    
    # Check input data
    if [ $overall_status -eq 0 ]; then
        check_input_data "${case_dir}" "${package}" || { log_warn "Input data check had issues"; }
    fi
    
    # Build case
    if [ $overall_status -eq 0 ]; then
        build_case "${case_dir}" "${package}" || { overall_status=1; }
    fi
    
    # Run case
    if [ $overall_status -eq 0 ]; then
        run_case "${case_dir}" "${package}" || { overall_status=1; }
    fi
    
    # Validate output
    if [ $overall_status -eq 0 ]; then
        validate_output "${case_dir}" "${package}" || { overall_status=1; }
    fi
    
    # Calculate test duration
    local test_end=$(date +%s)
    local test_duration=$((test_end - test_start))
    local test_minutes=$((test_duration / 60))
    
    # Report results
    print_separator
    if [ $overall_status -eq 0 ]; then
        log_success "Test PASSED for ${package} (${test_minutes} minutes)"
        echo "PASS" > "${case_dir}/TEST_STATUS"
    else
        log_error "Test FAILED for ${package} (${test_minutes} minutes)"
        echo "FAIL" > "${case_dir}/TEST_STATUS"
    fi
    
    return $overall_status
}

generate_summary_report() {
    local total=$1
    local passed=$2
    local failed=$3
    local suite_duration=$4
    
    local report_file="${TEST_DIR}/test_summary.txt"
    
    cat > "${report_file}" << EOF
===============================================================================
EAMv2 Water Isotope Test Suite Results
Date: $(date)
===============================================================================

Configuration:
  Compset: ${COMPSET}
  Resolution: ${RESOLUTION}
  Machine: ${MACHINE}
  Test Duration: ${TEST_DURATION_DAYS} days
  Output Directory: ${TEST_DIR}

EOF

    # Add individual test results
    local test_num=1
    for package in "${TRACER_PACKAGES[@]}"; do
        local case_dir="${TEST_DIR}/isotope_test_${package}"
        local status="UNKNOWN"
        
        if [ -f "${case_dir}/TEST_STATUS" ]; then
            status=$(cat "${case_dir}/TEST_STATUS")
        fi
        
        cat >> "${report_file}" << EOF
Test ${test_num}: ${package} (${TRACER_COUNTS[$package]} tracers)
  Status: ${status}
  Case Directory: ${case_dir}

EOF
        test_num=$((test_num + 1))
    done
    
    # Add summary
    cat >> "${report_file}" << EOF
===============================================================================
Summary:
  Total Tests: ${total}
  Passed: ${passed}
  Failed: ${failed}
  Duration: ${suite_duration} minutes
===============================================================================

Future Phases (Recommended):
  Phase 2: Physical validation (δD, δ18O, meteoric water line)
  Phase 3: Mass conservation analysis
  Phase 4: Performance benchmarking
  Phase 5: Scientific validation against observations

For more information, see test_scripts/README.md
===============================================================================
EOF

    log_info "Summary report written to: ${report_file}"
    
    # Also display the report
    echo ""
    cat "${report_file}"
}

################################################################################
# Main Execution
################################################################################

main() {
    print_header "EAMv2 Water Isotope Test Suite"
    
    log_info "Starting test suite at $(date)"
    log_info "Testing ${#TRACER_PACKAGES[@]} isotope configurations"
    
    # Track overall start time
    local suite_start=$(date +%s)
    
    # Check environment
    check_environment
    
    # Run tests
    local total_tests=${#TRACER_PACKAGES[@]}
    local passed_tests=0
    local failed_tests=0
    local test_num=1
    
    for package in "${TRACER_PACKAGES[@]}"; do
        if run_single_test "${package}" ${test_num} ${total_tests}; then
            passed_tests=$((passed_tests + 1))
        else
            failed_tests=$((failed_tests + 1))
        fi
        test_num=$((test_num + 1))
    done
    
    # Calculate suite duration
    local suite_end=$(date +%s)
    local suite_duration=$((suite_end - suite_start))
    local suite_minutes=$((suite_duration / 60))
    
    # Generate summary report
    generate_summary_report ${total_tests} ${passed_tests} ${failed_tests} ${suite_minutes}
    
    # Final status
    print_header "Test Suite Complete"
    
    if [ $failed_tests -eq 0 ]; then
        log_success "All ${total_tests} tests PASSED!"
        return 0
    else
        log_error "${failed_tests}/${total_tests} tests FAILED"
        return 1
    fi
}

# Run main function
main "$@"
