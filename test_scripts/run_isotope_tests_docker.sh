#!/bin/bash
################################################################################
# Docker Test Runner for EAMv2 Water Isotope Tests
#
# This wrapper script makes it easy to run the isotope test suite in Docker.
# It handles mounting the correct directories and running the tests.
#
# Usage:
#   ./run_isotope_tests_docker.sh [OPTIONS]
#
# Options:
#   --interactive, -i    Start interactive shell instead of running tests
#   --help, -h          Show this help message
#
# Author: Auto-generated test suite
# Date: 2026-04-07
################################################################################

set -e
set -u

################################################################################
# Configuration
################################################################################

# Docker image
DOCKER_IMAGE="${E3SM_DOCKER_IMAGE:-rfiorella/model-containers:e3sm-openmpi-dev-latest}"

# Auto-detect repository path (script location)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DEFAULT_REPO_PATH="$( cd "${SCRIPT_DIR}/.." && pwd )"

# Host paths - can be overridden by environment variables or command-line args
REPO_PATH="${E3SM_REPO_PATH:-${DEFAULT_REPO_PATH}}"
INPUT_DATA_PATH="${E3SM_INPUT_DATA:-/data/inputdata}"
OUTPUT_PATH="${E3SM_OUTPUT:-/data/eamv2}"

# Container paths (fixed)
CONTAINER_REPO="/home/e3smuser/EAMv2-wiso"
CONTAINER_INPUT="/home/e3smuser/inputdata"
CONTAINER_OUTPUT="/home/e3smuser/output"

# Test script path
TEST_SCRIPT="${CONTAINER_REPO}/test_scripts/test_isotopes.sh"

# Runtime options
DRY_RUN=false
DEBUG_COMPILE=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

################################################################################
# Functions
################################################################################

print_header() {
    echo -e "${BLUE}"
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo -e "${NC}"
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

show_help() {
    cat << EOF
EAMv2 Water Isotope Test Suite - Docker Runner

Usage:
    $0 [OPTIONS]

Options:
    --interactive, -i           Start interactive shell in container
    --input-data PATH          Set input data directory (default: ${INPUT_DATA_PATH})
    --output PATH              Set output directory (default: ${OUTPUT_PATH})
    --repo PATH                Set repository path (default: auto-detected)
    --image IMAGE              Set Docker image (default: ${DOCKER_IMAGE})
    --debug                    Enable debug compilation flags (slower but easier to debug)
    --dry-run                  Show what would be done without executing
    --help, -h                 Show this help message

Environment Variables:
    E3SM_INPUT_DATA            Input data directory (overrides default)
    E3SM_OUTPUT                Output directory (overrides default)
    E3SM_REPO_PATH             Repository path (overrides auto-detection)
    E3SM_DOCKER_IMAGE          Docker image to use

Description:
    This script runs the EAMv2 water isotope test suite in a Docker container.
    It automatically mounts the required directories and executes the tests.

Requirements:
    - Docker must be installed and running
    - Docker memory limit: 4-8 GB recommended (tests use ~1-2 GB with ne4pg2)
    - Sufficient disk space for input data (~20-50 GB) and output (~2-5 GB per test)

Directory Structure:
    Repository:  ${REPO_PATH}
                 → Mounted as: ${CONTAINER_REPO}
    
    Input Data:  ${INPUT_DATA_PATH}
                 → Mounted as: ${CONTAINER_INPUT}
                 → Created automatically if it doesn't exist
    
    Output:      ${OUTPUT_PATH}
                 → Mounted as: ${CONTAINER_OUTPUT}
                 → Created automatically if it doesn't exist

Test Duration (ne4pg2 functional tests):
    - First run: 30-90 minutes (includes input data download)
    - Subsequent runs: 15-30 minutes (reuses cached input data)
    - Note: Fast testing enabled by coarse ne4pg2 resolution

Examples:
    # Run the full test suite with defaults
    $0

    # Use custom input/output directories
    $0 --input-data ~/e3sm_data/input --output ~/e3sm_data/output

    # Start an interactive session for debugging
    $0 --interactive

    # Dry run to check configuration
    $0 --dry-run

    # Use environment variables
    E3SM_INPUT_DATA=~/data/input E3SM_OUTPUT=~/data/output $0

Resolution Note:
    These tests use ne4pg2 resolution (~100km, 384 columns) for fast
    functional testing. This validates isotope physics and tracers with
    minimal resource requirements. For production-quality validation,
    use ne30pg2 or higher resolution in separate production runs.

For more information, see: ${REPO_PATH}/test_scripts/README.md

EOF
}

check_docker() {
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        log_error "Docker daemon is not running"
        exit 1
    fi
    
    log_success "Docker is available"
}

check_directories() {
    log_info "Checking directories..."
    
    # Check repository
    if [ ! -d "${REPO_PATH}" ]; then
        log_error "Repository not found: ${REPO_PATH}"
        log_error "Use --repo PATH to specify the correct location"
        exit 1
    fi
    
    if [ ! -f "${REPO_PATH}/test_scripts/test_isotopes.sh" ]; then
        log_error "Test script not found in repository: ${REPO_PATH}/test_scripts/test_isotopes.sh"
        log_error "Make sure you're pointing to the EAMv2-wiso repository root"
        exit 1
    fi
    log_success "Repository found: ${REPO_PATH}"
    
    # Create input data directory if needed
    if [ ! -d "${INPUT_DATA_PATH}" ]; then
        log_warn "Input data directory not found: ${INPUT_DATA_PATH}"
        log_info "Creating directory..."
        
        # Try to create without sudo first
        if mkdir -p "${INPUT_DATA_PATH}" 2>/dev/null; then
            log_success "Created: ${INPUT_DATA_PATH}"
        elif command -v sudo &>/dev/null; then
            # Fall back to sudo if needed
            log_info "Trying with sudo..."
            sudo mkdir -p "${INPUT_DATA_PATH}"
            sudo chown -R $USER:$(id -gn) "${INPUT_DATA_PATH}"
            log_success "Created with sudo: ${INPUT_DATA_PATH}"
        else
            log_error "Cannot create directory: ${INPUT_DATA_PATH}"
            log_error "Please create it manually or specify a different path with --input-data"
            exit 1
        fi
    else
        log_success "Input data directory found: ${INPUT_DATA_PATH}"
    fi
    
    # Create output directory if needed
    if [ ! -d "${OUTPUT_PATH}" ]; then
        log_warn "Output directory not found: ${OUTPUT_PATH}"
        log_info "Creating directory..."
        
        # Try to create without sudo first
        if mkdir -p "${OUTPUT_PATH}" 2>/dev/null; then
            log_success "Created: ${OUTPUT_PATH}"
        elif command -v sudo &>/dev/null; then
            # Fall back to sudo if needed
            log_info "Trying with sudo..."
            sudo mkdir -p "${OUTPUT_PATH}"
            sudo chown -R $USER:$(id -gn) "${OUTPUT_PATH}"
            log_success "Created with sudo: ${OUTPUT_PATH}"
        else
            log_error "Cannot create directory: ${OUTPUT_PATH}"
            log_error "Please create it manually or specify a different path with --output"
            exit 1
        fi
    else
        log_success "Output directory found: ${OUTPUT_PATH}"
    fi
}

pull_image() {
    log_info "Checking Docker image..."
    
    if docker image inspect "${DOCKER_IMAGE}" &> /dev/null; then
        log_success "Docker image found: ${DOCKER_IMAGE}"
    else
        log_warn "Docker image not found locally: ${DOCKER_IMAGE}"
        log_info "Pulling image (this may take a few minutes)..."
        docker pull "${DOCKER_IMAGE}"
        log_success "Image pulled successfully"
    fi
}

run_tests() {
    print_header "Starting EAMv2 Water Isotope Test Suite"
    
    log_info "Running tests in Docker container..."
    log_info "This will take approximately 2-6 hours"
    log_info ""
    
    docker run -it \
        -v "${REPO_PATH}:${CONTAINER_REPO}" \
        -v "${INPUT_DATA_PATH}:${CONTAINER_INPUT}" \
        -v "${OUTPUT_PATH}:${CONTAINER_OUTPUT}" \
        -e DEBUG_COMPILE="${DEBUG_COMPILE}" \
        --hostname docker \
        "${DOCKER_IMAGE}" \
        bash "${TEST_SCRIPT}"
    
    local exit_code=$?
    
    echo ""
    if [ $exit_code -eq 0 ]; then
        log_success "Test suite completed successfully!"
    else
        log_error "Test suite failed with exit code: ${exit_code}"
    fi
    
    return $exit_code
}

run_interactive() {
    print_header "Starting Interactive Docker Session"
    
    log_info "Starting container in interactive mode..."
    log_info "To run tests manually, execute:"
    log_info "    bash ${TEST_SCRIPT}"
    log_info ""
    
    docker run -it \
        -v "${REPO_PATH}:${CONTAINER_REPO}" \
        -v "${INPUT_DATA_PATH}:${CONTAINER_INPUT}" \
        -v "${OUTPUT_PATH}:${CONTAINER_OUTPUT}" \
        -e DEBUG_COMPILE="${DEBUG_COMPILE}" \
        --hostname docker \
        "${DOCKER_IMAGE}" \
        bash
}

################################################################################
# Main
################################################################################

main() {
    local interactive=false
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--interactive)
                interactive=true
                shift
                ;;
            --input-data)
                INPUT_DATA_PATH="$2"
                shift 2
                ;;
            --output)
                OUTPUT_PATH="$2"
                shift 2
                ;;
            --repo)
                REPO_PATH="$2"
                shift 2
                ;;
            --image)
                DOCKER_IMAGE="$2"
                shift 2
                ;;
            --debug)
                DEBUG_COMPILE=true
                shift
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    print_header "EAMv2 Water Isotope Test Suite - Docker Runner"
    
    # Pre-flight checks
    check_docker
    check_directories
    
    if [ "$DRY_RUN" = true ]; then
        print_header "Dry Run - Configuration Summary"
        echo "Docker Image:    ${DOCKER_IMAGE}"
        echo "Repository:      ${REPO_PATH}"
        echo "                 → ${CONTAINER_REPO}"
        echo "Input Data:      ${INPUT_DATA_PATH}"
        echo "                 → ${CONTAINER_INPUT}"
        echo "Output:          ${OUTPUT_PATH}"
        echo "                 → ${CONTAINER_OUTPUT}"
        echo "Test Script:     ${TEST_SCRIPT}"
        echo "Debug Compile:   ${DEBUG_COMPILE}"
        echo ""
        echo "Docker command that would be run:"
        echo "docker run -it \\"
        echo "  -v \"${REPO_PATH}:${CONTAINER_REPO}\" \\"
        echo "  -v \"${INPUT_DATA_PATH}:${CONTAINER_INPUT}\" \\"
        echo "  -v \"${OUTPUT_PATH}:${CONTAINER_OUTPUT}\" \\"
        echo "  -e DEBUG_COMPILE=\"${DEBUG_COMPILE}\" \\"
        echo "  --hostname docker \\"
        echo "  \"${DOCKER_IMAGE}\" \\"
        echo "  bash \"${TEST_SCRIPT}\""
        echo ""
        log_success "Dry run complete - everything looks good!"
        exit 0
    fi
    
    pull_image
    
    echo ""
    
    # Run tests or interactive session
    if [ "$interactive" = true ]; then
        run_interactive
    else
        run_tests
        exit $?
    fi
}

main "$@"
