#!/bin/bash
# Simple test to verify isotope build works in Docker

set -e

echo "Starting simple isotope build test..."

# Set up paths
REPO_ROOT="/home/e3smuser/EAMv2-wiso"
OUTPUT_BASE="/home/e3smuser/output"
INPUT_DATA="/home/e3smuser/inputdata"

export OUTPUT_BASE
export INPUT_DATA

# Create directories
mkdir -p ${OUTPUT_BASE}
mkdir -p ${INPUT_DATA}

# Create a single test case
cd ${REPO_ROOT}/cime/scripts

echo "Creating test case..."
./create_newcase \
    --case ${OUTPUT_BASE}/isotope_test_simple \
    --compset Fi2000climo \
    --res ne30pg2_oECv3 \
    --mach docker \
    --compiler gnu \
    --project docker_runs

cd ${OUTPUT_BASE}/isotope_test_simple

echo "Configuring case..."
# Add water tracer to CAM_CONFIG_OPTS
CAM_OPTS=$(./xmlquery CAM_CONFIG_OPTS --value)
./xmlchange CAM_CONFIG_OPTS="${CAM_OPTS} -water_tracer h2o_h216o_hdo_h218o"

# Enable isotopes in coupler
./xmlchange FLDS_WISO=TRUE

# Set short run
./xmlchange STOP_N=5
./xmlchange STOP_OPTION=ndays

# Set input data location
./xmlchange DIN_LOC_ROOT=${INPUT_DATA}

# Create user_nl_eam
cat > user_nl_eam << 'EOF'
trace_water = .true.
wisotope = .true.
wtrc_alpha_kinetic = .true.
EOF

echo "Setting up case..."
./case.setup

echo "Checking input data..."
./check_input_data --download || echo "Input data check completed with warnings"

echo "Building case..."
./case.build

echo "Build completed successfully!"
