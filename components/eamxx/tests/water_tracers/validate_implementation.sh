#!/bin/bash
#
# Validation script for spec 2026-05-28-extend-qv-tracer-2a-infrastructure
#
# This script must be run on a properly configured E3SM machine with:
# - CMake 3.18+
# - Kokkos (via EKAT)
# - MPI
# - Build environment loaded (use ./scripts/test-all-eamxx --preserve-env -m MACHINE)
#
# Usage:
#   cd components/eamxx
#   ./tests/water_tracers/validate_implementation.sh

set -e

echo "=== Spec 2a Infrastructure Validation ==="
echo ""

EAMXX_ROOT=$(pwd)
BUILD_DIR="$EAMXX_ROOT/build/spec2a-validation"

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf "$BUILD_DIR"

echo ""
echo "=== Success Criterion: compile-infrastructure ==="
echo "Configuring with SCREAM_NUM_TRACERS=1..."
cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DSCREAM_NUM_TRACERS=1

echo ""
echo "Building field infrastructure..."
cmake --build "$BUILD_DIR" --target scream_share -j4

if [ $? -eq 0 ]; then
  echo "SUCCESS: Infrastructure compiles with TRACER FieldTag"
else
  echo "FAIL: Infrastructure compilation failed"
  exit 1
fi

echo ""
echo "=== Success Criterion: unit-test-tracer-layout ==="
ctest --test-dir "$BUILD_DIR" -R test_tracer_field_infrastructure --output-on-failure

if [ $? -eq 0 ]; then
  echo "SUCCESS: Tracer layout tests pass"
else
  echo "FAIL: Tracer layout tests failed"
  exit 1
fi

echo ""
echo "=== Success Criterion: explicit-indexing-bfb (BLOCKING GATE) ==="
echo "Building with SCREAM_TRACER_ACCESS=EXPLICIT..."

BUILD_EXPLICIT="$EAMXX_ROOT/build/spec2a-explicit"
rm -rf "$BUILD_EXPLICIT"

cmake -S . -B "$BUILD_EXPLICIT" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSCREAM_NUM_TRACERS=1 \
  -DSCREAM_TRACER_ACCESS=EXPLICIT

cmake --build "$BUILD_EXPLICIT" -j4

ctest --test-dir "$BUILD_EXPLICIT" -R test_tracer_access_patterns --output-on-failure

if [ $? -eq 0 ]; then
  echo "SUCCESS: Explicit indexing passes BFB"
  echo "EXPLICIT" > /tmp/spec2a_bfb_winner.txt
else
  echo "FAIL: Explicit indexing BFB failure - this is a BLOCKING gate"
  exit 1
fi

echo ""
echo "=== Success Criterion: subview-bfb (ADVISORY GATE) ==="
echo "Building with SCREAM_TRACER_ACCESS=SUBVIEW..."

BUILD_SUBVIEW="$EAMXX_ROOT/build/spec2a-subview"
rm -rf "$BUILD_SUBVIEW"

cmake -S . -B "$BUILD_SUBVIEW" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSCREAM_NUM_TRACERS=1 \
  -DSCREAM_TRACER_ACCESS=SUBVIEW

cmake --build "$BUILD_SUBVIEW" -j4

ctest --test-dir "$BUILD_SUBVIEW" -R test_tracer_access_patterns --output-on-failure

if [ $? -eq 0 ]; then
  echo "SUCCESS: Subview accessor passes BFB"
  echo "SUBVIEW" > /tmp/spec2a_bfb_winner.txt
  echo ""
  echo "=== DECISION: Use SUBVIEW accessor for specs 2b-2g ==="
else
  echo "ADVISORY GATE: Subview BFB failure - falling back to EXPLICIT indexing"
  echo ""
  echo "=== DECISION: Use EXPLICIT indexing for specs 2b-2g ==="
fi

echo ""
echo "=== All validation complete ==="
echo "Recommended accessor pattern: $(cat /tmp/spec2a_bfb_winner.txt)"
