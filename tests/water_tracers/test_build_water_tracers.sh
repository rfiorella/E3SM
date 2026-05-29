#!/bin/bash
# Test script to validate water tracer metadata headers compile
# Usage: bash tests/water_tracers/test_build_water_tracers.sh

set -e

echo "=========================================="
echo "Water Tracer Metadata Compilation Test"
echo "=========================================="
echo

# Navigate to water_tracers directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WATER_TRACERS_DIR="$REPO_ROOT/components/eamxx/src/physics/water_tracers"

echo "Repository root: $REPO_ROOT"
echo "Water tracers directory: $WATER_TRACERS_DIR"
echo

# Create a minimal test file that includes all headers
TEST_FILE=$(mktemp /tmp/test_water_tracer_headers.XXXXXX.cpp)
trap "rm -f $TEST_FILE" EXIT

cat > $TEST_FILE << 'EOF'
// Minimal test that all water tracer headers compile
#include <iostream>

// Include water tracer headers
#include "components/eamxx/src/physics/water_tracers/water_tracer_metadata.hpp"
#include "components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp"

int main() {
  using namespace scream::water_tracers;

  // Test metadata struct construction
  WaterTracerMetadata bulk("bulk_H2O", "Bulk Water", WaterTracerKind::bulk);
  WaterTracerMetadata test("test_tracer", "Test Tracer", WaterTracerKind::evaporation);

  std::cout << "WaterTracerMetadata compilation: PASS\n";

  // Test registry singleton
  auto& registry = WaterTracerRegistry::instance();
  std::cout << "WaterTracerRegistry instantiation: PASS\n";

  // Test registry operations
  int count = registry.count();
  std::cout << "Registry tracer count: " << count << " (should be 1 for bulk_H2O)\n";

  if (count == 1) {
    auto bulk_meta = registry.get(0);
    std::cout << "Bulk tracer name: " << bulk_meta.name << "\n";
    std::cout << "Bulk tracer longname: " << bulk_meta.longname << "\n";
  }

  // Test adding a tracer
  WaterTracerMetadata new_tracer("test", "Test", WaterTracerKind::evaporation);
  registry.register_tracer(new_tracer);

  count = registry.count();
  std::cout << "Registry tracer count after registration: " << count << " (should be 2)\n";

  std::cout << "\nAll water tracer header tests: PASS\n";
  return 0;
}
EOF

echo "Created test file: $TEST_FILE"
echo

# Compile test
echo "Compiling test file..."
cd "$REPO_ROOT"

if g++ -std=c++17 -Wall -Wextra -I. -o /tmp/test_water_tracer_exe $TEST_FILE; then
  echo "Compilation: SUCCESS"
  echo

  # Run test
  echo "Running test executable..."
  if /tmp/test_water_tracer_exe; then
    echo
    echo "=========================================="
    echo "All tests PASSED"
    echo "=========================================="
    rm -f /tmp/test_water_tracer_exe
    exit 0
  else
    echo "ERROR: Test executable failed"
    rm -f /tmp/test_water_tracer_exe
    exit 1
  fi
else
  echo "ERROR: Compilation failed"
  exit 1
fi
