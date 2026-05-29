---
spec_id: 2026-05-28-extend-qv-tracer-2a-infrastructure
spec_type: model-e3sm
spec_version: 1
title: "Extend qv to tracer dimension - Part 2a: Field Infrastructure"
created: 2026-05-29T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

dependencies:
  - 2026-05-28-water-tracer-metadata-and-gate

inputs:
  source_files:
    - specs/2026-05-28-water-tracer-metadata-and-gate.md
    - components/eamxx/src/share/field/field_layout.hpp
    - components/eamxx/src/share/field/field_layout.cpp
    - components/eamxx/src/share/grid/abstract_grid.hpp
    - components/eamxx/src/share/grid/grids_manager.hpp
    - components/eamxx/src/share/field/field_manager.hpp
    - components/eamxx/src/share/field/field_manager.cpp
  data: []
  baseline: null

deliverables:
  - components/eamxx/src/share/field/field_layout.hpp
  - components/eamxx/src/share/field/field_layout.cpp
  - components/eamxx/src/share/grid/abstract_grid.hpp
  - components/eamxx/src/share/grid/abstract_grid.cpp
  - components/eamxx/src/share/grid/point_grid.hpp
  - components/eamxx/src/share/grid/point_grid.cpp
  - components/eamxx/src/share/field/field_tracer_access.hpp
  - tests/water_tracers/test_tracer_field_infrastructure.cpp
  - tests/water_tracers/test_tracer_access_patterns.cpp
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-infrastructure
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2a -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 && cmake --build build/pr2a --target field_layout field_manager -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/share/field/field_layout.hpp
      - deliverable: components/eamxx/src/share/grid/abstract_grid.hpp
      - claim: "TRACER FieldTag and get_3d_tracer_layout() compile without errors"

  - id: unit-test-tracer-layout
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2a -R test_tracer_field_infrastructure --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_tracer_field_infrastructure.cpp
      - claim: "FieldLayout with TRACER dimension creates correct 3D layout (tracer, col, lev)"

  - id: explicit-indexing-bfb
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2a-explicit -DCMAKE_BUILD_TYPE=Release -DSCREAM_NUM_TRACERS=1 -DSCREAM_TRACER_ACCESS=EXPLICIT && cmake --build build/pr2a-explicit -j && ctest --test-dir build/pr2a-explicit -R test_tracer_access_patterns --output-on-failure"
    expect: exit_zero
    phase: testing
    gate: blocking
    on_fail: "Document BFB failure and halt. Investigate field layout or indexing arithmetic."
    resolution_notes: "If explicit indexing fails BFB, the fundamental layout is wrong. Must fix before subview testing."
    verifies:
      - claim: "Explicit slot-0 indexing qv(0, icol, ilev) produces BFB-identical results vs scalar baseline"

  - id: subview-bfb
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2a-subview -DCMAKE_BUILD_TYPE=Release -DSCREAM_NUM_TRACERS=1 -DSCREAM_TRACER_ACCESS=SUBVIEW && cmake --build build/pr2a-subview -j && ctest --test-dir build/pr2a-subview -R test_tracer_access_patterns --output-on-failure"
    expect: exit_zero
    phase: testing
    gate: advisory
    on_fail: "Document subview BFB failure. Fall back to explicit indexing for specs 2b-2g. Update resolved_decisions."
    resolution_notes: "If subview fails BFB but explicit passed, use EXPLICIT for remaining specs. Subview can be revisited post-campaign."
    verifies:
      - claim: "Kokkos subview accessor produces BFB-identical results vs scalar baseline"

  - id: performance-overhead-explicit
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark_field_access.sh --pattern=EXPLICIT --tracers=1 --output=bench-explicit.json && python3 scripts/compare_benchmarks.py bench-baseline.json bench-explicit.json --threshold=2.0"
    expect: exit_zero
    phase: testing
    gate: advisory
    verifies:
      - claim: "Explicit indexing overhead < 2% vs scalar baseline"

  - id: performance-overhead-subview
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark_field_access.sh --pattern=SUBVIEW --tracers=1 --output=bench-subview.json && python3 scripts/compare_benchmarks.py bench-baseline.json bench-subview.json --threshold=2.0"
    expect: exit_zero
    phase: testing
    gate: advisory
    verifies:
      - claim: "Subview accessor overhead < 2% vs scalar baseline"

model_specific:
  validation_tier: 1
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  machine: docker-local
  baseline_tag: pre-group1
  test_suite: null

---

# Purpose

Extend EAMxx field infrastructure to support water tracer dimensions. This is the foundational infrastructure for extending `qv` and all other water species from `(col, lev)` to `(tracer, col, lev)` layout.

This spec establishes:
1. `TRACER` as a first-class `FieldTag`
2. `get_3d_tracer_layout(ncols, nlevs, ntracers)` method in grid API
3. Two accessor patterns (explicit indexing and Kokkos subview) with BFB validation for both
4. The accessor pattern decision for specs 2b-2g (based on which passes BFB)

**This spec does NOT modify any physics process code**. It only establishes the infrastructure.

# Approach

## Phase 1: Add TRACER FieldTag

Extend `FieldTag` enum in `field_layout.hpp`:
```cpp
enum class FieldTag {
  // ... existing tags ...
  TRACER,  // Water tracer constituent dimension (bulk + isotopes)
};
```

Update `to_string()` and `from_string()` helpers.

## Phase 2: Implement get_3d_tracer_layout()

Add method to `AbstractGrid` and concrete implementations:
```cpp
// abstract_grid.hpp
virtual FieldLayout get_3d_tracer_layout(
  const int ntracers,
  const std::string& name = "") const = 0;

// point_grid.cpp (and other grid types)
FieldLayout PointGrid::get_3d_tracer_layout(
  const int ntracers, 
  const std::string& name) const {
  using namespace ShortFieldTagsNames;
  return FieldLayout({TRACER, COL, LEV}, {ntracers, m_num_local_cols, m_num_levs}, name);
}
```

## Phase 3: Accessor Pattern Helpers

Create `field_tracer_access.hpp` with both patterns:
```cpp
// Explicit indexing helper (always use slot 0)
template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk_explicit(const ViewT& field) {
  // Returns functor that does field(0, icol, ilev) access
  // Preserves 3D view semantics but enforces slot-0
}

// Subview helper
template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk_subview(const ViewT& field) {
  return Kokkos::subview(field, 0, Kokkos::ALL, Kokkos::ALL);
  // Returns 2D view into slot-0
}
```

Controlled by CMake: `-DSCREAM_TRACER_ACCESS=[EXPLICIT|SUBVIEW]`

## Phase 4: Unit Tests

**test_tracer_field_infrastructure.cpp**:
- Verify `FieldLayout` with `TRACER` dimension creates correct shape
- Verify `get_3d_tracer_layout(1, ...)` for bulk-only case
- Verify `get_3d_tracer_layout(3, ...)` for multi-tracer case

**test_tracer_access_patterns.cpp**:
- Create a simple kernel that reads/writes a tracer field
- Run with both EXPLICIT and SUBVIEW patterns
- Compare output against known scalar baseline (BFB check)
- Measure performance overhead

# Implementation Pattern

This spec establishes infrastructure only. No physics process modifications.

Specs 2b-2g will use the accessor pattern that passes BFB in this spec:
- If both pass BFB → use SUBVIEW (cleaner, slight perf win)
- If only EXPLICIT passes → use EXPLICIT (safer)
- If neither passes → HALT campaign, architectural redesign needed

# Resolved Decisions

## Decision 1: Compile-time tracer count
Use compile-time constant `SCREAM_NUM_TRACERS` set at CMake configure time.

**Rationale**: PR 1's BFB gate passed using compile-time constant. Simpler implementation, easier to maintain BFB guarantee across 40+ files. Runtime flexibility not required for initial isotope implementation.

**Implications**: Changing tracer count requires reconfigure + rebuild. Acceptable for research use case.

## Decision 2: New grid layout method
Add new `get_3d_tracer_layout()` method instead of repurposing `get_3d_vector_layout()`.

**Rationale**: Tracer dimension has different semantics than vector dimension (CMP, coordinate components). Clear separation aids documentation and future multi-dimensional fields (e.g., tracer + vertical levels per aerosol mode).

**Implications**: One new virtual method in `AbstractGrid`, implemented by all grid types (PointGrid, SEGrid, etc.).

## Decision 3: Hybrid accessor pattern validation
Implement both explicit indexing and Kokkos subview patterns. Test both for BFB in this spec. Use the pattern that passes BFB for specs 2b-2g.

**Rationale**: Subviews offer cleaner code and potentially better GPU performance, but BFB is non-negotiable. Testing both in isolated infrastructure spec minimizes risk.

**Implications**: 
- If subview fails BFB: Fall back to explicit indexing, document why in campaign notes
- If both pass: Use subview by default (performance + readability)
- If neither passes: Fundamental layout issue, halt campaign

# Ask Before

1. If both accessor patterns fail BFB in this spec, halt and consult before proceeding
2. If performance overhead exceeds 5% for either pattern, halt and consult
3. If grid implementations beyond PointGrid fail (SEGrid, etc.), document and decide whether to expand scope

# Out of Scope

- Any physics process modifications (deferred to specs 2b-2g)
- Field manager registration changes for qv (deferred to spec 2b)
- Multi-tracer testing (ntracers > 1) beyond unit tests
- Integration with water_tracer_registry from PR 1 (deferred to spec 2b)

# Notes

## Test Execution

Run unit tests on docker-local (CPU). BFB tests use Release build. Performance tests use optimized build with timing instrumentation.

## Success Gate Strategy

- `explicit-indexing-bfb` is **blocking**: If this fails, layout is fundamentally wrong
- `subview-bfb` is **advisory**: If this fails, fall back to explicit indexing
- Performance tests are **advisory**: Document if >2% overhead, but don't block

## Baseline for BFB Comparison

Use scalar `qv(col, lev)` baseline from pre-campaign tag. The test will:
1. Run simple kernel with scalar layout
2. Run same kernel with tracer layout + accessor pattern
3. Compare outputs bit-for-bit

## Files Modified But Not Listed as Deliverables

- CMakeLists.txt files for test registration
- Potentially grid_manager.cpp if factory needs updates
- Test utility headers (field_test_utils.hpp, etc.)

These are implementation details, not primary deliverables.
