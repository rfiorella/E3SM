---
spec_id: 2026-05-28-extend-qv-tracer-2b-shoc
spec_type: model-e3sm
spec_version: 1
title: "Extend qv to tracer dimension - Part 2b: SHOC Process"
created: 2026-05-29T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

dependencies:
  - 2026-05-28-water-tracer-metadata-and-gate
  - 2026-05-28-extend-qv-tracer-2a-infrastructure

inputs:
  source_files:
    - specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.hpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
    - components/eamxx/src/physics/shoc/eti/shoc_main.cpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
  - components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
  - components/eamxx/src/physics/shoc/impl/shoc_assumed_pdf_impl.hpp
  - components/eamxx/src/physics/shoc/impl/shoc_pblintd_impl.hpp
  - tests/water_tracers/test_qv_shoc.cpp
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-shoc-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2b-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 -DSCREAM_TRACER_ACCESS=SUBVIEW && cmake --build build/pr2b-n1 --target shoc -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
      - claim: "SHOC process compiles with tracer-aware qv field"

  - id: unit-test-shoc-qv-access
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2b-n1 -R test_qv_shoc --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_qv_shoc.cpp
      - claim: "SHOC kernels correctly access qv slot-0"

  - id: shoc-standalone-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr2b-n1 -R shoc_stand_alone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "SHOC standalone tests pass BFB vs pre-campaign baseline"
    on_fail: "Halt. Investigate which SHOC kernel broke BFB. Check accessor pattern usage."
    resolution_notes: "If BFB fails, bisect SHOC kernels to find culprit. Verify accessor pattern from spec 2a is used correctly."

  - id: shoc-integration-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr2b-n1 -R dp_.*shoc --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "SHOC integration tests (doubly-periodic) pass BFB"
    on_fail: "Halt. Check field registration and boundary exchange for tracer dimension."
    resolution_notes: "Integration test failure suggests field manager or boundary handling issue, not kernel issue."

  - id: compile-shoc-n3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2b-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr2b-n3 --target shoc -j"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "SHOC compiles with multiple tracers"

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  machine: docker-local
  baseline_tag: pre-group1
  dp_proxy: true

---

# Purpose

Extend SHOC (planetary boundary layer turbulence scheme) to use tracer-aware qv field. This is the first physics process modification and establishes the pattern for remaining processes in specs 2c-2g.

SHOC is chosen first because:
1. It's a major qv consumer (PBL mixing directly affects water vapor)
2. Has standalone test suite for isolated BFB validation
3. Relatively contained scope (~3 kernel impl files)

# Approach

Use the accessor pattern validated in spec 2a (either EXPLICIT or SUBVIEW based on 2a's BFB results).

## Step 1: Update Field Registration

In `eamxx_shoc_process_interface.cpp`:
```cpp
// Old
add_field<Required>("qv", scalar3d_mid, kg/kg, grid);

// New (use accessor pattern from spec 2a)
auto qv_layout = grid->get_3d_tracer_layout(SCREAM_NUM_TRACERS);
add_field<Required>("qv", qv_layout, kg/kg, grid);
```

## Step 2: Update Kernel Calls

In `shoc_main_impl.hpp` and other kernel implementation files:

**If spec 2a selected EXPLICIT pattern:**
```cpp
// Old
auto qv_view = qv.get_view<Real**>();
shoc_main(qv_view, ...);  // kernel expects (col, lev)

// New
auto qv_view = qv.get_view<Real***>();  // Now (tracer, col, lev)
shoc_main(qv_view, ...);  // kernel updated to expect (tracer, col, lev)

// Inside kernel:
qv(0, icol, ilev)  // Explicit slot-0 access
```

**If spec 2a selected SUBVIEW pattern:**
```cpp
// Old
auto qv_view = qv.get_view<Real**>();
shoc_main(qv_view, ...);

// New
auto qv_view = qv.get_view<Real***>();
auto qv_bulk = get_tracer_bulk_subview(qv_view);  // Helper from 2a
shoc_main(qv_bulk, ...);  // kernel unchanged, still sees (col, lev)
```

## Step 3: Update All SHOC Kernels

Kernels to update:
- `shoc_main_impl.hpp` - main driver
- `shoc_assumed_pdf_impl.hpp` - cloud fraction calculation using qv
- `shoc_pblintd_impl.hpp` - PBL depth calculation using qv

Each kernel either:
- Gets explicit `(0, icol, ilev)` indexing, OR
- Receives a 2D subview and remains unchanged

## Step 4: Testing

**Unit test** (`test_qv_shoc.cpp`):
- Verify SHOC can read tracer qv field
- Verify slot-0 access matches expected values
- Simple smoke test with known input

**Standalone tests** (existing):
- SHOC's existing standalone test suite runs with SCREAM_NUM_TRACERS=1
- Must pass BFB vs pre-campaign baseline

**Integration tests** (existing):
- Doubly-periodic test cases with SHOC enabled
- Must pass BFB vs pre-campaign baseline

# Implementation Pattern

This spec establishes the pattern for specs 2c-2g:
1. Update field registration to use tracer layout
2. Update kernel calls to use accessor pattern from spec 2a
3. Update kernel implementations consistently
4. Validate with standalone + integration tests
5. Enforce BFB requirement

# Resolved Decisions

## Accessor pattern
Use the pattern that passed BFB in spec 2a. The CMake variable `{{PATTERN}}` in success_criteria will be replaced with either `EXPLICIT` or `SUBVIEW` based on spec 2a results.

## Kernel modification strategy
All SHOC kernels updated atomically in this spec. Not incremental per-kernel because:
1. SHOC kernels are tightly coupled (main → assumed_pdf → pblintd)
2. Standalone tests require full SHOC path to be functional
3. Small scope (~3 impl files) makes atomic update low-risk

## Testing strategy
Use existing SHOC test infrastructure. No new test cases needed beyond unit test for field access. BFB requirement enforced via existing baseline comparison.

# Ask Before

1. If SHOC standalone tests fail BFB and investigation reveals architectural issue (not implementation bug)
2. If accessor pattern from spec 2a doesn't apply cleanly to SHOC (suggests pattern needs refinement)
3. If SHOC requires special handling different from pattern established in spec 2a

# Out of Scope

- Other physics processes (P3, RRTMGP, ZM) - deferred to specs 2c-2f
- Dynamics (HOMME) - deferred to spec 2f
- Surface coupling - deferred to spec 2g
- Field manager changes beyond field registration (already handled in 2a)
- Multi-tracer physics (ntracers > 1 actual physics) - deferred to future campaigns

# Notes

## Why SHOC First?

SHOC is strategically chosen as the first process because:
1. **Standalone test suite**: Can validate BFB in isolation before touching other processes
2. **Moderate complexity**: Not trivial (like surface coupling) but not massive (like P3)
3. **Clear qv usage**: SHOC uses qv for PBL diagnostics and mixing, with well-defined input/output
4. **Pattern validation**: Success here proves the infrastructure from spec 2a works in a real process

## Success Criteria Dependencies

The `compile-shoc-n1` criterion uses `{{PATTERN}}` placeholder which should be replaced with:
- `EXPLICIT` if spec 2a's explicit-indexing-bfb passed and subview-bfb failed
- `SUBVIEW` if spec 2a's subview-bfb passed

This is documented in spec 2a's resolved_decisions and should be carried forward to all subsequent specs.

## BFB Strategy

Two-level BFB validation:
1. **Standalone**: Isolates SHOC-specific issues
2. **Integration**: Catches field manager / boundary exchange issues

Both must pass. Standalone failure → kernel bug. Integration failure → infrastructure bug.
