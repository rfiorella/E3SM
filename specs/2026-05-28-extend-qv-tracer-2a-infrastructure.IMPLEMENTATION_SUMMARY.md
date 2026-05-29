# Implementation Summary: Spec 2a - Tracer Field Infrastructure

**Status**: IMPLEMENTATION COMPLETE - AWAITING MACHINE ENVIRONMENT FOR TESTING  
**Date**: 2026-05-29  
**Spec ID**: 2026-05-28-extend-qv-tracer-2a-infrastructure

## Overview

All code deliverables for tracer field infrastructure have been implemented. The implementation establishes:

1. TRACER as a first-class FieldTag
2. get_3d_tracer_layout() method in grid API
3. Two accessor patterns (EXPLICIT and SUBVIEW) with helpers
4. Comprehensive unit tests for BFB validation

## Deliverables Implemented

### Phase 1: TRACER FieldTag

**File**: `components/eamxx/src/share/field/field_tag.hpp`
- Added `Tracer` to `FieldTag` enum
- Added `TRACER` constant in `ShortFieldTagsNames` namespace
- Updated `e2str()` helper to return "tracer" for `FieldTag::Tracer`

### Phase 2: Grid Layout Methods

**Files Modified**:
- `components/eamxx/src/share/grid/abstract_grid.hpp`
  - Added pure virtual method: `get_3d_tracer_layout(const int ntracers, const std::string& name = "")`
  
- `components/eamxx/src/share/grid/point_grid.hpp`
  - Added override declaration
  
- `components/eamxx/src/share/grid/point_grid.cpp`
  - Implemented method returning layout: `{TRACER, COL, LEV}, {ntracers, ncols, nlevs}`
  
- `components/eamxx/src/share/grid/se_grid.hpp`
  - Added override declaration
  
- `components/eamxx/src/share/grid/se_grid.cpp`
  - Implemented method returning layout: `{TRACER, EL, GP, GP, LEV}, {ntracers, nelem, ngp, ngp, nlevs}`

### Phase 3: Accessor Pattern Helpers

**File**: `components/eamxx/src/share/field/field_tracer_access.hpp` (NEW)

Implements two accessor patterns:

1. **EXPLICIT INDEXING**:
   ```cpp
   auto accessor = get_tracer_bulk_explicit(qv_view);
   Real value = accessor(icol, ilev);  // internally: qv_view(0, icol, ilev)
   ```

2. **SUBVIEW**:
   ```cpp
   auto qv_bulk = get_tracer_bulk_subview(qv_view);
   Real value = qv_bulk(icol, ilev);  // 2D subview of slot 0
   ```

3. **CONDITIONAL ACCESSOR** (compile-time selection):
   ```cpp
   auto qv_bulk = get_tracer_bulk(qv_view);
   // Pattern selected by SCREAM_TRACER_ACCESS=[EXPLICIT|SUBVIEW]
   ```

### Phase 4: Unit Tests

**Files Created**:

1. `tests/water_tracers/CMakeLists.txt`
   - Registers two test executables using EkatCreateUnitTest

2. `tests/water_tracers/test_tracer_field_infrastructure.cpp`
   - Tests TRACER FieldTag recognition
   - Validates PointGrid tracer layouts (bulk-only and multi-tracer)
   - Validates SEGrid tracer layouts
   - Verifies dimension ordering consistency

3. `tests/water_tracers/test_tracer_access_patterns.cpp`
   - BFB validation: scalar baseline vs explicit indexing
   - BFB validation: scalar baseline vs subview accessor
   - Tests representative atmospheric kernel operations
   - Validates accessor helper functions

4. `tests/CMakeLists.txt`
   - Added `add_subdirectory(water_tracers)`

5. `tests/water_tracers/validate_implementation.sh` (NEW)
   - End-to-end validation script
   - Runs all success criteria in sequence
   - Documents recommended accessor pattern

## Architecture Decisions

### Tracer Dimension Ordering

**Point Grid**: `(tracer, col, lev)`
- Tracer dimension first for cache-friendly access to bulk water
- Column dimension partitioned across MPI ranks
- Level dimension innermost for vertical operations

**SE Grid**: `(tracer, elem, gp, gp, lev)`
- Tracer dimension first (consistent with Point Grid)
- Element dimension partitioned across MPI ranks
- GLL point dimensions middle
- Level dimension innermost

### Accessor Pattern Design

Both patterns access **slot 0** (bulk water) from tracer-dimension arrays:

- **EXPLICIT**: Direct 3D indexing `field(0, icol, ilev)`
  - Pros: Simple, no view creation overhead
  - Cons: Slightly more verbose in physics code
  
- **SUBVIEW**: Kokkos subview collapsing tracer dimension
  - Pros: Cleaner physics code, potential GPU optimization
  - Cons: View creation overhead (may impact BFB)

**Decision Point**: BFB tests in success criteria will determine which pattern to use in specs 2b-2g.

## Testing Requirements

### Success Criteria Mapping

1. **compile-infrastructure**: Verify TRACER FieldTag and layouts compile
   - Command: `cmake -B build -DSCREAM_NUM_TRACERS=1 && cmake --build build --target scream_share`

2. **unit-test-tracer-layout**: Validate layout correctness
   - Command: `ctest -R test_tracer_field_infrastructure`

3. **explicit-indexing-bfb** (BLOCKING GATE):
   - Command: `cmake -DSCREAM_TRACER_ACCESS=EXPLICIT && ctest -R test_tracer_access_patterns`
   - If fails: Fundamental layout issue - must fix before proceeding

4. **subview-bfb** (ADVISORY GATE):
   - Command: `cmake -DSCREAM_TRACER_ACCESS=SUBVIEW && ctest -R test_tracer_access_patterns`
   - If fails: Fall back to EXPLICIT for specs 2b-2g

### Environment Requirements

Tests require:
- CMake 3.18+
- Kokkos (provided via EKAT external)
- Catch2 (provided by EAMxx test infrastructure)
- C++17 compiler
- MPI (for grid partitioning tests)

**Recommended**: Use `./scripts/test-all-eamxx -m MACHINE` workflow or docker container.

## Files Not Modified (Out of Scope)

Per spec, the following are intentionally NOT modified:
- Any physics process code (deferred to specs 2b-2g)
- Field manager registration for qv (deferred to spec 2b)
- water_tracer_registry integration (deferred to spec 2b)
- FieldLayout type system (Scalar3D, Vector3D, etc.) - may need future extension

## Next Steps

1. **Run validation script** on configured E3SM machine:
   ```bash
   cd components/eamxx
   ./tests/water_tracers/validate_implementation.sh
   ```

2. **Document BFB results** in progress ledger:
   - Which accessor pattern passed?
   - Performance overhead measurements
   - Any unexpected failures

3. **Update campaign resolved_decisions** with:
   - BFB validation outcome
   - Chosen accessor pattern for specs 2b-2g
   - Any fallback strategies if both patterns fail

4. **Proceed to spec 2b** (qv field registration) using validated accessor pattern

## Known Limitations

1. **Grid Type Coverage**: Only PointGrid and SEGrid implemented
   - Other grid types (if any) will need tracer layout methods added
   - Unlikely to be needed for EAMxx water tracers

2. **LayoutType System**: TRACER dimension not yet integrated into LayoutType enum
   - Current system: Scalar3D, Vector3D, Tensor3D
   - Tracer fields don't fit this categorization
   - May need future extension if FieldLayout::type() queries are needed

3. **CMake Integration**: SCREAM_TRACER_ACCESS flag not yet in main CMakeLists.txt
   - Tests assume flag is passed manually
   - Should be added to machine files or project defaults after BFB validation

## Risk Mitigation

- **BFB Failure**: Both accessor patterns tested independently
  - If EXPLICIT fails: Layout design issue - halt and redesign
  - If SUBVIEW fails: Known risk - fall back to EXPLICIT
  - If both fail: Consult user per spec's "Ask Before" clause

- **Performance Overhead**: Advisory gates track >2% overhead
  - Will not block spec completion
  - Document if observed, consider optimizations post-campaign

## References

- Spec: `specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md`
- Progress: `specs/2026-05-28-extend-qv-tracer-2a-infrastructure.progress.md`
- Campaign: `campaigns/wiso.campaign.md`
- Dependency: Spec 1 (water-tracer-metadata-and-gate)
