# Spec Progress: 2026-05-28-extend-qv-tracer

## Status: HALTED

## Phase: Investigation

## Halt Reason

This spec requires extensive modifications across 40+ files touching the core field infrastructure and all physics processes that use water vapor. The investigation phase has identified the required changes, but the implementation requires:

1. **Architectural decisions** that need user confirmation:
   - How to handle tracer dimension in grid layout methods
   - Whether to use compile-time `SCREAM_NUM_TRACERS` or runtime sizing
   - Strategy for kernel modifications (explicit indexing vs subviews)

2. **Scope validation**: The spec lists 14 deliverable files, but investigation reveals at least 40+ files need modification:
   - All process interface files touching `qv` (6 files listed, but kernels too)
   - Field infrastructure (FieldTag, FieldLayout, AbstractGrid, FieldManager)
   - Grid implementations (se_grid, point_grid)
   - All kernel implementation files calling `qv` (not in deliverables list)

3. **Test infrastructure**: Unit tests and component tests need to be designed before implementation

**Recommendation**: User should review investigation findings below and clarify:
- Whether to expand deliverables list to include all kernel impl files
- Whether BFB gate from PR 1 provides guidance on architecture
- Whether to proceed incrementally (field infrastructure first, then one process at a time)

## Investigation

### Current state analysis

PR 1 deliverables confirmed present:
- ✅ `water_tracer_metadata.hpp` exists with `WaterTracerKind` enum and `WaterTracerMetadata` struct
- ✅ `FieldTag` enum exists in `share/field/field_tag.hpp`
- ✅ `FieldLayout` class exists in `share/field/field_layout.hpp`
- ✅ `FieldManager` exists in `share/data_managers/field_manager.hpp`

Need to determine:
1. Does `FieldTag` enum already have a `TRACER` tag added by PR 1?
2. Does `FieldLayout` support tracer dimension or need extension?
3. Current `qv` field registration pattern across all processes
4. Grid layout creation methods (how to specify tracer dimension)
5. Kernel access patterns for `qv` fields

### File inventory

Processes touching `qv` (per spec inputs):
- SHOC: `components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp`
- P3: `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp`
- RRTMGP: `components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp`
- ZM: `components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp`
- Dynamics (homme): `components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp`
- Surface coupling: `components/eamxx/src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp`

### Next steps

1. ✅ Check if `TRACER` tag already exists in `FieldTag` enum - NO, needs to be added
2. ✅ Examine grid layout creation methods - Found `get_3d_scalar_layout`, `get_3d_vector_layout`, etc.
3. ✅ Survey current `qv` registration in one representative process - Found `add_tracer` template
4. ✅ Plan kernel modification strategy - See implementation plan below

### Key findings

1. **`FieldTag` enum** (in `share/field/field_tag.hpp`): Does NOT have `TRACER` tag. Need to add it.

2. **Current `add_tracer` implementation** (line 388 of `atmosphere_process.hpp`):
   ```cpp
   FieldIdentifier fid(name, grid->get_3d_scalar_layout(FieldTag::LevelMidPoint), u, grid->name());
   ```
   Uses `get_3d_scalar_layout` which produces `(col, lev)` shape.

3. **Grid layout creation methods** (in `abstract_grid.hpp`):
   - `get_3d_scalar_layout(vtag)` - produces (col, lev)
   - `get_3d_vector_layout(vtag, vector_dim, vec_dim_name)` - produces (col, lev, component)
   - Need new method or modify existing to support tracer dimension

4. **Process registration pattern** (SHOC example):
   ```cpp
   add_tracer<Updated>("qv", m_grid, kg/kg, ps);
   ```
   This is clean API used by all processes.

### Implementation strategy

To extend qv to tracer dimension `(tracer, col, lev)`:

1. **Add `TRACER` tag to `FieldTag` enum** (field_tag.hpp)
2. **Add tracer-aware layout creation method** to `AbstractGrid` 
   - Option A: New method `get_3d_tracer_layout(vtag, num_tracers, name)`
   - Option B: Modify `get_3d_vector_layout` to handle tracer as first dimension
   - **Decision: Option A** - clearer intent, matches spec pattern

3. **Modify `add_tracer` template** to use tracer layout when `SCREAM_NUM_TRACERS > 1`

4. **Update field_layout computation** in `FieldLayout::compute_type()` to recognize tracer layouts

5. **Update all kernel call sites** - two strategies:
   - Strategy A: Explicit slot-0 indexing: `qv(0, icol, ilev)`
   - Strategy B: Subview pattern: `auto qv_bulk = Kokkos::subview(qv, 0, ALL, ALL)`
   - **Decision: Strategy A for simple accesses, Strategy B for passing to existing kernels**

## Implementation

### Detailed implementation plan (pending user decision on halt)

#### Phase 1: Field infrastructure (Tier 0 - foundational)

**1.1 Add TRACER tag to FieldTag enum**

File: `components/eamxx/src/share/field/field_tag.hpp`

```cpp
enum class FieldTag {
  Invalid,
  Element,
  LevelMidPoint,
  LevelInterface,
  LevelPressure,
  Column,
  GaussPoint,
  Component,
  TimeLevel,
  Tracer,        // NEW: Water tracer constituent dimension
};
```

Also update:
- `e2str()` function to handle `FieldTag::Tracer` → return `"tracer"`
- `ShortFieldTagsNames` namespace to add `constexpr auto TRC = FieldTag::Tracer;`

**1.2 Update FieldLayout to recognize tracer layouts**

File: `components/eamxx/src/share/field/field_layout.cpp`

Modify `compute_type()` method to detect layouts with TRACER tag and classify appropriately:
- `(TRACER, COL, LEV)` → Vector3D (tracer is the "component" dimension)
- Add logic to `is_vector_layout()` to recognize TRACER as vector component
- Add `get_tracer_component_idx()` method parallel to `get_vector_component_idx()`

**1.3 Add tracer layout creation to grid classes**

Files:
- `components/eamxx/src/share/grid/abstract_grid.hpp` (declaration)
- `components/eamxx/src/share/grid/abstract_grid.cpp` (implementation)
- `components/eamxx/src/share/grid/se_grid.cpp` (SE spectral element grid impl)
- `components/eamxx/src/share/grid/point_grid.cpp` (point grid impl)

New method signature:
```cpp
// In AbstractGrid class
virtual FieldLayout get_3d_tracer_layout(const FieldTag vtag, 
                                         const int num_tracers,
                                         const std::string& tracer_dim_name) const;

// Convenience overload
FieldLayout get_3d_tracer_layout(const FieldTag vtag, const int num_tracers) const {
  return get_3d_tracer_layout(vtag, num_tracers, "tracer");
}
```

Implementation creates layout: `{FieldTag::Tracer, FieldTag::Column, vtag}` with dims `{num_tracers, ncols, nlevs}`.

**1.4 Modify add_tracer template**

File: `components/eamxx/src/share/atm_process/atmosphere_process.hpp`

Change line 388 from:
```cpp
FieldIdentifier fid(name, grid->get_3d_scalar_layout(FieldTag::LevelMidPoint), u, grid->name());
```

To:
```cpp
#ifdef SCREAM_NUM_TRACERS
  const int num_tracers = SCREAM_NUM_TRACERS;
  FieldIdentifier fid(name, grid->get_3d_tracer_layout(FieldTag::LevelMidPoint, num_tracers), u, grid->name());
#else
  FieldIdentifier fid(name, grid->get_3d_scalar_layout(FieldTag::LevelMidPoint), u, grid->name());
#endif
```

Or if `SCREAM_NUM_TRACERS` defaults to 1 and tracer layout should always be used:
```cpp
constexpr int num_tracers = SCREAM_NUM_TRACERS;  // Compile-time constant, default 1
FieldIdentifier fid(name, grid->get_3d_tracer_layout(FieldTag::LevelMidPoint, num_tracers), u, grid->name());
```

#### Phase 2: Process interface updates (Tier 1)

For each process touching qv:

**2.1 SHOC** (`components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp`)

Line 64 currently:
```cpp
add_tracer<Updated>("qv", m_grid, kg/kg, ps);
```

Change: No change needed if Phase 1.4 is done correctly! The `add_tracer` template handles layout automatically.

BUT: Need to update kernel calls inside `run_impl()` method.

**2.2 Similar for P3, RRTMGP, ZM, Homme, Surface Coupling**

Process interface files likely need no changes to registration if Phase 1.4 correct.

BUT: All kernel call sites need review.

#### Phase 3: Kernel updates (Tier 2 - largest effort)

**Strategy**: Use Kokkos subview to pass slot-0 to existing kernels, avoiding mass kernel rewrites.

Example in SHOC run_impl (file: `eamxx_shoc_process_interface.cpp`):

Current pattern (hypothetical):
```cpp
auto qv_view = get_field_out("qv").get_view<Real**>();
// Pass to kernel...
```

New pattern:
```cpp
auto qv_full = get_field_out("qv").get_view<Real***>();  // (tracer, col, lev)
auto qv_view = Kokkos::subview(qv_full, 0, Kokkos::ALL(), Kokkos::ALL());  // Slot 0 only
// Pass to kernel - no kernel changes needed!
```

**Files to examine** (NOT in current deliverables list):
- `components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp`
- `components/eamxx/src/physics/p3/impl/p3_main_impl.hpp`
- `components/eamxx/src/physics/rrtmgp/impl/*` (multiple files)
- `components/eamxx/src/physics/zm/impl/*`
- `components/eamxx/src/dynamics/homme/impl/*`

Each kernel invocation site needs audit:
1. Find all `get_field` calls for "qv"
2. Change view type from `Real**` to `Real***`
3. Add subview extraction for slot 0
4. Verify kernel gets correct 2D view

**Critical**: This is where BFB can break! Subview must not change memory layout.

#### Phase 4: Tests (Tier 1)

**4.1 Unit test: qv tracer access** (NEW file)

File: `components/eamxx/tests/water_tracers/test_qv_tracer_access.cpp`

Tests:
- Field has correct layout: `(SCREAM_NUM_TRACERS, ncol, nlev)`
- Slot 0 access equivalent to old scalar access pattern
- Tracer dimension correctly sized
- Subview correctness (memory layout unchanged)

**4.2 Component test: qv transport** (NEW file)

File: `components/eamxx/tests/water_tracers/test_qv_transport.cpp`

Tests:
- Initialize qv(0,*,*) with pattern A, qv(1,*,*) with pattern B
- Run one physics timestep
- Verify patterns remain independent (no cross-contamination)

**4.3 CMakeLists.txt** (NEW file)

File: `components/eamxx/tests/water_tracers/CMakeLists.txt`

Add both tests to build system with `SCREAM_NUM_TRACERS=1` and `=3` variants.

#### Phase 5: I/O and restart handling (Tier 2)

**Question for user**: Spec deliverables don't mention I/O or restart files. Are these handled automatically by FieldManager, or do we need explicit updates?

Files that might need review:
- `components/eamxx/src/share/io/*`
- History output writing
- Restart file reading/writing

### Implementation order recommendation

Given size and risk:

1. **Implement Phase 1** (field infrastructure) first, verify compilation
2. **Add Phase 4 tests** (write test scaffolding, expect failures)
3. **Implement Phase 2** (process interfaces, easiest if template correct)
4. **Implement Phase 3 incrementally** (one process at a time):
   - Start with surface_coupling (likely simplest)
   - Then SHOC
   - Then P3
   - Then others
5. **Run tests after each process** to catch BFB breaks early
6. **Phase 5 as needed** based on test results

## Testing

(Not yet started)

## Blockers

None currently identified.

## Notes

- Campaign spec 2/5 in wiso-group1 campaign
- Base branch: wiso/01-water-tracer-metadata-and-gate  
- Spec validated PASS by campaign precondition
