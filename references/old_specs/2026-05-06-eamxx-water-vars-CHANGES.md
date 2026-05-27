# EAMxx qv Tracer Dimension - Implementation Changes Summary

## Overview

This document summarizes all code changes made to implement the qv tracer dimension feature (spec: 2026-05-06-eamxx-water-vars-add-tracer-dim.md).

## Design Decisions Implemented

1. **Compile-time count**: CMake preprocessor define `SCREAM_NUM_WATER_TRACERS` → `WTRC_MAX_CNST`
2. **Hook form**: Process-specific function pointer table (condensation, evaporation, deposition, sublimation, freezing, melting)
3. **qv always rank-3**: No `#ifdef` on field rank; (COL, CMP=1, LEV) when OFF, (COL, CMP=N, LEV) when ON

## Files Modified

### 1. CMake Configuration
**File:** `components/eamxx/src/physics/water_tracers/CMakeLists.txt`

- Added `option(SCREAM_TRACE_WATER "Enable water isotope tracer support" OFF)`
- Added `SCREAM_NUM_WATER_TRACERS` variable (default 1 when ON)
- Configure-time validation: rejects N≠1 when OFF
- Propagates `WTRC_MAX_CNST=${SCREAM_NUM_WATER_TRACERS}` as preprocessor define

### 2. Water Tracers Header
**File:** `components/eamxx/src/physics/water_tracers/water_tracers.hpp`

- Changed `WTRC_MAX_CNST` from `constexpr int = 1` to CMake-provided preprocessor define
- Added `#ifndef WTRC_MAX_CNST #error` guard
- Added `get_bulk_water_subview()` template function:
  - Extracts (COL, LEV) view from rank-3 (COL, CMP, LEV) field at CMP=0
  - Returns `Kokkos::subview(water_field_rank3, Kokkos::ALL(), 0, Kokkos::ALL())`
  - Compiled in all builds, zero cost with WTRC_MAX_CNST=1

### 3. Hook Interface Definition
**File (NEW):** `components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp`

Defines process-specific function pointer table:
- 6 hook types: `CondensationHookFn`, `EvaporationHookFn`, `DepositionHookFn`, `SublimationHookFn`, `FreezingHookFn`, `MeltingHookFn`
- Each has signature: `void(*)(view_3d qv, view_3d qc/qi, view_2d T, view_2d p, view_2d tendency, int ncol, int nlev)`
- Inline no-op defaults: `default_noop_condensation()`, etc. (eliminated by optimizer)
- Function pointers initialized to no-ops: `condensation_hook = default_noop_condensation`, etc.
- Declaration of `initialize_water_tracer_hooks()` (only when `SCREAM_TRACE_WATER` defined)

### 4. Hook Implementation
**File:** `components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp`

Added (under `#ifdef SCREAM_TRACE_WATER`):
- Passive-copy implementations for all 6 hooks
  - When N=1: immediate return (no-op)
  - When N>=2: loop over CMP indices 1..N-1, copy bulk tendency to each tracer
- `initialize_water_tracer_hooks()`: sets function pointers to passive-copy implementations

### 5. Core Tracer Allocation Change
**File:** `components/eamxx/src/share/atm_process/atmosphere_process.hpp`

- Added `#include "physics/water_tracers/water_tracers.hpp"`
- Modified `add_tracer()` template method (line ~388):
  - Changed from: `grid->get_3d_scalar_layout(FieldTag::LevelMidPoint)`
  - Changed to: `grid->get_3d_vector_layout(FieldTag::LevelMidPoint, scream::WaterTracers::WTRC_MAX_CNST, "water_tracer")`
  - Added comment explaining always-rank-3 design
  - Affects all `add_tracer<>("qv", ...)` calls throughout codebase

### 6. Consumer Updates (13 files)

All files updated to use `get_bulk_water_subview()` accessor:

**P3 Microphysics:**
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp` (2 sites: qv, qv_prev)

**SHOC Turbulence:**
- `components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp` (2 sites)

**ZM Convection:**
- `components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp` (1 site)

**Homme Dynamics:**
- `components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp` (1 site)

**IOP Forcing:**
- `components/eamxx/src/physics/iop_forcing/eamxx_iop_forcing_process_interface.cpp` (1 site)

**COSP Diagnostics:**
- `components/eamxx/src/physics/cosp/eamxx_cosp.cpp` (2 sites: device and host views)

**TMS:**
- `components/eamxx/src/physics/tms/eamxx_tms_process_interface.cpp` (1 site)

**MAM Aerosols:**
- `components/eamxx/src/physics/mam/eamxx_mam_generic_process_interface.cpp` (1 site)

**Surface Coupling:**
- `components/eamxx/src/control/atmosphere_surface_coupling_exporter.cpp` (1 site)

**Pattern for all consumers:**
```cpp
// Before:
const auto& qv = get_field_out("qv").get_view<Pack**>();

// After:
// qv is now rank-3 (COL, CMP, LEV); extract bulk water at CMP=0
const auto  qv_rank3 = get_field_out("qv").get_view<Pack***>();
const auto& qv       = scream::WaterTracers::get_bulk_water_subview(qv_rank3);
```

## Testing Strategy

### Build Configurations to Test
1. **Default (OFF)**: `SCREAM_TRACE_WATER=OFF` → WTRC_MAX_CNST=1
2. **ON N=1**: `SCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=1`
3. **ON N=2**: `SCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=2`

### Expected Behavior
- **Default vs ON N=1**: Bit-for-bit identical (BFB test)
  - Both use WTRC_MAX_CNST=1
  - Same field layout: (COL, CMP=1, LEV)
  - Same data path
  - Only difference: hook registration (but hooks are no-op at N=1)

- **ON N=2 passive-copy test**:
  - Tracer[0] = bulk water (physics-driven)
  - Tracer[1] initialized as copy of tracer[0]
  - After physics: tracer[1] == tracer[0] to machine epsilon
  - Validates end-to-end CMP dimension plumbing

## What's NOT in This Implementation

Per spec's out-of-scope section:
- Actual isotope fractionation physics (future spec)
- Hook call sites in P3/SHOC (not required for compilation; future spec)
- qc, qi, qr, qm rank lift (follow-up specs)
- Number concentration tracers (nc, ni, nr, bm)
- IO output of per-tracer names
- Restart/checkpoint
- Tag-tracer source attribution
- Surface fluxes for isotopes

## Follow-Up Specs Needed

1. **qc/qi/qr/qm rank lift**: Apply same pattern to condensed phases
2. **Tracer registry**: Build-time configuration of which isotopologues to track
3. **Fractionation physics**: Real isotope-aware hook implementations using `AlphaEq*` from water_isotopes.hpp
4. **In-substep hooks**: Insert hook calls at P3/SHOC phase-change sites
5. **Number concentration lift**: nc, ni, nr, bm to (COL, CMP, LEV)
6. **IO/restart**: NetCDF dimension labels, checkpoint support

## Notes for Reviewers

- The rank change affects ALL builds, not just isotope-enabled ones
- With WTRC_MAX_CNST=1, the CMP dimension is unit-stride and adds no runtime cost
- All consumers use the same accessor pattern (mechanical change)
- Hook interface is ready for future fractionation physics without touching call sites
- This is infrastructure only - no science changes yet
