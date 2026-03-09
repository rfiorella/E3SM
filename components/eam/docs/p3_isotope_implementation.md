# P3 Microphysics Water Isotope Implementation — Development History and Status

## Overview

This document records the development history of the water isotope capability for
EAMv3 with P3 microphysics, implemented on branch `rfiorella/eam-v3/wiso`. The work
was done in March 2026, building on a 2022 port (commit `2f332ef044`) of CAM6 water
isotope infrastructure into EAMv3.

### Problem Statement

The 2022 port brought CAM6 water isotope infrastructure into EAMv3, but that code
coupled to MG2 microphysics via `wtrc_mg_inter()`. EAMv3 uses P3 microphysics, which
has no isotope support and uses a single unified ice category (no separate snow),
requiring a new approach to ice isotope fractionation.

### Goal

Create a working water isotope capability in EAMv3 with P3 microphysics, initially
supporting HDO and H2-18O.

---

## Modified Files

### 1. `src/physics/cam/water_tracers.F90` (+422 lines)

**New subroutine: `wtrc_p3_inter()`** (~300 lines)

The core P3 isotope interface, called after `p3_main` returns. It:

1. **Builds a process rate matrix** from P3's `tend_out` array, mapping P3 process
   rates to the existing isotope framework:

   | P3 tend_out Index | Process | Rate Matrix Mapping |
   |-------------------|---------|---------------------|
   | 17 (`qidep`) + 19 (`qinuc`) | Vapor deposition + nucleation | vapor -> ice (with Rayleigh fractionation) |
   | 23 (`qi2qv_sublim`) | Ice sublimation | ice -> vapor (proportional) |
   | 11 (`qr2qv_evap`) | Rain evaporation | strat_rain -> vapor (proportional) |
   | 3 (`qc2qr_autoconv`) + 2 (`qc2qr_accret`) | Autoconversion + accretion | liquid -> strat_rain (proportional) |
   | 28 (`qc2qi_hetero_freeze`) + 15 (`qccol`) | Heterogeneous freezing + collection | liquid -> ice (proportional, no fractionation) |
   | 18 (`qrcol`) + 29 (`qr2qi_immers_freeze`) | Rain collection + immersion freezing | strat_rain -> ice (proportional, no fractionation) |
   | 24 (`qi2qr_melt`) + 33 (`qc2qr_ice_shed`) | Ice melting + shedding | ice -> strat_rain (proportional) |
   | 50 (`qiberg`) | Bergeron (WBF) process | liquid -> ice via vapor (with fractionation) |

2. **Applies rates iteratively** (`wtrc_niter` iterations, default 10) with:
   - **Rayleigh distillation** for vapor -> ice deposition: `fd = qloc/qloc0`,
     `qloc_iso = qloc0_iso * fd^alpha`
   - **Bergeron handling**: Encoded as `(iwtliq, iwtliq, iwtliq)` in the rate matrix.
     `wtrc_apply_rates` interprets this triple as liquid -> vapor -> ice with
     fractionation applied during the vapor -> ice step.
   - **Proportional transfer** for all other processes (collection, freezing,
     melting, etc.): `delta_iso = (q_iso/q_bulk) * delta_bulk`

3. **Applies sedimentation proportionally** using P3 tend_out indices:
   - Index 36: cloud liquid sedimentation tendency (kg/kg/s)
   - Index 38: cloud ice sedimentation tendency (kg/kg/s)
   - Index 40: rain sedimentation tendency (kg/kg/s)
   - Formula: `delta_q_iso_sed = (q_iso / q_bulk) * delta_q_bulk_sed * dt`

4. **Updates `ptend%q`** for all isotope tracers.

**New subroutine: `wtrc_p3_set_precip()`** (~50 lines)

Sets isotope surface precipitation in physics buffer fields. Uses the isotope ratio
at the lowest model level to partition bulk precipitation into isotope precipitation:

```
precip_iso = (q_iso_bottom / q_bulk_bottom) * precip_bulk
```

This populates `wtrc_srfpcp_indices(water_type, water_set)` pbuf fields that are
read by `physpkg.F90` for surface coupling.

### 2. `src/physics/p3/eam/micro_p3.F90` (+5/-1 lines)

- Increased `p3_tend_out` dimension from 49 to 50 in two declaration sites
- Added `qiberg` (Bergeron/WBF process rate) as `p3_tend_out(k,50)` after the
  existing index 35 block in the `ice_deposition_sublimation()` subroutine

### 3. `src/physics/p3/eam/micro_p3_interface.F90` (+59/-6 lines)

- Added `use water_tracer_vars` and `use water_tracers` imports
- Updated `tend_out` dimension from 49 to 50
- Added isotope tracer registration in `micro_p3_register()` via `wtrc_register()`
- Added isotope tracer detection in `micro_p3_implements_cnst()` via `wtrc_implements_cnst()`
- Added isotope IC initialization in `micro_p3_init_cnst()` (set to 0)
- Added `wtrc_init()` call at end of `micro_p3_init()`
- Added isotope tracer indices to `lq` array before `physics_ptend_init`
- Added rime fraction computation from `qm/qi` (clamped to [0,1])
- Added `wtrc_p3_inter()` call after bulk P3 tendency computation
- Added `P3_qiberg` history field registration and output

---

## Unmodified Files (Verified Compatible)

These files were reviewed and confirmed to need no changes:

| File | Reason |
|------|--------|
| `src/physics/cam/microp_driver.F90` | P3 dispatch already routes isotope registration/init correctly |
| `src/physics/cam/clubb_intr.F90` | Existing `trace_water` guards sufficient; commented-out process rate code already inactive; CLUBB transports isotope tracers as passive scalars |
| `src/physics/cam/zm_conv_intr.F90` | 2022 port code compatible with current EAMv3 interfaces |
| `src/physics/cam/convect_shallow.F90` | 2022 port code compatible |
| `src/physics/cam/uwshcu.F90` | 2022 port code compatible |
| `src/control/camsrfexch.F90` | 2022 port code handles isotope surface exchange |
| `src/cpl/atm_import_export.F90` | 2022 port code handles isotope coupler fields |
| `src/physics/cam/cam_diagnostics.F90` | 2022 port code handles isotope diagnostics |
| `share/util/water_isotopes.F90` | Fractionation routines reused as-is |
| `share/util/water_types.F90` | Type indices reused as-is |
| `src/physics/cam/water_tracer_vars.F90` | Configuration variables reused as-is |

---

## Key Design Decisions

### 1. Single Ice Category (No Separate Snow)

P3 uses a single unified ice category with a rime mass fraction (`qm/qi`) that
continuously varies from 0 (fully vapor-deposited, cloud ice-like) to 1 (fully
rimed, graupel-like). This differs from MG2, which has separate cloud ice and snow
categories.

**Approach**: The rime fraction is implicitly tracked through isotope ice composition.
Deposition processes (vapor -> ice) apply Rayleigh fractionation, while
collection/riming processes (liquid -> ice, rain -> ice) inherit source composition
with no fractionation. The resulting isotope ice composition naturally reflects the
rime history without needing explicit rime-fraction weighting during rate application.

### 2. Bergeron Process Encoding

The Bergeron (Wegener-Bergeron-Findeisen) process transfers mass from liquid to ice
via the vapor phase. It is encoded in the process rate matrix as
`(iwtliq, iwtliq, iwtliq)` — destination, source, and rate-type all set to liquid.
The `wtrc_apply_rates()` routine recognizes this triple and applies two-step
fractionation: liquid -> vapor (inverse fractionation) then vapor -> ice (forward
fractionation).

### 3. Proportional Sedimentation

Sedimentation does not fractionate isotopes. The isotope sedimentation tendency is
computed proportionally:

```
delta_q_iso = (q_iso / q_bulk) * delta_q_bulk * dt
```

This uses P3's sedimentation tendencies from `tend_out` indices 36 (liquid), 38 (ice),
and 40 (rain), which are in units of kg/kg/s.

### 4. IC Initialization

Isotope tracers not found on the IC file are initialized to zero. This means the model
needs a spinup period for isotope fields to reach physically reasonable values. Proper
IC files with isotope fields at natural abundance ratios are needed for production runs.

---

## Implementation Status

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Build system & compilation | **Complete** |
| 1 | `wtrc_p3_inter()` core subroutine | **Complete** |
| 2 | Expose P3 Bergeron rate (`qiberg`) | **Complete** |
| 3 | Wire up in `micro_p3_interface.F90` | **Complete** |
| 4 | CLUBB verification (no changes needed) | **Complete** |
| 5 | Convection verification (no changes needed) | **Complete** |
| 6 | Surface coupling & diagnostics verification | **Complete** |
| 7 | Documentation (`adding_isotope_species.md`) | **Complete** |
| 8 | Testing & validation | **Not started** |

---

## Remaining Work: Phase 8 — Testing & Validation

### Compilation Test
- Build with `wisotope=.true.` and `microp_scheme='P3'`
- Verify clean compilation with no errors or warnings related to isotope code

### BFB (Bit-for-Bit) Verification
- Run with `wisotope=.false.` and compare against a baseline without any isotope
  code changes
- Must be identical — isotope code should have zero impact when disabled

### Smoke Test
- Run for 1 day (SCAM or short global) with `wisotope=.true.`
- Verify no crashes, NaNs, or negative mixing ratios in isotope tracers
- Check that isotope surface precipitation fields are populated

### Physical Validation (Aquaplanet or F-compset)
Expected ranges for a physically reasonable simulation:

| Diagnostic | Expected Range | Notes |
|------------|---------------|-------|
| delta-18O in precipitation | -30 to 0 per mil (mid-latitudes) | Temperature-dependent depletion |
| delta-D in precipitation | -240 to 0 per mil | ~8x delta-18O |
| Deuterium excess (d = dD - 8*d18O) | -10 to +30 per mil | Kinetic fractionation signal |
| delta-18O in vapor | -40 to -10 per mil | More depleted than precipitation |

### Known Risks and Potential Issues

1. **CLUBB passive transport**: Isotope tracers are transported by CLUBB as passive
   scalars without fractionation. This means turbulent mixing and boundary layer
   processes do not fractionate, which is physically correct for turbulent transport
   but means any CLUBB-specific condensation/evaporation is not captured for isotopes.
   This is acceptable for a first implementation.

2. **Convection code vintage**: The convection isotope code is from the 2022 CAM6
   port. It was verified compatible with current EAMv3 interfaces but has not been
   tested end-to-end with P3 in the loop. Convection is independent of the
   microphysics scheme, so this should work, but runtime verification is needed.

3. **Process rate completeness**: The P3 tend_out indices used in `wtrc_p3_inter()`
   cover the major microphysical processes but may not capture every minor pathway.
   If mass conservation checks show significant imbalances, additional tend_out
   indices may need to be mapped.

4. **Sedimentation timing**: The proportional sedimentation is applied after the
   process rate iteration. If P3's internal sedimentation occurs at a different
   point in the time step than assumed, this could introduce small biases.

---

## File Quick Reference

```
components/eam/
  src/physics/cam/
    water_tracers.F90          -- wtrc_p3_inter(), wtrc_p3_set_precip() [MODIFIED]
    water_tracer_vars.F90      -- shared config variables [unchanged]
    microp_driver.F90          -- dispatch to P3 [unchanged]
    clubb_intr.F90             -- CLUBB interface [unchanged]
    zm_conv_intr.F90           -- ZM deep convection [unchanged]
    convect_shallow.F90        -- shallow convection dispatch [unchanged]
    uwshcu.F90                 -- UW shallow convection [unchanged]
    cam_diagnostics.F90        -- diagnostics [unchanged]
  src/physics/p3/eam/
    micro_p3.F90               -- qiberg in tend_out(k,50) [MODIFIED]
    micro_p3_interface.F90     -- isotope registration, init, wtrc_p3_inter call [MODIFIED]
  src/control/
    camsrfexch.F90             -- surface exchange [unchanged]
  src/cpl/
    atm_import_export.F90      -- coupler interface [unchanged]
  share/util/
    water_isotopes.F90         -- species definitions, fractionation [unchanged]
    water_types.F90            -- water type indices [unchanged]
  docs/
    adding_isotope_species.md  -- how to add species beyond HDO+H218O [NEW]
    p3_isotope_implementation.md -- this file [NEW]
```
