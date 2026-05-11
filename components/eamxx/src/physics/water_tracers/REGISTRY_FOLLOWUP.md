# Water Tracer Registry - Follow-up Items

## Completed with this Implementation

✅ **Water tracer registry infrastructure** (spec: 2026-05-06-eamxx-water-tracer-isotopologue-registry.md)
- Compile-time registry mapping tracer index → isotopologue catalog index
- Per-tracer metadata: name, catalog_idx, is_tag
- CMake configuration via `add_water_tracer()` in config files
- Configure-time validation
- Query API: `tracer_isotopologue()`, `tracer_name()`, `tracer_is_tag()`, `find_tracer_by_name()`

## Unblocked - Ready for Implementation

The registry now provides all metadata needed for the following downstream work:

### 1. **Fractionation Physics Implementation** (HIGH PRIORITY)
**Status:** UNBLOCKED - registry provides tracer → isotopologue → fractionation factor pathway

The registry enables fractionation-aware physics via this pattern:
```cpp
// In a WaterTracerHook implementation
for (int i = 0; i < WTRC_MAX_CNST; ++i) {
  auto cat_idx = tracer_isotopologue(i);
  if (cat_idx == 0) continue;  // Skip H2O (no self-fractionation)
  if (tracer_is_tag(i)) continue;  // Skip passive tags
  
  // Get isotopologue name from catalog
  auto iso_name = WaterIsotopologueNames[cat_idx];
  
  // Compute fractionation factor
  auto alpha_liq_vap = AlphaEqLiquidVapor(iso_name, temperature);
  auto alpha_ice_vap = AlphaEqIceVapor(iso_name, temperature);
  
  // Apply to tendencies...
}
```

**Next spec should cover:**
- WaterTracerHook implementation for equilibrium fractionation
- Application points: condensation, evaporation, deposition, sublimation
- Integration with P3 microphysics phase changes
- Unit tests with known fractionation factors

### 2. **In-Substep Fractionation Hooks** (DEFERRED - temperature trajectory issue)
**Status:** Architecturally ready, physically blocked

The registry is compatible with in-substep hooks, but implementation is deferred pending resolution of the temperature-trajectory issue documented in spec `2026-05-06-eamxx-water-condensates-add-tracer-dim.md`.

**Blocker:** Mid-substep phase changes in P3/SHOC need instantaneous temperature, but the atm_process interface only provides column-mean T. Fractionation factors are temperature-dependent, so applying them at the wrong T introduces systematic error.

**When unblocked:** The registry's query API works unchanged - `tracer_isotopologue()` returns the catalog index needed to look up fractionation parameters.

### 3. **Tag Tracer Semantics** (NEW CAPABILITY)
**Status:** Infrastructure present, semantics undefined

The `is_tag` flag is wired up but has no physics associated with it yet. Tag tracers (bulk H2O that tracks source attribution) are distinguishable from true isotopologues.

**Next spec should define:**
- How tag tracers are initialized (boundary conditions, surface fluxes)
- Whether tags follow identical advection/diffusion to bulk or have custom treatment
- How tags evolve through convection (ZM) and microphysics (P3)
- Output and diagnostics

### 4. **Number Concentration Tracers** (DEFERRED)
**Status:** Not yet started

The registry established here applies equally to number concentration fields (nc, nr, ni, bm) once they are lifted to (COL, CMP, LEV). The same catalog-index mapping and query API will work.

**Deferred because:** Number concentrations add complexity (aerosol interactions, size-distribution coupling) beyond the scope of mass-based isotope tracers.

### 5. **IO Output of Per-Tracer Names** (DEFERRED)
**Status:** Registry provides names, but output plumbing not yet implemented

The registry exposes `tracer_name(i)` as a string_view, but NetCDF output dimension labels and per-slice metadata are not yet wired up.

**Next spec should cover:**
- NetCDF dimension naming for CMP axis (e.g., "water_tracer")
- Per-slice coordinate variable with tracer names
- Metadata attributes (isotopologue, is_tag)

### 6. **Restart / Checkpoint** (DEFERRED)
**Status:** Not yet addressed

Restarting a simulation with a different tracer configuration (different number or ordering) is not yet supported. The registry is rebuilt from the config file each time, but checkpoint files contain no tracer metadata.

**Design question:** Should checkpoints store tracer names and validate on restart? Or should users be responsible for maintaining config-file consistency?

## Out of Scope (Not Registry-Related)

The following items are orthogonal to the registry and not unlocked by it:

- **Surface fluxes** (isotope-aware ocean/land interface)
- **EAM legacy code** (components/eam/ Fortran)
- **Compset / coupler changes** (cime_config/)
- **Water phase enum** (WaterType in water_types.hpp - orthogonal to isotopologue)

## Implementation Notes for Downstream Specs

1. **Catalog index is the authoritative isotopologue identifier.** Do not hard-code tracer indices - always query `tracer_isotopologue(i)` and use the returned catalog index to look up physical parameters.

2. **The catalog ordering in `WaterIsotopologueNames` is fixed.** Appending new isotopologues to the end is allowed; reordering is a breaking change that invalidates all existing tracer config files.

3. **Multiple tracers can share a catalog index.** Tag tracers (passive H2O) and bulk H2O both map to catalog index 0. Use `tracer_is_tag(i)` to distinguish them.

4. **Registry queries are constexpr and device-callable.** Use them freely in Kokkos kernels - there is no runtime overhead.

5. **`find_tracer_by_name()` is host-only.** For device code, iterate and compare `tracer_name(i)` directly.

---

**Last updated:** 2026-05-07  
**Registry spec:** specs/2026-05-06-eamxx-water-tracer-isotopologue-registry.md
