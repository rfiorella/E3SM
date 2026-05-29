# EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension

## Campaign Overview

This campaign extends all EAMxx water reservoir fields to support multiple tracer
constituents along a new leading dimension. At completion:

- All water mass fields (`qv`, `qc`, `qi`, `qr`, `qs`, `qm`, `rain`, `snow`) have shape `(tracer, col, lev)` instead of `(col, lev)`
- Users register water tracers via `add_water_tracer(NAME <name> LONGNAME <desc> UNITS <units>)` in CMake
- **Slot 0 is canonical bulk water** and must remain bit-for-bit identical to pre-campaign scalar baseline
- Slots 1+ hold additional tracers registered via `add_water_tracer`
- A validation test demonstrates tracer transport correctness: surface evaporation scaled by 0.5 for tracer 1 produces precipitation 0.5x that of tracer 0

**Critical success criterion:** Existing EAMxx regression tests pass bit-for-bit (BFB) when `SCREAM_NUM_TRACERS=1` (bulk water only in slot 0).

**Campaign scope:** 5 PRs, structural changes only, no isotope physics, includes tracer ratio utility for validation.

**Target timeline:** 8-10 weeks from PR 1 approval to Group 1 boundary baseline validation.

## Binding Design Constraints

### Tracer dimension

- **Compile-time constant:** `SCREAM_NUM_TRACERS` (default 1)
  - Total slots allocated, including bulk water in slot 0
  - Minimum value: 1 (bulk only)
  - Maximum value: 8
  - Incremented automatically by each `add_water_tracer()` call
- **Dimension ordering:** Leading tracer dimension
  - New shape: `qv(tracer, col, lev)` instead of scalar `qv(col, lev)`
  - Rationale: Kokkos vectorization patterns favor column-innermost
  - All water reservoirs use this layout

### Slot-0 semantics (binding)

- **Slot 0 = canonical bulk water** (always present)
- Slot 0 contains the total physical water amount used by legacy EAMxx physics
- Existing physics code reads/writes slot 0 only
- Bulk water is **never** reconstructed by summing tracer slots
- **BFB requirement:** When `SCREAM_NUM_TRACERS=1`, slot 0 must be bit-for-bit identical to pre-campaign scalar baseline for all existing tests
- This is a **hard gate**. Group 1 cannot complete until BFB is proven.

### Tracer slots 1+

- Slots 1, 2, ... hold additional tracers registered via `add_water_tracer()`
- Each tracer applies to **all** water reservoirs (`qv`, `qc`, `qi`, `qr`, `qs`, `qm`, `rain`, `snow`)
- Tracers are independently transported
- No physics hooks in Group 1 — tracers are passively advected
- Future campaigns add fractionation physics

### Performance gate

- **Acceptable overhead:** `SCREAM_NUM_TRACERS=1` (bulk-only mode) vs scalar baseline < 2% runtime
- **Measurement:** Full atmosphere component timestep, 10-run median
- **Gate location:** PR 1 prototype phase
- **Fallback:** If >2%, PR 1 must document architectural alternative and get approval before PR 2

### Water reservoirs in scope

The following fields require tracer dimension extension:

| Field | Current shape | New shape | PR |
|-------|--------------|-----------|-----|
| `qv` (vapor) | `(col, lev)` | `(tracer, col, lev)` | 2 |
| `qc` (cloud liquid) | `(col, lev)` | `(tracer, col, lev)` | 3 |
| `qi` (cloud ice) | `(col, lev)` | `(tracer, col, lev)` | 3 |
| `qm` (cloud mixed-phase) | `(col, lev)` | `(tracer, col, lev)` | 3 |
| `qr` (rain) | `(col, lev)` | `(tracer, col, lev)` | 4 |
| `qs` (snow) | `(col, lev)` | `(tracer, col, lev)` | 4 |
| `rain` (precip rate) | `(col, lev)` | `(tracer, col, lev)` | 4 |
| `snow` (precip rate) | `(col, lev)` | `(tracer, col, lev)` | 4 |

**Out of scope:** Cloud fraction, turbulent kinetic energy, and other diagnostic fields that are not water mass.

### CMake integration: `add_water_tracer`

Signature:

```cmake
add_water_tracer(
  NAME <tracer_name>
  LONGNAME <descriptive_name>
)
```

Effect:
- Registers tracer metadata in `scream::water_tracers::registry`
- Increments `SCREAM_NUM_TRACERS` at configure time
- Ensures field manager allocates sufficient slots across all water reservoirs
- Tracer applies to **all** water reservoirs automatically (no RESERVOIRS argument needed)
- Units are implicitly kg/kg (water mass mixing ratio) for all water tracers
- Example: Calling `add_water_tracer` once sets `SCREAM_NUM_TRACERS=2` (slot 0=bulk, slot 1=tracer)

Does **not** automatically add fractionation physics (that's future campaigns).

## Pre-Campaign Bootstrap

Before PR 1 opens, the following must exist:

### 1. `wiso-dev` long-lived branch

Created off `master` at campaign start:

```bash
git checkout master
git pull origin master
git checkout -b wiso-dev
git push -u origin wiso-dev
```

All Group 1 PRs merge into `wiso-dev`. After Group 1 completes and baselines pass, `wiso-dev` is proposed as a PR into `master`.

### 2. Pre-campaign Tier-2 baseline

Generate baseline at the SHA where `wiso-dev` branches from `master`.

**Procedure:**

```bash
cd components/eamxx

# Create baseline case
./scripts/test-all-eamxx -m <MACHINE> --baseline-dir /data/baselines/pre-group1 -g

# Record SHA, date, and config in /data/baselines/pre-group1/BASELINE.txt
```

Required baseline configurations:
- Machine: `docker-local` (or your HPC machine)
- Test types: `dbg`, `sp`, `fpe`

**Deliverable:** `/data/baselines/pre-group1/` with test outputs.

### 3. Directory structure

```bash
mkdir -p components/eamxx/src/physics/water_tracers
mkdir -p components/eamxx/tests/water_tracers
mkdir -p components/eamxx/tests/water_tracers/data
mkdir -p components/eamxx/docs/wiso
```

Commit this structure to `wiso-dev` before PR 1.

## PR Structure

### PR 1: Water-tracer metadata, types, and BFB feasibility gate

**Branch:** `wiso-01-tracer-metadata-and-gate`  
**Base:** `wiso-dev`  
**Tier:** 0 (no run required for merge, but prototype must run)  
**Estimated size:** 12 files, 600 lines

#### Dependencies

- Pre-campaign baseline exists at `/data/baselines/pre-group1/`
- `wiso-dev` branch exists and is up-to-date with `master`

#### Deliverables

1. **Types and enums** (`src/physics/water_tracers/water_tracer_metadata.hpp`):
   
   ```cpp
   namespace scream {
   namespace water_tracers {
   
   enum class WaterTracerKind {
     bulk,                  // Slot 0, canonical total water
     evaporation,          // Tracer with modified surface evaporation flux (Group 1 testing)
     // Future campaigns will add:
     // isotope, tagged_partition, tagged_diagnostic, auxiliary
     // See wiso_campaign_plan_revision1.md for full taxonomy
   };
   
   struct WaterTracerMetadata {
     std::string name;
     std::string longname;
     WaterTracerKind kind;
     // Units implicitly kg/kg for all water tracers in Group 1
     
     // Fields below reserved for future campaigns
     // See wiso_campaign_plan_revision1.md for full metadata specification
     bool conserved_independently;
     std::string ratio_standard;  // For isotope tracers
     std::string partition_group_id;  // For tagged_partition tracers
   };
   
   class WaterTracerRegistry {
     // Singleton, compile-time registration
     static WaterTracerRegistry& instance();
     void register_tracer(const WaterTracerMetadata& meta);
     const WaterTracerMetadata& get(const std::string& name) const;
     int count() const;
   };
   
   } // namespace water_tracers
   } // namespace scream
   ```

2. **CMake function** (`src/physics/water_tracers/CMakeLists.txt`):
   
   ```cmake
   function(add_water_tracer)
     set(options "")
     set(oneValueArgs NAME LONGNAME KIND RATIO_STANDARD PARTITION_GROUP_ID)
     set(multiValueArgs "")
     cmake_parse_arguments(TRACER "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
     
     # Validate required arguments
     if(NOT TRACER_NAME OR NOT TRACER_LONGNAME)
       message(FATAL_ERROR "add_water_tracer requires NAME and LONGNAME")
     endif()
     
     # Default KIND to evaporation for Group 1 test tracers
     if(NOT TRACER_KIND)
       set(TRACER_KIND "evaporation")
     endif()
     
     # Increment SCREAM_NUM_TRACERS
     math(EXPR SCREAM_NUM_TRACERS "${SCREAM_NUM_TRACERS} + 1")
     set(SCREAM_NUM_TRACERS ${SCREAM_NUM_TRACERS} PARENT_SCOPE)
     
     # Register metadata (units implicitly kg/kg)
     # (implementation generates water_tracer_config.hpp)
   endfunction()
   ```
   
   Sets `SCREAM_NUM_TRACERS` (default 1, meaning bulk-only).
   Generates `water_tracer_config.hpp` with tracer count and metadata.
   Units are implicitly kg/kg for all water tracers.

3. **Config header template** (`src/physics/water_tracers/water_tracer_config.hpp.in`):
   
   ```cpp
   #ifndef WATER_TRACER_CONFIG_HPP
   #define WATER_TRACER_CONFIG_HPP
   
   namespace scream {
   namespace water_tracers {
   
   constexpr int NUM_TRACERS = @SCREAM_NUM_TRACERS@;
   
   // Generated tracer metadata array
   // ...
   
   } // namespace water_tracers
   } // namespace scream
   
   #endif
   ```

4. **Design document** (`docs/wiso/tracer_data_model.md`):
   
   Must define:
   - Slot-0 semantics (canonical bulk, never reconstructed from sum)
   - `SCREAM_NUM_TRACERS` meaning: total slots including bulk in slot 0
   - Array indexing: slot 0 = bulk, slot 1+ = additional tracers
   - BFB feasibility result (see gate below)
   - Field layout decision (leading tracer dim)
   - Performance impact (<2% or fallback justification)
   - Tracer application scope (all reservoirs)
   - H216O implicit/explicit decision (deferred to isotope campaign)

5. **Prototype** (`src/physics/water_tracers/prototype/qv_extension_test.cpp`, not merged into production):
   
   Minimal Kokkos kernel extending `qv(col,lev)` → `qv(tracer,col,lev)`:
   
   ```cpp
   // Prototype test structure
   // 1. Allocate qv_scalar(ncol, nlev) and qv_tracer(1, ncol, nlev)
   // 2. Initialize with same values
   // 3. Run equivalent operations on both
   // 4. Verify bit-for-bit identical results
   // 5. Benchmark runtime for both versions
   ```
   
   Measures:
   - Runtime (target <2% overhead)
   - Memory bandwidth
   - Vectorization efficiency

#### BFB Feasibility Gate (blocking)

PR 1 includes a **mandatory checkpoint** before PRs 2-5 can start.

**Task:** Run the `qv_extension_test` prototype on target hardware and measure:

1. **Performance:** Compare `SCREAM_NUM_TRACERS=1` multidimensional access vs scalar
   - Pass: <2% overhead
   - Fail: Document why, propose alternative (e.g., template specialization for scalar path)

2. **BFB preservation:** Compare prototype output with known-good scalar qv kernel
   - Pass: Bit-identical
   - Fail: Document source of difference (packing? ordering? FP contraction?)

**Decision point:**

- If **both pass**: Update `tracer_data_model.md` with results, PRs 2-5 proceed with current design, BFB is a hard gate
- If **performance fails**: Update `tracer_data_model.md` with fallback (e.g., `if constexpr (SCREAM_NUM_TRACERS==1)` scalar path), get human approval, then proceed
- If **BFB fails**: Pause campaign, diagnose root cause, update design doc, get human approval for relaxed BFB tolerance or architectural change

**Approval:** A human must sign off on the design document and gate results before PR 1 merges.

#### Success Criteria

- [ ] `water_tracer_metadata.hpp` compiles without warnings
- [ ] `add_water_tracer(NAME test LONGNAME "Test Tracer")` generates valid config header
- [ ] `qv_extension_test` prototype compiles and runs
- [ ] Prototype shows <2% overhead OR `tracer_data_model.md` documents approved fallback
- [ ] Prototype output is BFB vs scalar OR `tracer_data_model.md` documents approved relaxation
- [ ] `docs/wiso/tracer_data_model.md` reviewed and approved by human
- [ ] No existing tests are modified (this PR is headers/docs/prototype only)

#### Out of Scope

- Actual field registration changes (that's PR 2)
- Isotope-specific logic (future campaign)
- Fractionation functions (future campaign)
- Multi-tracer validation test (PR 5)

---

### PR 2: Extend `qv` to tracer dimension

**Branch:** `wiso-02-extend-qv`  
**Base:** `wiso-01-tracer-metadata-and-gate`  
**Tier:** 2 (full Tier-2 validation required)  
**Estimated size:** 35 files, 900 lines

#### Dependencies

- PR 1 merged with approved `tracer_data_model.md`
- BFB feasibility gate passed

#### Deliverables

1. **Field layout extension** (`src/share/field/field_layout.hpp`):
   
   Add `TRACER` dimension tag:
   
   ```cpp
   enum FieldTag {
     // ... existing tags
     TRACER,  // New: tracer dimension
   };
   
   // New layout type
   FieldLayout tracer3d_mid_layout(int num_tracers, int ncols, int nlevs) {
     return FieldLayout({TRACER, COL, LEV}, {num_tracers, ncols, nlevs});
   }
   ```

2. **Field registration** in all processes that touch `qv`:
   
   Modify `add_field<Required>("qv", ...)` calls to use tracer-aware layout:
   
   ```cpp
   // Old
   add_field<Required>("qv", scalar3d_mid, kg/kg, grid);
   
   // New
   const int num_tracers = water_tracers::NUM_TRACERS;
   auto layout = grid.get_3d_vector_layout(true, num_tracers, "tracer");
   add_field<Required>("qv", layout, kg/kg, grid);
   ```
   
   Example locations:
   - `src/physics/shoc/eamxx_shoc_process_interface.cpp`
   - `src/physics/p3/eamxx_p3_process_interface.cpp`
   - `src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp`
   - `src/physics/zm/eamxx_zm_process_interface.cpp`
   - `src/dynamics/homme/eamxx_homme_process_interface.cpp`
   - `src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp`

3. **Kernel updates** to access `qv` with tracer slice:
   
   Pattern: `qv(icol, ilev)` → `qv(0, icol, ilev)` for all slot-0 bulk access
   
   Use Kokkos subview where entire tracer dimension is needed:
   
   ```cpp
   auto qv_bulk = Kokkos::subview(qv, 0, Kokkos::ALL, Kokkos::ALL);
   // Pass qv_bulk to existing kernels that expect (col, lev)
   ```
   
   Affected files (~25 files across `src/physics/`):
   - All parameterization interface files
   - All kernel implementation files that read/write `qv`
   - Diagnostic computation files

4. **Field manager integration** (`src/share/field/field_manager.cpp`, `src/share/field/field_manager.hpp`):
   
   - Update field allocation to respect `SCREAM_NUM_TRACERS`
   - Update restart I/O to handle tracer dimension (write all slots)
   - Update history I/O to handle tracer dimension (write all slots)
   - Add helper functions for tracer field access

5. **Unit tests** (`tests/water_tracers/test_qv_tracer_access.cpp`):
   
   ```cpp
   TEST_CASE("qv_tracer_dimension") {
     // Verify slot-0 access is equivalent to old scalar access
     // Verify tracer dimension is correctly sized (SCREAM_NUM_TRACERS)
     // When SCREAM_NUM_TRACERS=1, verify only 1 slot exists
     // Test field registration with tracer layout
     // Test Kokkos subview access patterns
   }
   ```

6. **Component test** (`tests/water_tracers/test_qv_transport.cpp`):
   
   Single-process test that initializes qv slots with different values and verifies they remain independent through a physics timestep.

#### Success Criteria

- [ ] All files compile without warnings with `-DSCREAM_NUM_TRACERS=1`
- [ ] All files compile without warnings with `-DSCREAM_NUM_TRACERS=3`
- [ ] Standalone `qv` unit test passes
- [ ] **BFB requirement:** All existing EAMxx tests pass BFB vs pre-campaign baseline with `SCREAM_NUM_TRACERS=1`
  - Verified via `cprnc` on all 3D/2D output fields
  - Zero tolerance: rtol=0, atol=0
  - Command: `./scripts/test-all-eamxx -m <MACHINE> --baseline-dir /data/baselines/pre-group1 -c`
- [ ] Performance regression: <2% runtime increase on full atmosphere timestep
- [ ] No segfaults, no MPI hangs, no NaN/Inf in output
- [ ] Component test `test_qv_transport` passes

#### Out of Scope

- Other water species (qc, qi, qr, qs, qm, rain, snow) — those are PRs 3-4
- Multi-tracer validation test (`SCREAM_NUM_TRACERS > 1`) — that's PR 5
- Isotope-specific kernels

---

### PR 3: Extend `qc`, `qi`, `qm` to tracer dimension

**Branch:** `wiso-03-extend-cloud`  
**Base:** `wiso-02-extend-qv`  
**Tier:** 2  
**Estimated size:** 40 files, 1000 lines

#### Dependencies

- PR 2 merged
- PR 2 BFB validation passed

#### Deliverables

Same pattern as PR 2, applied to cloud water species:

1. **Field registration** for `qc`, `qi`, `qm`:
   
   Update in:
   - `src/physics/p3/` (primary microphysics)
   - `src/physics/shoc/` (cloud fraction, turbulence)
   - `src/physics/cld_fraction/`
   - `src/physics/rrtmgp/` (radiation needs cloud water)
   - `src/physics/mam/` (aerosol activation, if present)

2. **Kernel updates**:
   
   All `qc(icol, ilev)` → `qc(0, icol, ilev)`
   All `qi(icol, ilev)` → `qi(0, icol, ilev)`
   All `qm(icol, ilev)` → `qm(0, icol, ilev)` if qm exists
   
   Use subviews for passing to existing kernels:
   
   ```cpp
   auto qc_bulk = Kokkos::subview(qc, 0, Kokkos::ALL, Kokkos::ALL);
   auto qi_bulk = Kokkos::subview(qi, 0, Kokkos::ALL, Kokkos::ALL);
   ```

3. **Derived fields** that compute from qc/qi:
   
   Fields like `inv_qc_relvar`, `eff_radius_qc`, `eff_radius_qi` may remain scalar (no tracer dim) if they're diagnostic-only and computed solely from slot-0 bulk values.

4. **Unit tests** (`tests/water_tracers/test_cloud_tracer_access.cpp`):
   
   ```cpp
   TEST_CASE("cloud_tracer_dimension") {
     // Verify qc, qi, qm have correct tracer dimension
     // Verify slot-0 access patterns
     // Verify subview access
   }
   ```

5. **Component test** (`tests/water_tracers/test_cloud_transport.cpp`):
   
   Similar to qv transport test, verifies cloud water tracers remain independent.

#### Success Criteria

- [ ] Compiles without warnings with `SCREAM_NUM_TRACERS=1` and `SCREAM_NUM_TRACERS=3`
- [ ] Unit tests pass
- [ ] **BFB requirement:** All existing tests pass BFB vs pre-campaign baseline with `SCREAM_NUM_TRACERS=1`
- [ ] Performance: <2% cumulative overhead (PR 2 + PR 3) vs pre-campaign baseline
- [ ] P3 microphysics single-process test passes BFB
- [ ] SHOC single-process test passes BFB
- [ ] Component test passes
- [ ] No segfaults, no MPI hangs, no NaN/Inf in output

#### Out of Scope

- Precipitation species (qr, qs, rain, snow) — that's PR 4
- Phase-change fractionation logic
- Multi-tracer validation test — that's PR 5

---

### PR 4: Extend `qr`, `qs`, `rain`, `snow` to tracer dimension

**Branch:** `wiso-04-extend-precip`  
**Base:** `wiso-03-extend-cloud`  
**Tier:** 2  
**Estimated size:** 35 files, 800 lines

#### Dependencies

- PR 3 merged
- PR 3 BFB validation passed

#### Deliverables

Same pattern as PRs 2-3, applied to precipitation:

1. **Field registration** for `qr`, `qs`, `rain`, `snow`:
   
   Update in:
   - `src/physics/p3/` (handles rain/snow microphysics)
   - Any sedimentation/transport code

2. **Kernel updates**:
   
   - `qr(icol, ilev)` → `qr(0, icol, ilev)`
   - `qs(icol, ilev)` → `qs(0, icol, ilev)`
   - `rain(icol, ilev)` → `rain(0, icol, ilev)`
   - `snow(icol, ilev)` → `snow(0, icol, ilev)`
   
   Use subviews:
   
   ```cpp
   auto qr_bulk = Kokkos::subview(qr, 0, Kokkos::ALL, Kokkos::ALL);
   auto qs_bulk = Kokkos::subview(qs, 0, Kokkos::ALL, Kokkos::ALL);
   auto rain_bulk = Kokkos::subview(rain, 0, Kokkos::ALL, Kokkos::ALL);
   auto snow_bulk = Kokkos::subview(snow, 0, Kokkos::ALL, Kokkos::ALL);
   ```

3. **P3 process-specific updates**:
   
   - Sedimentation must handle tracer dimension (passive advection for slots 1+)
   - Autoconversion, collection, etc. remain bulk-only (slot 0)
   - For slots 1+, apply same process rates as slot 0 (passive tracer behavior)

4. **Unit tests** (`tests/water_tracers/test_precip_tracer_access.cpp`):
   
   ```cpp
   TEST_CASE("precip_tracer_dimension") {
     // Verify qr, qs, rain, snow have correct tracer dimension
     // Verify slot-0 access patterns
   }
   ```

5. **Component test** (`tests/water_tracers/test_precip_transport.cpp`):
   
   Similar to previous transport tests.

#### Success Criteria

- [ ] Compiles without warnings with `SCREAM_NUM_TRACERS=1` and `SCREAM_NUM_TRACERS=3`
- [ ] Unit tests pass
- [ ] **BFB requirement:** All existing tests pass BFB vs pre-campaign baseline with `SCREAM_NUM_TRACERS=1`
- [ ] Performance: <2% cumulative overhead (PRs 2+3+4) vs pre-campaign baseline
- [ ] P3 microphysics with precipitation passes BFB
- [ ] Full atmosphere multi-process test passes BFB
- [ ] Component test passes
- [ ] No segfaults, no MPI hangs, no NaN/Inf in output

#### Out of Scope

- Isotope-aware sedimentation (future campaign)
- Multi-tracer validation test — that's PR 5

---

### PR 5: Tracer ratio utility and validation test

**Branch:** `wiso-05-tracer-validation`  
**Base:** `wiso-04-extend-precip`  
**Tier:** 1 (doubly-periodic test + component validation)  
**Estimated size:** 15 files, 500 lines

#### Dependencies

- PR 4 merged
- PR 4 BFB validation passed
- All water reservoirs now have tracer dimension

#### Deliverables

1. **Tracer ratio utility** (`src/physics/water_tracers/water_tracer_ratio.hpp`):
   
   ```cpp
   namespace scream {
   namespace water_tracers {
   
   // Compute ratio of tracer N to tracer 0 (bulk)
   // Returns tracer_N / tracer_0, with safeguards for division by zero
   template<typename ScalarT>
   KOKKOS_INLINE_FUNCTION
   ScalarT compute_tracer_ratio(
     const ScalarT& tracer_amount,  // Slot N amount
     const ScalarT& bulk_amount,    // Slot 0 amount
     const ScalarT& min_bulk = 1e-20  // Minimum bulk for valid ratio
   ) {
     if (bulk_amount < min_bulk) {
       return 0.0;  // Or NaN, depending on use case
     }
     return tracer_amount / bulk_amount;
   }
   
   // Kokkos-parallel version for fields
   template<typename ViewT>
   void compute_tracer_ratio_field(
     const ViewT& tracer_field,  // Shape: (tracer, col, lev)
     int tracer_idx,             // Which tracer to compute ratio for
     ViewT& ratio_field          // Output: (col, lev)
   );
   
   } // namespace water_tracers
   } // namespace scream
   ```

2. **Test case configuration** (`tests/water_tracers/test_tracer_scaling.cpp`):
   
   Doubly-periodic test case that:
   - Registers one additional water tracer via `add_water_tracer(NAME test_tracer LONGNAME "Test Evaporation Tracer")`
   - Tracer has KIND=evaporation, indicating its evaporation flux is modified
   - Sets `SCREAM_NUM_TRACERS=2` (slot 0=bulk, slot 1=test_tracer)
   - Initializes test_tracer to exactly 0.5 * bulk water at t=0 for all reservoirs
   - Scales surface evaporation flux for slot 1 by exactly 0.5
   - Runs for sufficient timesteps for precipitation to develop and reach quasi-steady state
   - Verifies tracer ratios remain at 0.5 to within numerical precision (rtol=1e-12)

3. **Surface flux modification and initialization** (test harness code):
   
   ```cpp
   // Initial conditions: set tracer 1 = 0.5 * tracer 0 for all reservoirs
   void initialize_test_tracers(FieldManager& field_mgr) {
     const auto qv = field_mgr.get("qv");
     const auto qc = field_mgr.get("qc");
     const auto qi = field_mgr.get("qi");
     const auto qr = field_mgr.get("qr");
     const auto qs = field_mgr.get("qs");
     
     const int ncols = qv.extent(1);
     const int nlevs = qv.extent(2);
     
     Kokkos::parallel_for(ncols, KOKKOS_LAMBDA(int icol) {
       for (int ilev = 0; ilev < nlevs; ++ilev) {
         // Initialize tracer 1 to exactly 0.5 * bulk
         qv(1, icol, ilev) = 0.5 * qv(0, icol, ilev);
         qc(1, icol, ilev) = 0.5 * qc(0, icol, ilev);
         qi(1, icol, ilev) = 0.5 * qi(0, icol, ilev);
         qr(1, icol, ilev) = 0.5 * qr(0, icol, ilev);
         qs(1, icol, ilev) = 0.5 * qs(0, icol, ilev);
       }
     });
   }
   
   // Surface flux: scale tracer 1 by exactly 0.5
   void apply_scaled_surface_flux(FieldManager& field_mgr) {
     auto qv_surf_flux = field_mgr.get("qv_surf_flux");
     const int ncols = qv_surf_flux.extent(1);
     
     // Slot 0: normal bulk flux (computed by physics)
     // Slot 1: scaled by exactly 0.5
     Kokkos::parallel_for(ncols, KOKKOS_LAMBDA(int icol) {
       qv_surf_flux(1, icol) = 0.5 * qv_surf_flux(0, icol);
     });
   }
   ```

4. **Validation script** (`tests/water_tracers/validate_tracer_scaling.py`):
   
   ```python
   #!/usr/bin/env python3
   import netCDF4 as nc
   import numpy as np
   
   def validate_tracer_scaling(output_file, rtol=1e-12, atol=1e-15):
       """
       Validate that tracer 1 / tracer 0 ratios are preserved to numerical precision.
       
       For passive tracers with no fractionation physics, the ratio should be
       preserved exactly except for floating-point roundoff error.
       
       Args:
           output_file: Path to history output NetCDF
           rtol: Relative tolerance (default 1e-12, suitable for double precision)
           atol: Absolute tolerance for near-zero values (default 1e-15)
       
       Returns:
           True if validation passes, raises AssertionError otherwise
       """
       ds = nc.Dataset(output_file)
       
       reservoirs = ['qv', 'qc', 'qi', 'qr', 'qs', 'rain', 'snow']
       expected_ratio = 0.5
       
       print(f"\nValidating tracer ratio preservation (expected ratio = {expected_ratio})")
       print(f"Tolerance: rtol={rtol}, atol={atol}\n")
       
       all_passed = True
       
       for var in reservoirs:
           if var not in ds.variables:
               continue
           
           # Shape: (time, tracer, col, lev) or (time, tracer, col)
           data = ds.variables[var][:]
           
           # Get final timestep
           bulk = data[-1, 0, ...]  # Slot 0
           tracer = data[-1, 1, ...]  # Slot 1
           
           # Compute ratio where bulk > threshold
           # Use threshold based on atol to avoid division by near-zero
           threshold = max(1e-20, atol * 100)
           mask = bulk > threshold
           
           if not np.any(mask):
               print(f"{var}: No significant bulk water (all values < {threshold})")
               continue
           
           ratio = tracer[mask] / bulk[mask]
           
           # Check that ratio is within tolerance of expected value
           # Use numpy.allclose logic: |ratio - expected| <= atol + rtol * expected
           ratio_error = np.abs(ratio - expected_ratio)
           tolerance = atol + rtol * expected_ratio
           
           passed = np.all(ratio_error <= tolerance)
           
           # Compute statistics
           mean_ratio = np.mean(ratio)
           max_error = np.max(ratio_error)
           max_rel_error = max_error / expected_ratio if expected_ratio > 0 else max_error
           
           status = "PASS" if passed else "FAIL"
           
           print(f"{var}:")
           print(f"  Mean ratio: {mean_ratio:.15e}")
           print(f"  Max absolute error: {max_error:.3e}")
           print(f"  Max relative error: {max_rel_error:.3e}")
           print(f"  Points checked: {np.sum(mask)}")
           print(f"  Status: {status}")
           
           if not passed:
               all_passed = False
               # Find worst violation
               worst_idx = np.unravel_index(np.argmax(ratio_error), ratio_error.shape)
               print(f"  Worst violation at index {worst_idx}:")
               print(f"    bulk = {bulk[mask][np.argmax(ratio_error)]:.15e}")
               print(f"    tracer = {tracer[mask][np.argmax(ratio_error)]:.15e}")
               print(f"    ratio = {ratio[np.argmax(ratio_error)]:.15e}")
           
           print()
       
       ds.close()
       
       if not all_passed:
           raise AssertionError("Tracer ratio validation failed for one or more reservoirs")
       
       print("All reservoirs passed tracer ratio validation\n")
       return True
   
   if __name__ == '__main__':
       import sys
       output_file = sys.argv[1]
       
       # Allow override of tolerances from command line
       rtol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-12
       atol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-15
       
       validate_tracer_scaling(output_file, rtol, atol)
   ```

5. **Test integration** (add to EAMxx test suite):
   
   ```cmake
   # In tests/water_tracers/CMakeLists.txt
   
   add_test(
     NAME water_tracer_scaling_dp
     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_tracer_scaling_test.sh
   )
   ```
   
   Shell script `run_tracer_scaling_test.sh`:
   
   ```bash
   #!/bin/bash
   # Run DP test with tracer scaling
   # Verify output with validation script
   
   set -e
   
   # Setup test case
   TEST_DIR=${CMAKE_BINARY_DIR}/tests/water_tracers/scaling_test
   mkdir -p ${TEST_DIR}
   cd ${TEST_DIR}
   
   # Create test with add_water_tracer
   ${EAMXX_ROOT}/scripts/create-test-case.sh \
     --case-name tracer_scaling \
     --compset F2010-SCREAMv1-DP \
     --res ne4pg2_ne4pg2 \
     --tracer-config "add_water_tracer(NAME test_tracer LONGNAME 'Test Evaporation Tracer')"
   
   # Build
   make -j4
   
   # Run
   ./run_test.sh
   
   # Validate
   python3 ${CMAKE_CURRENT_SOURCE_DIR}/validate_tracer_scaling.py output.nc
   ```

6. **Documentation update** (`docs/wiso/tracer_data_model.md`):
   
   Add section on tracer ratio computation and validation methodology.

#### Success Criteria

- [ ] `compute_tracer_ratio` utility compiles without warnings
- [ ] Unit test for `compute_tracer_ratio` passes
- [ ] Test case with `SCREAM_NUM_TRACERS=2` and `add_water_tracer` compiles
- [ ] Test case runs to completion without errors
- [ ] **Tracer ratio validation:** `validate_tracer_scaling.py` passes with rtol=1e-12, atol=1e-15
  - All water reservoirs (qv, qc, qi, qr, qs, rain, snow) maintain ratio of 0.5 to within numerical precision
  - Maximum relative error < 1e-12 for double precision arithmetic (< 1e-6 for single precision builds)
  - This validates that passive tracer transport preserves ratios exactly (no spurious numerical diffusion or corruption)
  - **This is the critical test that Group 1 is correct** - if ratios drift, there's a bug in tracer handling
- [ ] No segfaults, no MPI hangs, no NaN/Inf in output
- [ ] Performance: <2% overhead for `SCREAM_NUM_TRACERS=2` vs `SCREAM_NUM_TRACERS=1`
- [ ] Documentation updated

#### Out of Scope

- Isotope fractionation (future campaign)
- Multiple tracers beyond test_tracer
- Tagged tracers
- Production-ready surface flux infrastructure for tracers (basic test harness only)

---

## Group 1 Boundary: Baseline Regeneration and Final Validation

After PR 5 merges into `wiso-dev`, perform final Group 1 validation:

### 1. Regenerate baseline with `SCREAM_NUM_TRACERS=1`

```bash
cd components/eamxx
git checkout wiso-dev
git pull origin wiso-dev

# Generate new baseline with tracer dimension present, SCREAM_NUM_TRACERS=1
./scripts/test-all-eamxx -m <MACHINE> --baseline-dir /data/baselines/group1-boundary -g

# Record SHA, date, SCREAM_NUM_TRACERS=1 in /data/baselines/group1-boundary/BASELINE.txt
```

### 2. Compare against pre-campaign baseline

```bash
# Use cprnc to compare each test output
./scripts/compare-baselines.sh \
  /data/baselines/pre-group1 \
  /data/baselines/group1-boundary \
  --rtol 0 --atol 0
```

Expected result: **BFB identical** for all fields in all tests when `SCREAM_NUM_TRACERS=1`.

If not BFB:
1. Identify which field(s) differ
2. Bisect PRs 2-5 to find introduction point
3. Fix and regenerate baseline
4. Repeat until BFB

**Group 1 is not complete until this passes.**

### 3. Run tracer scaling validation test

```bash
cd tests/water_tracers
./run_tracer_scaling_test.sh
```

Expected result: Validation script passes with all ratios within tolerance.

### 4. Performance validation

Run full atmosphere timestep benchmark:

```bash
# Baseline (SCREAM_NUM_TRACERS=1)
./scripts/benchmark-eamxx.sh --tracers 1 > bench_1tracer.txt

# With additional tracer (SCREAM_NUM_TRACERS=2)
./scripts/benchmark-eamxx.sh --tracers 2 > bench_2tracers.txt

# Compare runtimes
python3 scripts/compare_benchmarks.py bench_1tracer.txt bench_2tracers.txt
```

Expected result: <5% overhead for `SCREAM_NUM_TRACERS=2` vs `SCREAM_NUM_TRACERS=1`.

---

## Branch and Baseline Model

### Branch structure

```
master
  └─ wiso-dev (long-lived)
       ├─ wiso-01-tracer-metadata-and-gate
       ├─ wiso-02-extend-qv (bases on wiso-01)
       ├─ wiso-03-extend-cloud (bases on wiso-02)
       ├─ wiso-04-extend-precip (bases on wiso-03)
       └─ wiso-05-tracer-validation (bases on wiso-04)
```

### Merge strategy

1. PR 1 opens against `wiso-dev`
2. After PR 1 merges, PR 2 opens against `wiso-dev` (rebased from `wiso-01` branch)
3. After PR 2 merges, PR 3 opens against `wiso-dev` (rebased from `wiso-02` branch)
4. After PR 3 merges, PR 4 opens against `wiso-dev` (rebased from `wiso-03` branch)
5. After PR 4 merges, PR 5 opens against `wiso-dev` (rebased from `wiso-04` branch)
6. After PR 5 merges and Group 1 boundary validation passes, open final PR: `wiso-dev → master`

### Baseline naming

- `pre-group1/` — generated from `master` before campaign starts
- `group1-boundary/` — generated from `wiso-dev` after PR 5 merges
- Optional intermediate baselines (if a PR introduces known non-BFB change):
  - Never needed if BFB gate is strictly enforced

**Policy:** Baselines are never overwritten. If regeneration is needed, append `-r2`, `-r3`, etc.

---

## Test Execution Guide

### Test Layers

EAMxx uses a tiered testing system:

| Tier | Description | Runtime | When to run |
|------|-------------|---------|-------------|
| 0 | Header-only, CMake config, unit tests | seconds-minutes | Every commit |
| 1 | Single-process component tests, DP integration | minutes | Before PR opens |
| 2 | Multi-process integration, full atmosphere | hours | Before PR merges |

### Running All Tests (Tier 2 validation)

**From EAMxx root:**

```bash
cd components/eamxx

# Run full test suite with baseline comparison
./scripts/test-all-eamxx -m <MACHINE> \
  --baseline-dir /data/baselines/pre-group1 -c
```

This runs ~50 tests covering:
- Build types: `sp` (single precision), `dbg` (debug), `fpe` (floating point exceptions)
- Single-process tests for each parameterization
- Multi-process integration tests
- Doubly-periodic dynamics tests

**Expected runtime:** 2-4 hours on 4 nodes.

### Running Single-Process Tests (Faster iteration)

**Test a specific parameterization:**

```bash
cd components/eamxx/scripts

# P3 microphysics (for PRs 3-4)
./create-test-suite.sh p3_standalone

cd ../ctest-build/<build_type>
ctest -R p3_standalone -VV

# SHOC turbulence (for PRs 2-3)
./create-test-suite.sh shoc_standalone
cd ../ctest-build/<build_type>
ctest -R shoc_standalone -VV
```

**Test list for Group 1:**

Priority tests (run these every PR):
- `p3_standalone` — PR 3, 4, 5 (cloud and precip)
- `shoc_standalone` — PR 2, 3 (qv and cloud)
- `rrtmgp_standalone` — PR 2, 3 (radiation needs qv, qc, qi)
- `surface_coupling` — PR 2 (qv exchange with coupler)

### Running Tracer Validation Test (PR 5)

```bash
cd components/eamxx/tests/water_tracers

# Run tracer scaling test
./run_tracer_scaling_test.sh

# Expected output:
# Validating tracer ratio preservation (expected ratio = 0.5)
# Tolerance: rtol=1e-12, atol=1e-15
#
# qv:
#   Mean ratio: 5.000000000000e-01
#   Max absolute error: 2.220e-16
#   Max relative error: 4.441e-16
#   Points checked: 3456
#   Status: PASS
#
# qc:
#   Mean ratio: 5.000000000000e-01
#   Max absolute error: 1.110e-16
#   Max relative error: 2.220e-16
#   Points checked: 2134
#   Status: PASS
# ...
# All reservoirs passed tracer ratio validation
```

### Interpreting Test Results

**BFB pass:**

```
Test p3_standalone: PASS (BFB)
  All fields match baseline within rtol=0, atol=0
```

**BFB fail:**

```
Test p3_standalone: FAIL (field mismatch)
  Field qc differs:
    Max abs diff: 1.2e-14
    Location: (col=45, lev=23)
```

**Action on BFB failure during Group 1:**

This is a **blocking failure**. Do not proceed:

1. Verify `SCREAM_NUM_TRACERS=1` is set correctly
2. Check that all slot-0 accesses use index 0 explicitly
3. Verify Kokkos subview does not change memory layout
4. Use debugger to find where bit-pattern diverges
5. If BFB is truly impossible, escalate to pause campaign and update design doc per PR 1 gate

**Tracer ratio test pass:**

All reservoirs have ratios preserved to numerical precision (max relative error < 1e-12).

**Tracer ratio test fail:**

```
qr:
  Mean ratio: 5.000142378492e-01
  Max absolute error: 1.424e-04
  Max relative error: 2.848e-04
  Points checked: 1234
  Status: FAIL
  Worst violation at index (45, 23):
    bulk = 1.234567890123e-05
    tracer = 6.190394561234e-06
    ratio = 5.014237849231e-01
```

**Action on tracer ratio failure:**

Ratio errors > 1e-12 indicate a bug in the implementation, not just numerical noise:

1. **Check initialization:** Verify tracer 1 is initialized to exactly 0.5 * tracer 0
2. **Check surface fluxes:** Verify surface flux scaling is applied correctly to slot 1
3. **Check passive transport:** Verify all physics processes apply the same rates to slots 0 and 1
4. **Check for tracer corruption:** Look for code that:
   - Inadvertently resets slot 1 to zero
   - Copies slot 0 to slot 1
   - Applies different physics to slot 1 vs slot 0
   - Uses wrong slot indexing
5. **Check sedimentation:** Verify sedimentation handles tracer dimension correctly
6. **Check phase changes:** Verify condensation/evaporation applies same rate to all slots

Common bugs:
- Hardcoded `qv(icol, ilev)` instead of `qv(0, icol, ilev)` (missing slot index)
- Loop over `itracer` but using hardcoded index 0 inside loop
- Conditional logic that only processes slot 0
- Field copies that don't preserve tracer dimension

### Performance Testing

**Measure runtime overhead:**

```bash
cd components/eamxx

# Run baseline (pre-campaign SHA)
git checkout <pre-campaign-sha>
./scripts/test-all-eamxx -m <MACHINE> -t sp --timing
# Record median runtime

# Run PR branch
git checkout wiso-02-extend-qv
./scripts/test-all-eamxx -m <MACHINE> -t sp --timing
# Compare runtimes

# Acceptable: <2% increase
```

---

## Skills Referenced

Each PR should reference the following EAMxx-specific skills during implementation:

- **`eamxx-cpp-conventions`** — C++17, Kokkos, EKAT patterns
- **`e3sm-build-and-test`** — CMake, test harness, cprnc usage
- **`e3sm-platforms`** — docker-local configuration, HPC machines
- **`scientific-modeling-conventions`** — provenance, units, field metadata
- **`regression-baseline`** — BFB comparison, baseline management
- **`tracer-conservation`** — invariants and tolerance conventions (read-only for Group 1)

PR 1 additionally references:
- **`spec-schema`** — design document structure

---

## Risks and Mitigations

### Risk 1: BFB preservation proves impossible

**Likelihood:** Low-Medium (depends on Kokkos layout flexibility)  
**Impact:** High (blocks campaign)  
**Mitigation:**
- PR 1 feasibility gate catches this early
- Fallback: Template specialization for `SCREAM_NUM_TRACERS=1` to use scalar path
- Ultimate fallback: Relax BFB requirement to tight tolerances (rtol=1e-14), document justification, get human approval

### Risk 2: Performance overhead exceeds 2%

**Likelihood:** Medium (memory layout changes can impact cache)  
**Impact:** Medium (affects acceptance into master)  
**Mitigation:**
- PR 1 prototype measures this upfront
- Fallback: Template specialization for `SCREAM_NUM_TRACERS=1` case
- Profile hot loops, consider SoA ↔ AoS layout changes
- Worst case: Accept up to 5% overhead if scientifically justified and approved

### Risk 3: Field manager doesn't support new dimension easily

**Likelihood:** Low (field manager designed for flexibility)  
**Impact:** High (PR 2 blocked)  
**Mitigation:**
- Study field manager API in PR 1 phase
- Consult EAMxx core team if `FieldLayout` extension is non-trivial
- Worst case: Extend field manager as part of PR 2 (increases PR size, may need to split PR 2 into 2a/2b)

### Risk 4: Restart/history I/O breaks with new dimension

**Likelihood:** Medium (I/O is complex)  
**Impact:** High (breaks production runs)  
**Mitigation:**
- Test restart roundtrip in PR 2: write restart, read back, verify BFB
- History output: Verify tracer dimension appears correctly in NetCDF
- Use `ncdump` to inspect output structure
- Fallback: If I/O complexity is high, defer full restart support to separate PR between PR 4 and PR 5

### Risk 5: Tracer ratio validation test fails

**Likelihood:** Medium (many code paths to check for passive tracer handling)  
**Impact:** Medium (PR 5 blocked, indicates bug in PRs 2-4)  
**Mitigation:**
- Start with simpler test: single timestep, verify ratios preserved after one physics step
- If test fails, bisect to find which process breaks ratio preservation
- Use debugger to track where ratio diverges
- Common issues: missing tracer loops, hardcoded slot-0 access, incorrect subview usage
- This is a **bug finder**, not a physics validation - ratio should be exact for passive tracers

### Risk 6: Merge conflicts with concurrent master changes

**Likelihood:** High (EAMxx actively developed)  
**Impact:** Medium (merge effort)  
**Mitigation:**
- Rebase `wiso-dev` onto `master` monthly (or after major EAMxx releases)
- Cascade rebases through PR stack: rebase `wiso-01`, then `wiso-02`, etc.
- Use merge instead of rebase if conflict resolution is too complex
- Keep PRs moving quickly to minimize divergence time

---

## Success Criteria for Group 1 Completion

Group 1 is **complete** when:

- [ ] All 5 PRs merged into `wiso-dev`
- [ ] Group 1 boundary baseline generated with `SCREAM_NUM_TRACERS=1`
- [ ] Group 1 boundary baseline is BFB vs pre-campaign baseline (rtol=0, atol=0)
- [ ] Performance overhead <2% for `SCREAM_NUM_TRACERS=1` vs pre-campaign
- [ ] Performance overhead <5% for `SCREAM_NUM_TRACERS=2` vs `SCREAM_NUM_TRACERS=1`
- [ ] `add_water_tracer(NAME test LONGNAME "Test Tracer")` compiles and runs
- [ ] Tracer ratio validation test passes: all ratios preserved to within rtol=1e-12 (double precision) or 1e-6 (single precision)
- [ ] `docs/wiso/tracer_data_model.md` exists and approved
- [ ] No regressions in any existing test
- [ ] Final PR `wiso-dev → master` opened and reviewed

---

## Post-Group-1: What Comes Next?

Group 1 delivers the **structural foundation** for water tracers. All water arrays now support multiple tracers, and passive tracer transport is verified.

Future campaigns build on this foundation:

- **Group 2:** Fractionation primitives (equilibrium, kinetic, net alpha functions)
- **Group 3:** Parameterization hooks (phase-change fractionation in SHOC, P3, etc.)
- **Group 4:** Tagged tracers (region-tagged, spherical harmonic decomp, etc.)
- **Group 5:** Production support (IC, restart, diagnostics)

Group 1 must be rock-solid because every future PR depends on this array structure and the tracer transport machinery.

---

## Appendix A: File Inventory (Estimated)

This is a rough guide; actual file counts determined during PR development.

### PR 1: ~12 files, 600 lines
- `src/physics/water_tracers/water_tracer_metadata.hpp`
- `src/physics/water_tracers/water_tracer_metadata.cpp`
- `src/physics/water_tracers/water_tracer_registry.hpp`
- `src/physics/water_tracers/water_tracer_registry.cpp`
- `src/physics/water_tracers/CMakeLists.txt`
- `src/physics/water_tracers/water_tracer_config.hpp.in`
- `src/physics/water_tracers/prototype/qv_extension_test.cpp`
- `cmake/add_water_tracer.cmake`
- `docs/wiso/tracer_data_model.md`
- `tests/water_tracers/CMakeLists.txt`
- Updates to `src/physics/CMakeLists.txt`
- Updates to top-level `CMakeLists.txt`

### PR 2: ~35 files, 900 lines
- Field layout: 3 files (`field_layout.hpp`, `field_layout.cpp`, `field_utils.cpp`)
- Field manager: 4 files (`field_manager.hpp`, `field_manager.cpp`, I/O updates)
- Process interfaces: 6 files (SHOC, P3, RRTMGP, ZM, dynamics, surface coupling)
- Kernel implementations: 15 files (all files that read/write `qv`)
- Tests: 5 files (unit tests, component tests, test harness)
- CMake: 2 files (test configuration, build updates)

### PR 3: ~40 files, 1000 lines
- Similar to PR 2, but more processes touch `qc`/`qi` (radiation, aerosols, cloud fraction)
- Additional files for `qm` if present

### PR 4: ~35 files, 800 lines
- Mostly P3 microphysics changes (precip-specific processes)
- Sedimentation code updates

### PR 5: ~15 files, 500 lines
- Tracer ratio utility: 2 files
- Test case: 3 files (C++ test, shell script, Python validation)
- Test data: 2 files (configuration, expected output)
- Documentation: 2 files (tracer_data_model.md update, test guide)
- CMake: 2 files (test integration)
- Miscellaneous: 4 files

**Total estimated: ~140 files changed across Group 1.**

---

## Appendix B: CMake Example Usage

After Group 1 completes, users can register water tracers like this:

### Example 1: Single tracer (Group 1 test case)

```cmake
# In a test CMakeLists.txt or project configuration

add_water_tracer(
  NAME test_tracer
  LONGNAME "Test Evaporation Tracer"
)

# This sets SCREAM_NUM_TRACERS=2 (slot 0=bulk, slot 1=test_tracer)
# KIND defaults to "evaporation" for Group 1
# Units are implicitly kg/kg
# Test harness will scale the evaporation flux for this tracer
```

### Example 2: Multiple tracers (future isotope campaign)

```cmake
add_water_tracer(
  NAME HDO
  LONGNAME "Hydrogen Deuterium Oxide (HDO)"
  KIND isotope
  RATIO_STANDARD VSMOW
)

add_water_tracer(
  NAME H2O18
  LONGNAME "H2^18O Water Isotope"
  KIND isotope
  RATIO_STANDARD VSMOW
)

# This sets SCREAM_NUM_TRACERS=3 (slot 0=bulk, slot 1=HDO, slot 2=H2O18)
# KIND, RATIO_STANDARD are optional and used by future campaigns
# See wiso_campaign_plan_revision1.md for complete tracer taxonomy
# Units remain kg/kg for all
```

### Generated code

The `add_water_tracer` function generates `water_tracer_config.hpp`:

```cpp
#ifndef WATER_TRACER_CONFIG_HPP
#define WATER_TRACER_CONFIG_HPP

namespace scream {
namespace water_tracers {

constexpr int NUM_TRACERS = 3;

constexpr const char* TRACER_NAMES[] = {
  "bulk_H2O",
  "HDO",
  "H2O18"
};

// Additional metadata arrays...

} // namespace water_tracers
} // namespace scream

#endif
```

---

## Appendix C: Docker Testing Guide

If using the `rfiorella/model-containers:e3sm-openmpi-dev-latest` container:

### Launch container

```bash
# Mount source directory
docker run -it --rm \
  -v /path/to/EAMXX-wiso:/code/E3SM/EAMXX-wiso \
  -v /path/to/baselines:/data/baselines \
  rfiorella/model-containers:e3sm-openmpi-dev-latest bash
```

### Inside container

```bash
cd /code/E3SM/EAMXX-wiso/components/eamxx

# Machine config auto-detected as docker-local
./scripts/test-all-eamxx -m docker-local -t sp

# Baselines stored in /data/baselines (persists across container restarts)
```

### Generate pre-campaign baseline

```bash
# From master branch, before starting Group 1
git checkout master
cd components/eamxx
./scripts/test-all-eamxx -m docker-local -t sp \
  --baseline-dir /data/baselines/pre-group1 -g
```

### Test a PR against baseline

```bash
# From PR branch
git checkout wiso-02-extend-qv
cd components/eamxx
./scripts/test-all-eamxx -m docker-local -t sp \
  --baseline-dir /data/baselines/pre-group1 -c
```

Expected output:
```
Test p3_standalone: PASS (BFB)
Test shoc_standalone: PASS (BFB)
...
All tests passed.
```

---

## Contact and Escalation

If you encounter issues during Group 1:

1. **BFB failure with `SCREAM_NUM_TRACERS=1`:** This is a blocking failure. Pause campaign, diagnose root cause, escalate if needed. Do not proceed with subsequent PRs until resolved.

2. **Performance regression >2%:** Profile hot loops, check memory layout. If >5%, escalate for architectural review.

3. **Tracer ratio test failure:** Diagnose which process is breaking ratio. May need physics team consultation.

4. **Technical questions:** Consult EAMxx core team or reference this campaign document.

**Campaign lead:** (Insert your name/contact here)

---

## Document Revision History

- **2026-05-28:** Initial Group 1 campaign revision created, following revision1 schema
- Subsequent revisions tracked in git history
