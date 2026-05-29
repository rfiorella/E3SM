# EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension

## Overview

This campaign extends EAMxx water reservoir fields to support multiple tracer
constituents along a new leading dimension. At completion, all water mass
fields (`qv`, `qc`, `qi`, `qr`, `qs`, `qm`, `rain`, `snow`) have a tracer
dimension, and users can register water tracers via `add_water_tracer()` in
CMake.

**Critical success criterion:** Existing EAMxx regression tests pass bit-for-bit
(BFB) with `SCREAM_WATER_TRACERS=0` (bulk water only in slot 0).

**Campaign scope:** 4 PRs, structural changes only, no isotope physics.

**Target timeline:** 6-8 weeks from PR 1 approval to Group 1 boundary baseline.

## Campaign Artifacts

This campaign produces:

1. **Types and metadata system** (`scream::water_tracers` namespace)
2. **CMake registration function** (`add_water_tracer`)
3. **Extended field dimensions** for 8 water reservoirs
4. **Group 1 boundary baseline** (`SCREAM_WATER_TRACERS=0`, BFB with pre-campaign)
5. **Test execution guide** for validating BFB preservation
6. **Design document** (`docs/wiso/tracer_data_model.md`)

## Design Constraints (Binding)

### Tracer dimension

- **Compile-time constant:** `SCREAM_WATER_TRACERS` (default 0)
  - Defines the number of *additional* water tracers beyond bulk water
  - Total slots allocated: `SCREAM_WATER_TRACERS + 1` (slot 0 = bulk, slots 1+ = tracers)
  - Maximum allowed value: 7 (for total of 8 slots)
- **Dimension ordering:** Leading tracer dimension
  - Example: `qv(tracer, col, lev)` instead of scalar `qv(col, lev)`
  - Rationale: Kokkos vectorization patterns favor column-innermost
  - Dimension size: Always `SCREAM_WATER_TRACERS + 1`, regardless of how many tracers are registered
  
### Performance gate

- **Acceptable overhead:** `SCREAM_WATER_TRACERS=0` (bulk-only mode) vs scalar baseline < 2% runtime
- **Measurement:** Full atmosphere component timestep, 10-run median
- **Gate location:** PR 1 prototype phase
- **Fallback:** If >2%, PR 1 must document architectural alternative before PR 2

### BFB preservation requirement

All existing EAMxx tests must pass BFB when comparing:
- **Baseline:** Pre-campaign SHA on `master` (scalar water fields)
- **Test:** Group 1 boundary SHA (tracer dimension present, `SCREAM_WATER_TRACERS=0`, only bulk water in slot 0)

This is a **hard gate**. Group 1 cannot complete until this is proven.

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

**Out of scope:** Cloud fraction, turbulent kinetic energy, and other
diagnostic fields that are not water mass.

### Slot-0 semantics

- **Slot 0 = canonical bulk water** (always present)
- When `SCREAM_WATER_TRACERS=0`, only slot 0 exists (bulk water only)
- When `SCREAM_WATER_TRACERS > 0`, slots 1, 2, ... hold additional tracers (isotopes, tagged water, etc.)
- Existing physics code reads/writes slot 0 only
- Bulk water is **never** reconstructed by summing tracer slots
- Future isotope campaigns will populate slots 1+ with HDO, H218O, etc.

### CMake integration: `add_water_tracer`

Signature (conceptual, PR 1 delivers exact API):

```cmake
add_water_tracer(
  NAME <tracer_name>
  LONGNAME <descriptive_name>
  UNITS <kg/kg or equivalent>
  RESERVOIRS qv qc qi qr qs qm rain snow
)
```

Effect:
- Registers tracer metadata in `scream::water_tracers::registry`
- Increments `SCREAM_WATER_TRACERS` at configure time
- Ensures field manager allocates sufficient slots in each reservoir
- Example: Calling `add_water_tracer` twice sets `SCREAM_WATER_TRACERS=2`, allocating 3 total slots (0=bulk, 1=first tracer, 2=second tracer)

Does **not** automatically add fractionation physics (that's future campaigns).

## Pre-Campaign Bootstrap

Before PR 1 opens, the following must exist:

### 1. `wiso-group1` long-lived branch

Created off `master` at campaign start:

```bash
git checkout master
git pull origin master
git checkout -b wiso-group1
git push -u origin wiso-group1
```

All Group 1 PRs merge into `wiso-group1`. After Group 1 completes and
baselines pass, `wiso-group1` is proposed as a PR into `master`.

### 2. Pre-campaign Tier-2 baseline

Generate baseline at the SHA where `wiso-group1` branches from `master`.

**Procedure:**

```bash
cd components/eamxx
git checkout wiso-group1

# Create baseline case
./scripts/test-all-eamxx -m <MACHINE> -c <COMPSET> -r <RESOLUTION> \
  --baseline-dir /data/baselines/pre-group1 -g

# This generates baseline without comparing
# Record SHA, date, and compset in /data/baselines/pre-group1/BASELINE.txt
```

Required baseline configurations:
- Machine: `docker-local` (or your HPC machine)
- Compset: `F2010-SCREAMv1`
- Resolution: `ne4pg2_ne4pg2`
- Test types: `dbg`, `sp`, `fpe`

**Deliverable:** `/data/baselines/pre-group1/` with ~50 test outputs.

### 3. Empty `water_tracers` directory structure

```bash
mkdir -p components/eamxx/src/physics/water_tracers
mkdir -p components/eamxx/tests/water_tracers
mkdir -p components/eamxx/docs/wiso
```

Commit this structure to `wiso-group1` before PR 1.

## PR Structure

### PR 1: Water-tracer metadata, types, and BFB feasibility gate

**Branch:** `wiso-01-tracer-metadata-and-gate`  
**Base:** `wiso-group1`  
**Tier:** 0 (no run required for merge, but prototype must run)  
**Estimated size:** 8 files, 400 lines

#### Deliverables

1. **Types and enums** (`water_tracers/water_tracer_metadata.hpp`):
   - `enum class WaterTracerKind { bulk, isotope, tagged_partition, tagged_diagnostic, auxiliary }`
   - `struct WaterTracerMetadata` with fields:
     - `std::string name`
     - `WaterTracerKind kind`
     - `std::vector<std::string> reservoirs` (subset of `qv,qc,qi,qr,qs,qm,rain,snow`)
     - `std::string units`
     - `bool conserved_independently`
   - `class WaterTracerRegistry` (singleton, compile-time registration)

2. **CMake function** (`water_tracers/CMakeLists.txt`):
   - `add_water_tracer()` macro
   - Sets `SCREAM_WATER_TRACERS` (default 0, meaning bulk-only)
   - Generates `water_tracer_config.hpp` with tracer count and metadata

3. **Design document** (`docs/wiso/tracer_data_model.md`):
   - Slot-0 semantics (canonical bulk, never reconstructed from sum)
   - `SCREAM_WATER_TRACERS` meaning: count of additional tracers, not including bulk
   - Array indexing: slot 0 = bulk, slot 1+ = additional tracers
   - BFB feasibility result (see gate below)
   - Field layout decision (leading tracer dim)
   - Performance impact (<2% or fallback justification)
   - H216O implicit/explicit decision (deferred to isotope campaign)

4. **Prototype** (`water_tracers/prototype/qv_extension_test.cpp`, not merged):
   - Minimal Kokkos kernel extending `qv(col,lev)` → `qv(tracer,col,lev)`
   - Benchmarks `SCREAM_WATER_TRACERS=0` (single slot) vs scalar equivalent
   - Measures:
     - Runtime (target <2% overhead)
     - Memory bandwidth
     - Vectorization efficiency (check assembly if needed)

#### BFB Feasibility Gate

PR 1 includes a **mandatory checkpoint** before PRs 2-4 can start.

**Task:** Run the `qv_extension_test` prototype on target hardware and measure:

1. **Performance:** Compare `SCREAM_WATER_TRACERS=0` multidimensional access vs scalar
   - Pass: <2% overhead
   - Fail: Document why, propose alternative (e.g., template specialization for scalar path)

2. **BFB preservation:** Compare prototype output with known-good scalar qv kernel
   - Pass: Bit-identical
   - Fail: Document source of difference (packing? ordering? FP contraction?)

**Decision point:**

- If **both pass**: PRs 2-4 proceed with current design, BFB is a hard gate
- If **performance fails**: Update `tracer_data_model.md` with fallback (e.g., `if constexpr (SCREAM_WATER_TRACERS==0)` scalar path), then proceed
- If **BFB fails**: Pause campaign, diagnose root cause, update design doc, then proceed with relaxed BFB tolerance (requires approval)

**Approval:** A human must sign off on the design document before PR 1 merges.

#### Success Criteria

- [ ] `water_tracer_metadata.hpp` compiles without warnings
- [ ] `add_water_tracer(NAME test LONGNAME "Test" UNITS kg/kg RESERVOIRS qv)` generates valid config header
- [ ] `qv_extension_test` prototype compiles and runs
- [ ] Prototype shows <2% overhead OR `tracer_data_model.md` documents approved fallback
- [ ] Prototype output is BFB vs scalar OR `tracer_data_model.md` documents approved relaxation
- [ ] `docs/wiso/tracer_data_model.md` reviewed and approved by human
- [ ] No existing tests are modified (this PR is headers/docs only for production code)

#### Out of Scope

- Actual field registration changes (that's PR 2)
- Isotope-specific logic (future campaign)
- Fractionation functions (future campaign)
- Test harness for multi-tracer mode (PR 2+)

---

### PR 2: Extend `qv` to tracer dimension

**Branch:** `wiso-02-extend-qv`  
**Base:** `wiso-01-tracer-metadata-and-gate`  
**Tier:** 2 (full Tier-2 validation required)  
**Estimated size:** 30 files, 800 lines

#### Dependencies

- PR 1 merged with approved `tracer_data_model.md`
- BFB feasibility gate passed

#### Deliverables

1. **Field registration** in all processes that touch `qv`:
   - Modify `add_field<Required>("qv", ...)` calls to use tracer-aware layout
   - Update layout from `scalar3d_mid` to `tracer3d_mid` (new layout type)
   - Example locations:
     - `physics/shoc/eamxx_shoc_process_interface.cpp`
     - `physics/p3/eamxx_p3_process_interface.cpp`
     - `physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp`
     - `physics/zm/eamxx_zm_process_interface.cpp`

2. **Kernel updates** to access `qv` with tracer slice:
   - Pattern: `qv(icol, ilev)` → `qv(0, icol, ilev)` for all slot-0 bulk access
   - Use Kokkos subview where entire tracer dimension is needed:
     ```cpp
     auto qv_bulk = Kokkos::subview(qv, 0, Kokkos::ALL, Kokkos::ALL);
     ```
   - Affected files: ~20 files across `physics/`, estimated from grep results

3. **Field manager integration**:
   - Extend `FieldLayout` to support `TRACER` dimension tag
   - Update field allocation to respect `SCREAM_WATER_TRACERS + 1` (for total slot count)
   - Update restart/history I/O to handle tracer dimension (write all slots, even when `SCREAM_WATER_TRACERS=0`)

4. **Unit tests**:
   - `tests/water_tracers/test_qv_tracer_access.cpp`
     - Verifies slot-0 access is equivalent to old scalar access
     - Verifies tracer dimension is correctly sized (`SCREAM_WATER_TRACERS + 1`)
     - When `SCREAM_WATER_TRACERS=0`, verifies only 1 slot exists

#### Success Criteria

- [ ] All files compile without warnings
- [ ] `add_water_tracer(NAME bulk RESERVOIRS qv)` successfully registers
- [ ] Standalone `qv` unit test passes
- [ ] **BFB requirement:** All existing EAMxx tests pass BFB vs pre-campaign baseline
  - Verified via `cprnc` on all 3D/2D output fields
  - Zero tolerance: rtol=0, atol=0
  - Test list: see "Test Execution Guide" section below
- [ ] Performance regression: <2% runtime increase on full atmosphere timestep
- [ ] No segfaults, no MPI hangs, no NaN/Inf in output

#### Out of Scope

- Other water species (qc, qi, qr, qs, qm, rain, snow) — those are PRs 3-4
- Multi-tracer registration (`SCREAM_WATER_TRACERS > 0`) — validated in later campaign
- Isotope-specific kernels

---

### PR 3: Extend `qc`, `qi`, `qm` to tracer dimension

**Branch:** `wiso-03-extend-cloud`  
**Base:** `wiso-02-extend-qv`  
**Tier:** 2  
**Estimated size:** 35 files, 900 lines

#### Dependencies

- PR 2 merged
- PR 2 BFB validation passed

#### Deliverables

Same pattern as PR 2, applied to cloud water species:

1. **Field registration** for `qc`, `qi`, `qm`:
   - `physics/p3/` (primary microphysics)
   - `physics/shoc/` (cloud fraction, turbulence)
   - `physics/cld_fraction/`
   - `physics/rrtmgp/` (radiation needs cloud water)
   - `physics/mam/` (aerosol activation)

2. **Kernel updates**:
   - All `qc(icol, ilev)` → `qc(0, icol, ilev)`
   - All `qi(icol, ilev)` → `qi(0, icol, ilev)`
   - All `qm(icol, ilev)` → `qm(0, icol, ilev)` if qm exists

3. **Derived fields** that compute from qc/qi:
   - `inv_qc_relvar` (cloud variance)
   - `eff_radius_qc`, `eff_radius_qi` (effective radius)
   - These may remain scalar (no tracer dim) if they're diagnostic-only

4. **Unit tests**:
   - `tests/water_tracers/test_cloud_tracer_access.cpp`

#### Success Criteria

- [ ] Compiles without warnings
- [ ] `add_water_tracer(NAME bulk RESERVOIRS qv qc qi qm)` registers successfully
- [ ] Unit tests pass
- [ ] **BFB requirement:** All existing tests pass BFB vs pre-campaign baseline
- [ ] Performance: <2% cumulative overhead (PR 2 + PR 3)
- [ ] P3 microphysics single-process test passes BFB
- [ ] SHOC single-process test passes BFB

#### Out of Scope

- Precipitation species (qr, qs, rain, snow) — that's PR 4
- Phase-change fractionation logic

---

### PR 4: Extend `qr`, `qs`, `rain`, `snow` to tracer dimension

**Branch:** `wiso-04-extend-precip`  
**Base:** `wiso-03-extend-cloud`  
**Tier:** 2  
**Estimated size:** 30 files, 700 lines

#### Dependencies

- PR 3 merged
- PR 3 BFB validation passed

#### Deliverables

Same pattern as PRs 2-3, applied to precipitation:

1. **Field registration** for `qr`, `qs`, `rain`, `snow`:
   - `physics/p3/` (handles rain/snow microphysics)
   - Any sedimentation/transport code

2. **Kernel updates**:
   - `qr(icol, ilev)` → `qr(0, icol, ilev)`
   - `qs(icol, ilev)` → `qs(0, icol, ilev)`
   - `rain(icol, ilev)` → `rain(0, icol, ilev)`
   - `snow(icol, ilev)` → `snow(0, icol, ilev)`

3. **P3 process-specific updates**:
   - Sedimentation must handle tracer dimension
   - Autoconversion, collection, etc. remain bulk-only (slot 0)

4. **Unit tests**:
   - `tests/water_tracers/test_precip_tracer_access.cpp`

#### Success Criteria

- [ ] Compiles without warnings
- [ ] `add_water_tracer(NAME bulk RESERVOIRS qv qc qi qr qs qm rain snow)` successful
- [ ] Unit tests pass
- [ ] **BFB requirement:** All existing tests pass BFB vs pre-campaign baseline
- [ ] Performance: <2% cumulative overhead (PRs 2+3+4)
- [ ] P3 microphysics with precipitation passes BFB
- [ ] Full atmosphere multi-process test passes BFB

#### Out of Scope

- Convective precipitation separate from stratiform (decide in PR 4 if needed)
- Isotope-aware sedimentation (future campaign)

---

## Group 1 Boundary: Baseline Regeneration

After PR 4 merges into `wiso-group1`, generate the **Group 1 boundary baseline**:

```bash
cd components/eamxx
git checkout wiso-group1
git pull origin wiso-group1

# Generate new baseline with tracer dimension present, SCREAM_WATER_TRACERS=0
./scripts/test-all-eamxx -m <MACHINE> -c <COMPSET> -r <RESOLUTION> \
  --baseline-dir /data/baselines/group1-boundary -g

# Record SHA, date, SCREAM_WATER_TRACERS=0 in BASELINE.txt
```

**Validation:**

Compare Group 1 boundary baseline against pre-campaign baseline:

```bash
./scripts/compare-baselines.sh \
  /data/baselines/pre-group1 \
  /data/baselines/group1-boundary \
  --rtol 0 --atol 0
```

Expected result: **BFB identical** for all fields in all tests.

If not BFB:
1. Identify which field(s) differ
2. Bisect PRs 2-4 to find introduction point
3. Fix and regenerate baseline
4. Repeat until BFB

**Group 1 is not complete until this passes.**

---

## Test Execution Guide

This section provides commands to run existing EAMxx tests and verify BFB
preservation during Group 1 development.

### Prerequisites

- Working EAMxx build environment
- Access to compute resources (tests require ~1-4 hours total)
- Baselines generated as described in "Pre-Campaign Bootstrap"

### Test Layers

EAMxx uses a tiered testing system:

| Tier | Description | Runtime | When to run |
|------|-------------|---------|-------------|
| 0 | Header-only, CMake config | seconds | Every commit |
| 1 | Unit tests (Catch2), single-process | minutes | Before PR opens |
| 2 | Multi-process integration, full atmosphere | hours | Before PR merges |

### Running All Tests (Tier 2 validation)

**From EAMxx root:**

```bash
cd components/eamxx

# Run full test suite
./scripts/test-all-eamxx -m <MACHINE> -t sp,dbg,fpe
```

This runs ~50 tests covering:
- Build types: `sp` (single precision), `dbg` (debug), `fpe` (floating point exceptions)
- Single-process tests for each parameterization
- Multi-process integration tests
- Doubly-periodic dynamics tests

**Expected runtime:** 2-4 hours on 4 nodes.

### Running Tests with Baseline Comparison

**During PR development:**

```bash
# Compare against pre-campaign baseline
./scripts/test-all-eamxx -m <MACHINE> -t sp \
  --baseline-dir /data/baselines/pre-group1 -c

# The -c flag compares; -g flag generates
```

**Output:** Pass/fail for each test. Failures report which fields differ.

### Running Single-Process Tests (Faster iteration)

**Test a specific parameterization:**

```bash
cd components/eamxx/scripts

# P3 microphysics
./create-test-suite.sh p3_standalone

cd ../ctest-build/<build_type>
ctest -R p3_standalone -VV

# SHOC turbulence
./create-test-suite.sh shoc_standalone
cd ../ctest-build/<build_type>
ctest -R shoc_standalone -VV
```

**Test list for Group 1:**

Priority tests (run these every PR):
- `p3_standalone` — PR 3 and 4 (cloud and precip)
- `shoc_standalone` — PR 2 and 3 (qv and cloud)
- `rrtmgp_standalone` — PR 2 and 3 (radiation needs qv, qc, qi)
- `zm_standalone` — PR 2 (convection needs qv)

Full test list (run before merge):
- All above
- `homme_shoc_cld_spa_p3_rrtmgp` (full dynamics-physics)
- `surface_coupling` (PR 2, qv exchange with coupler)

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

**Action on failure:**

1. Check if difference is in a field you modified (expected during development)
2. If in an unrelated field, investigate potential side effects
3. Use `cprnc` for detailed field-by-field comparison:

```bash
cprnc -m baseline_output.nc test_output.nc
```

### Performance Testing

**Measure runtime overhead:**

```bash
cd components/eamxx

# Run baseline (pre-campaign SHA)
git checkout <pre-campaign-sha>
./scripts/test-all-eamxx -m <MACHINE> -t sp --timing

# Record median runtime, store in baseline_timing.txt

# Run PR branch
git checkout wiso-02-extend-qv
./scripts/test-all-eamxx -m <MACHINE> -t sp --timing

# Compare runtimes
# Acceptable: <2% increase
```

**Detailed profiling (if needed):**

```bash
# Run with profiler
SCREAM_PROFILE=1 ./scripts/test-all-eamxx -m <MACHINE> -t sp

# Analyze with profiler tool (depends on machine):
# - NERSC: visit generated .json files in web browser
# - LCRC: use Intel VTune
# - Docker: use perf or gprof
```

### Docker-Specific Testing

If using the `rfiorella/model-containers:e3sm-openmpi-dev-latest` container:

```bash
# Launch container with source mounted
docker run -it --rm -v /path/to/EAMXX-wiso:/code/E3SM/EAMXX-wiso \
  rfiorella/model-containers:e3sm-openmpi-dev-latest bash

# Inside container
cd /code/E3SM/EAMXX-wiso/components/eamxx

# Machine config auto-detected as docker-local
./scripts/test-all-eamxx -m docker-local -t sp

# Baselines stored in /data/baselines (mount a volume for persistence)
```

### Continuous Integration

If your repository has GitHub Actions or similar CI:

1. Add workflow to run `test-all-eamxx` on every PR to `wiso-group1`
2. Cache baseline data as workflow artifact
3. Set `--baseline-dir` to cached baseline path
4. Fail PR if any test shows non-BFB differences

Example GitHub Actions snippet (adapt as needed):

```yaml
name: EAMxx Group 1 Tests
on:
  pull_request:
    branches: [wiso-group1]
jobs:
  test:
    runs-on: ubuntu-latest
    container: rfiorella/model-containers:e3sm-openmpi-dev-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run tests
        run: |
          cd components/eamxx
          ./scripts/test-all-eamxx -m docker-local -t sp \
            --baseline-dir /data/baselines/pre-group1 -c
```

---

## Branch and Baseline Model

### Branch structure

```
master
  └─ wiso-group1 (long-lived)
       ├─ wiso-01-tracer-metadata-and-gate
       ├─ wiso-02-extend-qv (bases on wiso-01)
       ├─ wiso-03-extend-cloud (bases on wiso-02)
       └─ wiso-04-extend-precip (bases on wiso-03)
```

### Merge strategy

1. PR 1 opens against `wiso-group1`
2. After PR 1 merges, PR 2 opens against `wiso-group1` (rebased from `wiso-01` branch)
3. After PR 2 merges, PR 3 opens against `wiso-group1` (rebased from `wiso-02` branch)
4. After PR 3 merges, PR 4 opens against `wiso-group1` (rebased from `wiso-03` branch)
5. After PR 4 merges and Group 1 boundary baseline passes, open final PR: `wiso-group1 → master`

### Baseline naming

- `pre-group1/` — generated from `master` before campaign starts
- `group1-boundary/` — generated from `wiso-group1` after PR 4 merges
- Optional intermediate baselines (if a PR introduces known non-BFB change):
  - `wiso-02-qv/` — if PR 2 needs to adjust BFB policy
  - etc.

**Policy:** Baselines are never overwritten. If regeneration is needed, append `-r2`, `-r3`, etc.

---

## PR Size Policy

Target limits:
- **Files changed:** <30 per PR
- **Lines changed:** <800 per PR

PRs 2-4 may slightly exceed these (estimated 30-35 files, 700-900 lines) due to
the pervasive nature of water field access. This is acceptable because:
- Changes are mechanical (pattern: `qv(i,j)` → `qv(0,i,j)`)
- BFB testing validates correctness
- Each PR focuses on a single water species group

If a PR exceeds 40 files or 1200 lines:
- Split into `Na` and `Nb` sub-PRs
- Example: PR 3a extends qc/qi, PR 3b extends qm

---

## Skills Referenced

Each PR spec should reference the following EAMxx-specific skills:

- **`eamxx-cpp-conventions`** — C++17, Kokkos, EKAT patterns
- **`e3sm-build-and-test`** — CMake, test harness, cprnc usage
- **`e3sm-platforms`** — docker-local configuration, HPC machines
- **`scientific-modeling-conventions`** — provenance, units, field metadata
- **`regression-baseline`** — BFB comparison, baseline management

PR 1 additionally references:
- **`spec-schema`** — design document structure
- **`tracer-conservation`** (read-only, to understand future requirements)

---

## Risks and Mitigations

### Risk 1: BFB preservation proves impossible

**Likelihood:** Low (similar changes done in CAM)  
**Impact:** High (blocks campaign)  
**Mitigation:**
- PR 1 feasibility gate catches this early
- Fallback: Relax BFB requirement to tight tolerances (rtol=1e-14), document justification
- Ultimate fallback: Scalar compatibility path via `if constexpr`

### Risk 2: Performance overhead exceeds 2%

**Likelihood:** Medium (depends on memory layout)  
**Impact:** Medium (affects acceptance into master)  
**Mitigation:**
- PR 1 prototype measures this upfront
- Fallback: Template specialization for `SCREAM_WATER_TRACERS=0` case
- Profile hot loops, consider SoA ↔ AoS layout changes

### Risk 3: Field manager doesn't support new dimension easily

**Likelihood:** Low (field manager designed for flexibility)  
**Impact:** High (PR 2 blocked)  
**Mitigation:**
- Study field manager API in PR 1 phase
- Consult EAMxx core team if `FieldLayout` extension is non-trivial
- Worst case: Extend field manager as part of PR 2 (increases PR size)

### Risk 4: Restart/history I/O breaks with new dimension

**Likelihood:** Medium (I/O is complex)  
**Impact:** High (breaks production runs)  
**Mitigation:**
- Test restart roundtrip in PR 2: write restart, read back, verify BFB
- History output: Verify tracer dimension appears correctly in NetCDF
- Fallback: Defer full restart support to later PR (out-of-scope for Group 1)

### Risk 5: Merge conflicts with concurrent master changes

**Likelihood:** High (EAMxx actively developed)  
**Impact:** Medium (merge effort)  
**Mitigation:**
- Rebase `wiso-group1` onto `master` monthly (or after major EAMxx releases)
- Cascade rebases through PR stack: rebase `wiso-01`, then `wiso-02`, etc.
- Use merge instead of rebase if conflict resolution is too complex

---

## Success Criteria for Group 1 Completion

Group 1 is **complete** when:

- [ ] All 4 PRs merged into `wiso-group1`
- [ ] Group 1 boundary baseline generated
- [ ] Group 1 boundary baseline is BFB vs pre-campaign baseline (rtol=0, atol=0)
- [ ] Performance overhead <2% for full atmosphere timestep
- [ ] `add_water_tracer(NAME test RESERVOIRS qv qc qi qr qs qm rain snow)` compiles and runs (sets `SCREAM_WATER_TRACERS=1`, allocates 2 slots)
- [ ] Building with `SCREAM_WATER_TRACERS=0` (default, bulk-only) passes all tests BFB
- [ ] `docs/wiso/tracer_data_model.md` exists and approved
- [ ] Test execution guide validated by independent user
- [ ] No regressions in any existing test
- [ ] Final PR `wiso-group1 → master` opened and reviewed

---

## Post-Group-1: What Comes Next?

Group 1 delivers the **structural foundation**. Future campaigns build on it:

- **Group 2:** Fractionation primitives (equilibrium, kinetic, net alpha functions)
- **Group 3:** Parameterization hooks (phase-change fractionation in SHOC, P3, etc.)
- **Group 4:** Tagged tracers (region-tagged, spherical harmonic decomp, etc.)

Group 1 must be rock-solid because every future PR depends on this array structure.

---

## Appendix: File Inventory (Estimated)

This is a rough guide; actual file counts determined during PR development.

### PR 1: ~8 files
- `src/physics/water_tracers/water_tracer_metadata.hpp`
- `src/physics/water_tracers/water_tracer_registry.hpp`
- `src/physics/water_tracers/water_tracer_registry.cpp`
- `src/physics/water_tracers/CMakeLists.txt`
- `src/physics/water_tracers/prototype/qv_extension_test.cpp`
- `cmake/add_water_tracer.cmake`
- `docs/wiso/tracer_data_model.md`
- `tests/water_tracers/CMakeLists.txt`

### PR 2: ~30 files
- Field registration: 5 files (`shoc`, `p3`, `rrtmgp`, `zm`, `surface_coupling`)
- Kernel updates: 15 files (across parameterizations)
- Field manager: 3 files (`FieldLayout`, allocation, I/O)
- Tests: 5 files (unit tests, test data)
- Docs: 2 files (update design doc, inline comments)

### PR 3: ~35 files
- Similar to PR 2, but more files touch `qc`/`qi` (radiation, aerosols, cloud fraction)

### PR 4: ~30 files
- Mostly P3 microphysics changes (precip-specific processes)

**Total estimated: ~100 files changed across Group 1.**

---

## Contact and Escalation

If you encounter issues during Group 1:

1. **BFB failure:** Check if it's a known limitation (e.g., vectorization order changed). If truly blocking, escalate to pause campaign and revisit design.
2. **Performance regression:** Profile hot loops, consider layout changes, escalate if >5% overhead.
3. **Technical questions:** Consult EAMxx core team or this campaign document.

**Campaign lead:** (Insert your name/contact here)

---

## Document Revision History

- **2026-MM-DD:** Initial Group 1 campaign plan created
- Subsequent revisions tracked in git history
