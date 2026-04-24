# P3 Water Isotope Implementation - Issue Templates

**Format**: GitHub Issues / Kanban-compatible JSON  
**Version**: 1.0  
**Date**: 2026-04-21

---

## Issues Overview

Total Issues: 50
- Phase 0: 1 issue (completed)
- Phase 1: 7 issues
- Phase 2: 8 issues  
- Phase 3: 6 issues
- Phase 4: 6 issues
- Phase 5: 6 issues
- Phase 6: 5 issues
- Phase 7: 5 issues
- Phase 8: 3 issues
- Phase 9: 3 issues
- Phase 10: 4 issues (optional)

---

## Issue Templates (Markdown Format)

### Phase 0: Infrastructure Setup

#### Issue #0.1: Project Planning and Specification ✓
**Status**: Completed  
**Priority**: Critical  
**Labels**: `phase-0`, `planning`, `documentation`  
**Assignees**: TBD  
**Dependencies**: None  

**Description**:
Document existing systems and create comprehensive implementation roadmap for porting water isotope infrastructure from EAMv2-wiso (MG2 microphysics) to EAMv3-wiso (P3 microphysics).

**Tasks**:
- [x] Catalog P3 microphysics functions requiring isotope support
- [x] Document EAMv2-wiso isotope infrastructure
- [x] Create detailed specification document
- [ ] Review specification with science team
- [ ] Finalize design decisions

**Acceptance Criteria**:
- Specification document approved by team
- Function mapping complete
- Phased implementation plan agreed upon

**Deliverables**:
- `P3_ISOTOPE_IMPLEMENTATION_SPEC.md`
- `P3_ISOTOPE_FUNCTION_MAPPING.md`
- Issue templates for all phases

---

### Phase 1: Core Infrastructure Port

#### Issue #1.1: Port water_isotopes.F90 Module
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-1`, `infrastructure`, `fractionation`  
**Assignees**: TBD  
**Dependencies**: #0.1  
**Estimated Effort**: 3-5 days  

**Description**:
Port the core isotope fractionation calculation module from EAMv2-wiso to EAMv3-wiso. This module contains fundamental fractionation physics for all isotope species.

**Source File**: `/Users/rfiorella/code/E3SM/EAMv2-wiso/share/util/water_isotopes.F90`  
**Destination**: `/Users/rfiorella/code/E3SM/EAMv3-wiso/share/util/water_isotopes.F90`

**Tasks**:
- [ ] Copy water_isotopes.F90 to EAMv3-wiso
- [ ] Update module dependencies for EAMv3 compatibility
- [ ] Verify all fractionation functions compile
- [ ] Test `wiso_alpl()` - liquid-vapor equilibrium fractionation
- [ ] Test `wiso_alpi()` - ice-vapor equilibrium fractionation
- [ ] Test `wiso_akel()` - kinetic fractionation (liquid)
- [ ] Test `wiso_akci()` - kinetic fractionation (ice)
- [ ] Verify isotope species constants (H2O, HDO, H218O, H217O, HTO)
- [ ] Update build system (Makefile/CMakeLists)
- [ ] Create unit tests

**Acceptance Criteria**:
- Module compiles without errors in EAMv3-wiso
- All fractionation functions return expected values
- Unit tests pass for all isotope species
- Documentation updated

**Key Functions**:
- `wiso_init()` - Initialize module
- `wiso_alpl()` - Liquid-vapor fractionation (Horita & Wesolowski 1994)
- `wiso_alpi()` - Ice-vapor fractionation (Merlivat & Nief 1967)
- `wiso_akel()` - Kinetic fractionation for evaporation
- `wiso_akci()` - Kinetic fractionation for ice
- `wiso_flxoce()` - Ocean evaporation isotopes
- `wiso_decay()` - Radioactive decay (HTO)

---

#### Issue #1.2: Port water_types.F90 Module
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-1`, `infrastructure`, `types`  
**Assignees**: TBD  
**Dependencies**: #1.1  
**Estimated Effort**: 1-2 days  

**Description**:
Port water phase type definitions module. Defines water types (vapor, liquid, ice, rain, snow) and master fractionation interface.

**Source File**: `/Users/rfiorella/code/E3SM/EAMv2-wiso/share/util/water_types.F90`  
**Destination**: `/Users/rfiorella/code/E3SM/EAMv3-wiso/share/util/water_types.F90`

**Tasks**:
- [ ] Copy water_types.F90 to EAMv3-wiso
- [ ] Verify water type constants (iwtvap, iwtliq, iwtice, etc.)
- [ ] Test `wtype_get_alpha()` - master fractionation lookup
- [ ] Update build system
- [ ] Create unit tests

**Acceptance Criteria**:
- Module compiles without errors
- All water types defined correctly
- Master fractionation interface works
- Unit tests pass

**Key Constants**:
- `iwtvap` - Water vapor type
- `iwtliq` - Cloud liquid type
- `iwtice` - Cloud ice type
- `iwtstrain` - Stratiform rain type
- `iwtstsnow` - Stratiform snow type
- `iwtcvrain` - Convective rain type
- `iwtcvsnow` - Convective snow type

---

#### Issue #1.3: Port water_tracer_vars.F90 Module
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-1`, `infrastructure`, `configuration`  
**Assignees**: TBD  
**Dependencies**: #1.2  
**Estimated Effort**: 1-2 days  

**Description**:
Port configuration variables and namelist parameters for water isotope tracers.

**Source File**: `/Users/rfiorella/code/E3SM/EAMv2-wiso/components/eam/src/physics/cam/water_tracer_vars.F90`  
**Destination**: `/Users/rfiorella/code/E3SM/EAMv3-wiso/components/eam/src/physics/cam/water_tracer_vars.F90`

**Tasks**:
- [ ] Copy water_tracer_vars.F90
- [ ] Define tracer indexing arrays
- [ ] Set up namelist parameters
- [ ] Update build system
- [ ] Document configuration options

**Acceptance Criteria**:
- Module compiles without errors
- Configuration variables accessible
- Namelist can be read

**Key Variables**:
- `wisotope` - Enable water isotopes flag
- `wtrc_alpha_kinetic` - Enable kinetic fractionation flag
- `wtrc_niter` - Number of iterations for solver
- `wtrc_names` - Tracer names array
- `wtrc_species_names` - Isotope species names

---

#### Issue #1.4: Port Core Functions from water_tracers.F90
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-1`, `infrastructure`, `core`  
**Assignees**: TBD  
**Dependencies**: #1.3  
**Estimated Effort**: 5-7 days  

**Description**:
Port core water tracer module with essential functions for initialization, rate tracking, and basic operations. This is a subset of the full module (~7,800 lines) - only core functions needed for Phase 1.

**Source File**: `/Users/rfiorella/code/E3SM/EAMv2-wiso/components/eam/src/physics/cam/water_tracers.F90`  
**Destination**: `/Users/rfiorella/code/E3SM/EAMv3-wiso/components/eam/src/physics/cam/water_tracers.F90`

**Tasks**:
- [ ] Create skeleton water_tracers.F90 module
- [ ] Port `wtrc_readnl()` - Read namelist
- [ ] Port `wtrc_init()` - Initialize module
- [ ] Port `wtrc_register()` - Register constituents
- [ ] Port `wtrc_init_cnst()` - Initialize constituent values
- [ ] Port `wtrc_init_rates()` - Initialize rate arrays
- [ ] Port `wtrc_add_rates()` - Add process rates
- [ ] Port `wtrc_get_alpha()` - Fractionation factor wrapper
- [ ] Port `wtrc_ratio()` - Calculate isotope ratios
- [ ] Port `wtrc_ratio_all()` - Calculate all tracer ratios
- [ ] Add module to build system
- [ ] Create unit tests for each function

**Acceptance Criteria**:
- Module compiles without errors
- All core functions work correctly
- Unit tests pass
- No runtime errors during initialization

**Functions Included** (Phase 1 subset):
- Initialization: `wtrc_init()`, `wtrc_register()`, `wtrc_init_cnst()`
- Rate handling: `wtrc_init_rates()`, `wtrc_add_rates()`, `wtrc_add_rate()`
- Utilities: `wtrc_get_alpha()`, `wtrc_ratio()`, `wtrc_ratio_all()`
- Queries: `wtrc_is_vap()`, `wtrc_is_liq()`, `wtrc_is_ice()`, etc.

**Functions Deferred** to later phases:
- `wtrc_apply_rates()` → Phase 2
- `stewart_isoevap()` → Phase 2
- `wtrc_sediment()` → Phase 7
- Convection functions → Phase 8

---

#### Issue #1.5: Extend physics_state for Isotope Tracers
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-1`, `infrastructure`, `state`  
**Assignees**: TBD  
**Dependencies**: #1.4  
**Estimated Effort**: 2-3 days  

**Description**:
Extend the physics state structure to include isotope tracer arrays. Ensure compatibility with existing physics infrastructure.

**File to Modify**: `components/eam/src/physics/cam/physics_types.F90`

**Tasks**:
- [ ] Review physics_state structure
- [ ] Add isotope tracer pointers/arrays
- [ ] Update physics_ptend for isotope tendencies
- [ ] Modify allocation routines
- [ ] Update initialization routines
- [ ] Test state structure changes
- [ ] Verify backward compatibility (isotopes disabled)

**Acceptance Criteria**:
- Physics state structure includes isotope tracers
- Code compiles with and without isotopes enabled
- No memory leaks
- Existing physics tests still pass

---

#### Issue #1.6: Register Isotope Constituents
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-1`, `infrastructure`, `constituents`  
**Assignees**: TBD  
**Dependencies**: #1.5  
**Estimated Effort**: 2-3 days  

**Description**:
Register all water isotope tracers as model constituents. Set up proper indexing and metadata.

**File to Modify**: `components/eam/src/physics/cam/physpkg.F90`

**Tasks**:
- [ ] Call `wtrc_register()` during constituent registration
- [ ] Set up constituent indices for each isotope species
- [ ] Configure advection for isotope tracers
- [ ] Set up vertical diffusion for isotope tracers
- [ ] Add surface flux fields
- [ ] Test constituent registration
- [ ] Verify tracer indices are correct

**Acceptance Criteria**:
- All isotope tracers registered as constituents
- Tracers properly advected
- Tracer indices accessible throughout code
- Initial conditions can be read/set

**Constituents to Register** (per water type):
- H2O (standard)
- H216O
- HDO
- H218O
- H217O
- HTO

**Water Types**:
- Vapor (always)
- Cloud liquid (always)
- Cloud ice (always)
- Stratiform rain (if enabled)
- Stratiform snow (if enabled)
- Convective rain (if enabled, Phase 8)
- Convective snow (if enabled, Phase 8)

---

#### Issue #1.7: Add Namelist Configuration for Isotopes
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-1`, `infrastructure`, `namelist`  
**Assignees**: TBD  
**Dependencies**: #1.6  
**Estimated Effort**: 1-2 days  

**Description**:
Create namelist entries for configuring water isotope tracers. Add to build system and documentation.

**Files to Modify**:
- `components/eam/bld/namelist_files/namelist_definition.xml`
- `components/eam/bld/namelist_files/namelist_defaults.xml`
- `components/eam/bld/configure`

**Tasks**:
- [ ] Add water_tracer_nl namelist group definition
- [ ] Define all isotope configuration parameters
- [ ] Set default values
- [ ] Create test namelists
- [ ] Add validation logic
- [ ] Document namelist options
- [ ] Create example use cases

**Acceptance Criteria**:
- Namelist can be read successfully
- Default configuration works
- All options documented
- Validation catches errors

**Key Namelist Parameters**:
- `wisotope` - Enable isotopes (logical)
- `wtrc_alpha_kinetic` - Enable kinetic fractionation (logical)
- `wtrc_niter` - Number of solver iterations (integer)
- `wtrc_qmin` - Minimum tracer threshold (real)
- `wtrc_add_stprecip` - Track stratiform precip (logical)
- `wtrc_add_cvprecip` - Track convective precip (logical)

---

### Phase 2: Stratiform Rain Evaporation

#### Issue #2.1: Modify evaporate_rain() for Rate Output
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-2`, `p3-modification`, `rain-evaporation`  
**Assignees**: TBD  
**Dependencies**: #1.7  
**Estimated Effort**: 2-3 days  

**Description**:
Modify P3's `evaporate_rain()` function to output evaporation rate in format compatible with isotope interface. This is the highest priority process for precipitation isotope composition.

**File to Modify**: `components/eam/src/physics/p3/eam/micro_p3.F90` (~line 2800)

**Tasks**:
- [ ] Review current `evaporate_rain()` implementation
- [ ] Add output for evaporation rate (kg/kg/s)
- [ ] Ensure rate is available per column and level
- [ ] Add diagnostic output for debugging
- [ ] Test rate calculation
- [ ] Verify mass conservation
- [ ] Document changes

**Acceptance Criteria**:
- Evaporation rate correctly calculated and output
- Rate units are kg/kg/s
- Mass balance maintained
- No change to existing P3 behavior when isotopes disabled
- Code compiles and runs

**Output Added**:
```fortran
real(rtype), intent(out) :: qr_evap_rate(ncol, nlev)  ! Rain evaporation rate (kg/kg/s)
```

---

#### Issue #2.2: Create wtrc_p3_inter_phase2() Interface
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-2`, `interface`, `new-code`  
**Assignees**: TBD  
**Dependencies**: #2.1  
**Estimated Effort**: 3-4 days  

**Description**:
Create simplified P3 isotope interface for Phase 2, handling only rain evaporation. This is modeled after `wtrc_mg_inter()` from EAMv2-wiso but simplified.

**New File**: `components/eam/src/physics/p3/eam/water_tracers_p3.F90`

**Tasks**:
- [ ] Create new module skeleton
- [ ] Define interface subroutine signature
- [ ] Add process rate collection (rain evaporation only)
- [ ] Initialize rate arrays using `wtrc_init_rates()`
- [ ] Call `wtrc_add_rates()` for evaporation
- [ ] Call `wtrc_apply_rates()` (simplified version)
- [ ] Pass temperature and cloud fraction
- [ ] Add to build system
- [ ] Create unit tests
- [ ] Document interface

**Acceptance Criteria**:
- Interface function compiles
- Can collect rain evaporation rate
- Correctly calls isotope framework
- Unit tests pass

**Function Signature**:
```fortran
subroutine wtrc_p3_inter_phase2(state, ptend, pbuf, precr, preci, top_lev, &
                                dtime, qr_evap_rate, T, p, qv, liqcldf, icecldf)
```

---

#### Issue #2.3: Port stewart_isoevap() Function
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-2`, `fractionation`, `stewart-model`  
**Assignees**: TBD  
**Dependencies**: #2.2  
**Estimated Effort**: 3-4 days  

**Description**:
Port the Stewart (1975) rain evaporation isotope model from EAMv2-wiso. This is the core physics for rain evaporation fractionation.

**Source**: `water_tracers.F90:3678` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `stewart_isoevap()` function
- [ ] Update dependencies
- [ ] Test with idealized inputs
- [ ] Verify fractionation factors
- [ ] Test Marshall-Palmer raindrop size calculation
- [ ] Verify equilibration calculation
- [ ] Test full Stewart equation
- [ ] Compare results to EAMv2-wiso
- [ ] Document function
- [ ] Create unit tests

**Acceptance Criteria**:
- Function compiles without errors
- Produces correct isotope fractionation
- Results match EAMv2-wiso for test cases
- Unit tests pass

**Key Physics**:
- Stewart equation: R_stw = ((R_r0 - γ·R_v0)·f^β) + γ·R_v0
- Equilibrium fractionation: γ = (α_e·h)/(1 - α_e·α_k·(1-h))
- Kinetic fractionation: β = (1 - α_e·α_k·(1-h))/(α_e·α_k·(1-h))
- Partial equilibration applied based on droplet size and fall time

---

#### Issue #2.4: Port wtrc_equil_time() Function
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-2`, `equilibration`, `utilities`  
**Assignees**: TBD  
**Dependencies**: #2.3  
**Estimated Effort**: 2-3 days  

**Description**:
Port function that calculates the fraction of rain that experiences equilibration during fall. Critical for realistic below-cloud fractionation.

**Source**: `water_tracers.F90:6580` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_equil_time()` function
- [ ] Test equilibration time calculation
- [ ] Verify droplet size dependence
- [ ] Test fall distance effects
- [ ] Verify diffusion time scale
- [ ] Create unit tests
- [ ] Document function

**Acceptance Criteria**:
- Function compiles and runs
- Equilibration fraction in range [0, 1]
- Results physically reasonable
- Unit tests pass

**Inputs**:
- Isotope species
- Temperature
- Pressure
- Raindrop radius
- Layer thickness (dz)
- Fractionation factor
- Diffusivity

**Output**:
- `fequil` - Fraction experiencing equilibration (0-1)

---

#### Issue #2.5: Port wtrc_liqvap_equil() Function
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-2`, `equilibration`, `solver`  
**Assignees**: TBD  
**Dependencies**: #2.4  
**Estimated Effort**: 2-3 days  

**Description**:
Port iterative liquid-vapor equilibration solver. Exchanges isotopes between rain and vapor phases.

**Source**: `water_tracers.F90:6381` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_liqvap_equil()` function
- [ ] Test iterative solver convergence
- [ ] Verify mass conservation
- [ ] Test with various fractionation factors
- [ ] Test partial equilibration
- [ ] Create unit tests
- [ ] Document function

**Acceptance Criteria**:
- Function compiles and runs
- Solver converges in reasonable iterations
- Mass is conserved
- Unit tests pass

**Algorithm**:
1. Calculate equilibrium distribution
2. Apply partial equilibration fraction
3. Exchange isotope mass between phases
4. Verify mass conservation

---

#### Issue #2.6: Port/Simplify wtrc_apply_rates() for Phase 2
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-2`, `core`, `apply-rates`  
**Assignees**: TBD  
**Dependencies**: #2.5  
**Estimated Effort**: 5-7 days  

**Description**:
Port and simplify the core isotope tendency calculator from EAMv2-wiso. Phase 2 version handles ONLY rain evaporation. Full version will be extended in later phases.

**Source**: `water_tracers.F90:2068` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_apply_rates()` function skeleton
- [ ] Simplify to handle only rain evaporation
- [ ] Implement iterative solver framework
- [ ] Call `stewart_isoevap()` for rain-vapor exchange
- [ ] Apply fractionation to tracers
- [ ] Calculate tendencies (ptend)
- [ ] Update precipitation isotope ratios
- [ ] Test with single-column case
- [ ] Verify mass conservation
- [ ] Document simplifications

**Acceptance Criteria**:
- Function compiles and runs
- Rain evaporation fractionation works correctly
- Mass is conserved
- Tendencies calculated correctly
- Single-column tests pass

**Simplified for Phase 2**:
- Only handles: iwtstrain → iwtvap (rain evaporation)
- Does NOT handle (deferred):
  - Condensation (Phase 3)
  - Ice processes (Phase 4)
  - Freezing/melting (Phase 5)
  - Collection (Phase 6)

---

#### Issue #2.7: Integrate Phase 2 Interface with P3
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-2`, `integration`, `p3-interface`  
**Assignees**: TBD  
**Dependencies**: #2.6  
**Estimated Effort**: 3-4 days  

**Description**:
Integrate Phase 2 isotope interface into P3 microphysics interface. Ensure proper calling sequence and data flow.

**File to Modify**: `components/eam/src/physics/p3/eam/micro_p3_interface.F90`

**Tasks**:
- [ ] Add call to `wtrc_p3_inter_phase2()` in `micro_p3_tend()`
- [ ] Pass rain evaporation rate from P3
- [ ] Pass thermodynamic state (T, p, qv)
- [ ] Pass cloud fractions
- [ ] Get isotope tendencies back
- [ ] Add diagnostic output
- [ ] Test integration
- [ ] Verify no errors when isotopes disabled

**Acceptance Criteria**:
- Interface called correctly
- Data flow works properly
- Isotope tendencies applied to tracers
- Code runs with and without isotopes
- No performance degradation when isotopes disabled

**Calling Sequence**:
```fortran
if (wisotope) then
  call wtrc_p3_inter_phase2(state, ptend, pbuf, precr, preci, &
                            top_lev, dtime, qr_evap_rate, &
                            T, p, qv, liqcldf, icecldf)
endif
```

---

#### Issue #2.8: Phase 2 Testing and Validation
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-2`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #2.7  
**Estimated Effort**: 5-7 days  

**Description**:
Comprehensive testing of Phase 2 implementation. Validate rain evaporation isotope fractionation against theory and observations.

**Tasks**:
- [ ] Create single-column test case (rain falling through dry layer)
- [ ] Run test with different humidity levels
- [ ] Verify δD and δ18O evolution
- [ ] Compare to Stewart model predictions
- [ ] Compare to EAMv2-wiso MG2 results
- [ ] Test mass conservation
- [ ] Check for numerical issues
- [ ] Performance testing
- [ ] Create validation report
- [ ] Document test cases

**Acceptance Criteria**:
- Rain evaporation reduces δD (more negative)
- Magnitude consistent with Stewart model
- Mass conserved to < 0.01%
- Results comparable to EAMv2-wiso
- No numerical instabilities
- Performance acceptable

**Test Cases**:
1. **Dry layer evaporation**: Rain falling through RH=50% layer
2. **Humidity sensitivity**: Test RH=30%, 50%, 70%, 90%
3. **Temperature sensitivity**: Test T=273K, 283K, 293K
4. **Rain rate sensitivity**: Light vs heavy rain
5. **Isotope closure**: Verify δD-δ18O relationship

**Expected Behavior**:
- Higher evaporation → more negative δD
- Lower RH → stronger fractionation
- Warmer T → weaker fractionation
- Larger drops → less equilibration

---

### Phase 3: Vapor-Liquid Phase Changes

#### Issue #3.1: Modify cloud_water_autoconversion()
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-3`, `p3-modification`, `autoconversion`  
**Assignees**: TBD  
**Dependencies**: #2.8  
**Estimated Effort**: 2-3 days  

**Description**:
Modify P3 autoconversion function to output rate for isotope tracking. Autoconversion itself doesn't fractionate (liquid → liquid) but must be tracked for mass conservation.

**File to Modify**: `micro_p3.F90` (~line 1500)

**Tasks**:
- [ ] Add autoconversion rate output
- [ ] Verify rate units (kg/kg/s)
- [ ] Test rate calculation
- [ ] Add to diagnostic output
- [ ] Document changes

**Acceptance Criteria**:
- Rate correctly calculated and output
- Mass balance maintained
- No change to P3 behavior when isotopes disabled

---

#### Issue #3.2: Modify cloud_rain_accretion()
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-3`, `p3-modification`, `accretion`  
**Assignees**: TBD  
**Dependencies**: #3.1  
**Estimated Effort**: 2-3 days  

**Description**:
Modify P3 accretion function to output rate for isotope tracking. No fractionation but mass tracking required.

**File to Modify**: `micro_p3.F90` (~line 1600)

**Tasks**:
- [ ] Add accretion rate output
- [ ] Verify rate units
- [ ] Test rate calculation
- [ ] Add to diagnostic output
- [ ] Document changes

**Acceptance Criteria**:
- Rate correctly output
- Mass conserved
- No change when isotopes disabled

---

#### Issue #3.3: Extract Condensation Rate from P3
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-3`, `p3-modification`, `condensation`  
**Assignees**: TBD  
**Dependencies**: #3.2  
**Estimated Effort**: 3-4 days  

**Description**:
Extract condensation rate from P3's saturation adjustment. This is implicit in P3 and may require calculation from before/after values.

**File to Modify**: `micro_p3.F90` (main scheme)

**Tasks**:
- [ ] Identify where supersaturation is removed
- [ ] Calculate condensation rate (change in qv, qc)
- [ ] Output rate in compatible format
- [ ] Handle cloud evaporation (negative condensation)
- [ ] Test rate extraction
- [ ] Document approach

**Acceptance Criteria**:
- Condensation rate correctly calculated
- Both condensation and evaporation handled
- Mass conserved

**Note**: P3 doesn't explicitly calculate condensation rate like MG2, so this may require:
- Option 1: Calculate from Δqv and Δqc
- Option 2: Modify P3 to output explicit rate
- Requires coordination with P3 developers

---

#### Issue #3.4: Port wtrc_vap_distil() for Condensation
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-3`, `fractionation`, `rayleigh`  
**Assignees**: TBD  
**Dependencies**: #3.3  
**Estimated Effort**: 2-3 days  

**Description**:
Port Rayleigh distillation routine for condensation fractionation. This handles fractionation when vapor condenses to liquid.

**Source**: `water_tracers.F90:2776` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_vap_distil()` function
- [ ] Test Rayleigh fractionation calculation
- [ ] Verify with analytical solutions
- [ ] Test various α values
- [ ] Create unit tests
- [ ] Document function

**Acceptance Criteria**:
- Function compiles and runs
- Rayleigh fractionation correct
- Results match analytical solutions
- Unit tests pass

**Physics**:
Rayleigh distillation: R_liquid / R_vapor = α^(1/(α-1)) · (f^(α-1))
where f = fraction of vapor remaining

---

#### Issue #3.5: Extend wtrc_apply_rates() for Condensation
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-3`, `core`, `apply-rates`  
**Assignees**: TBD  
**Dependencies**: #3.4  
**Estimated Effort**: 4-5 days  

**Description**:
Extend `wtrc_apply_rates()` to handle vapor-liquid phase changes (condensation and evaporation) in addition to rain evaporation.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Add condensation rate handling (iwtvap → iwtliq)
- [ ] Add cloud evaporation rate handling (iwtliq → iwtvap)
- [ ] Implement Rayleigh distillation for condensation
- [ ] Implement equilibrium fractionation for evaporation
- [ ] Handle autoconversion (iwtliq → iwtstrain)
- [ ] Handle accretion (iwtliq → iwtstrain)
- [ ] Test all processes together
- [ ] Verify mass conservation
- [ ] Document extensions

**Acceptance Criteria**:
- All vapor-liquid processes handled
- Fractionation applied correctly
- Mass conserved
- Tests pass

**Processes Added**:
1. Condensation: iwtvap → iwtliq (Rayleigh)
2. Cloud evaporation: iwtliq → iwtvap (equilibrium)
3. Autoconversion: iwtliq → iwtstrain (conservative)
4. Accretion: iwtliq → iwtstrain (conservative)

---

#### Issue #3.6: Phase 3 Testing and Validation
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-3`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #3.5  
**Estimated Effort**: 4-5 days  

**Description**:
Test and validate vapor-liquid phase change isotope fractionation.

**Tasks**:
- [ ] Create cloud formation test case
- [ ] Test condensation fractionation
- [ ] Test cloud evaporation fractionation
- [ ] Verify cloud water δD values
- [ ] Compare to EAMv2-wiso
- [ ] Test mass conservation
- [ ] Create validation report

**Acceptance Criteria**:
- Condensation depletes vapor (more negative δD)
- Cloud water enriched relative to vapor
- Fractionation factors correct
- Mass conserved
- Results comparable to EAMv2-wiso

**Test Cases**:
1. Adiabatic ascent with condensation
2. Cloud evaporation
3. Combined condensation and rain formation
4. Isotope closure verification

---

### Phase 4: Vapor-Ice Phase Changes

#### Issue #4.1: Modify ice_nucleation()
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-4`, `p3-modification`, `nucleation`  
**Assignees**: TBD  
**Dependencies**: #3.6  
**Estimated Effort**: 2-3 days  

**Description**:
Modify P3 ice nucleation function to output nucleation rate and supersaturation for isotope fractionation.

**File to Modify**: `micro_p3.F90` (~line 1200)

**Tasks**:
- [ ] Add ice nucleation rate output
- [ ] Output supersaturation (Si)
- [ ] Output temperature
- [ ] Test rate calculation
- [ ] Document changes

**Acceptance Criteria**:
- Nucleation rate correctly output
- Supersaturation provided for kinetic fractionation
- Mass conserved

---

#### Issue #4.2: Modify ice_deposition_sublimation()
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-4`, `p3-modification`, `deposition`  
**Assignees**: TBD  
**Dependencies**: #4.1  
**Estimated Effort**: 3-4 days  

**Description**:
Modify P3 ice deposition/sublimation function to output bidirectional rates separately. Critical for correct fractionation sign.

**File to Modify**: `micro_p3.F90` (~line 2200)

**Tasks**:
- [ ] Separate deposition and sublimation rates
- [ ] Output each rate separately
- [ ] Provide supersaturation for both directions
- [ ] Test rate extraction
- [ ] Document changes

**Acceptance Criteria**:
- Deposition and sublimation rates separated
- Correct sign for each process
- Supersaturation available
- Mass conserved

**Critical**: Must distinguish:
- Deposition (positive): iwtvap → iwtice (with kinetic fractionation)
- Sublimation (negative): iwtice → iwtvap (equilibrium only)

---

#### Issue #4.3: Modify ice_supersat_conservation()
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-4`, `p3-modification`, `supersaturation`  
**Assignees**: TBD  
**Dependencies**: #4.2  
**Estimated Effort**: 2-3 days  

**Description**:
Modify function that removes ice supersaturation to output deposition rate.

**File to Modify**: `micro_p3.F90` (~line 2300)

**Tasks**:
- [ ] Add supersaturation deposition rate output
- [ ] Provide Si before and after
- [ ] Test rate calculation
- [ ] Document changes

**Acceptance Criteria**:
- Rate correctly output
- Supersaturation adjustment tracked
- Mass conserved

---

#### Issue #4.4: Implement Ice-Vapor Fractionation
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-4`, `fractionation`, `ice`  
**Assignees**: TBD  
**Dependencies**: #4.3  
**Estimated Effort**: 4-5 days  

**Description**:
Implement ice-vapor fractionation physics in `wtrc_apply_rates()`. This is more complex than liquid-vapor due to supersaturation dependence.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Add ice nucleation handling (iwtvap → iwtice)
- [ ] Add deposition handling with kinetic fractionation
- [ ] Add sublimation handling (equilibrium only)
- [ ] Implement supersaturation-dependent kinetics
- [ ] Use `wiso_alpi()` for equilibrium fractionation
- [ ] Use `wiso_akci()` for kinetic fractionation
- [ ] Test ice processes
- [ ] Document implementation

**Acceptance Criteria**:
- Ice-vapor fractionation works correctly
- Supersaturation effects included
- Deposition and sublimation distinguished
- Tests pass

**Key Physics**:
- Equilibrium: `α_eq = wiso_alpi(T, ispec)`
- Kinetic: `α_k = wiso_akci(Si, ispec)`
- Total deposition: `α_total = α_eq × α_k`
- Sublimation: `α_total = α_eq` (no kinetic)

---

#### Issue #4.5: Extend wtrc_apply_rates() for Ice Processes
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-4`, `core`, `apply-rates`  
**Assignees**: TBD  
**Dependencies**: #4.4  
**Estimated Effort**: 5-6 days  

**Description**:
Fully integrate ice processes into `wtrc_apply_rates()`. Now handles vapor-liquid AND vapor-ice.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Integrate ice nucleation fractionation
- [ ] Integrate deposition fractionation (with kinetics)
- [ ] Integrate sublimation fractionation (equilibrium)
- [ ] Handle supersaturation adjustment
- [ ] Test all processes together
- [ ] Test vapor-liquid-ice transitions
- [ ] Verify mass conservation across all phases
- [ ] Document complete implementation

**Acceptance Criteria**:
- All vapor-ice processes handled
- Coexistence with vapor-liquid processes works
- Mass conserved across all phases
- Tests pass

**Processes Handled** (Phase 4):
- Vapor ↔ liquid (from Phase 3)
- Vapor ↔ ice (new)
- Liquid ↔ liquid (from Phase 3)

**Processes Deferred**:
- Liquid ↔ ice (Phase 5)
- Collection (Phase 6)

---

#### Issue #4.6: Phase 4 Testing and Validation
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-4`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #4.5  
**Estimated Effort**: 5-6 days  

**Description**:
Test and validate vapor-ice phase change isotope fractionation.

**Tasks**:
- [ ] Create cirrus cloud test case (pure ice)
- [ ] Create mixed-phase cloud test case
- [ ] Test ice nucleation fractionation
- [ ] Test depositional growth fractionation
- [ ] Test sublimation fractionation
- [ ] Verify supersaturation effects
- [ ] Test ice δD and δ18O values
- [ ] Compare to observations
- [ ] Compare to EAMv2-wiso
- [ ] Create validation report

**Acceptance Criteria**:
- Ice depleted relative to vapor
- Deposition shows kinetic effects
- Supersaturation enhances depletion
- Results physically reasonable
- Mass conserved

**Test Cases**:
1. Cirrus cloud formation (T < -40°C)
2. Mixed-phase cloud (0°C > T > -40°C)
3. Supersaturation sensitivity
4. Comparison to ice core data

---

### Phase 5: Freezing and Melting

#### Issue #5.1: Modify All Freezing Functions
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-5`, `p3-modification`, `freezing`  
**Assignees**: TBD  
**Dependencies**: #4.6  
**Estimated Effort**: 3-4 days  

**Description**:
Modify all P3 freezing functions to output rates: immersion freezing (cloud and rain) and homogeneous freezing.

**Files to Modify**: `micro_p3.F90`
- `cldliq_immersion_freezing()` (~line 1800)
- `rain_immersion_freezing()` (~line 1900)
- `homogeneous_freezing()` (~line 1750)

**Tasks**:
- [ ] Modify `cldliq_immersion_freezing()` for rate output
- [ ] Modify `rain_immersion_freezing()` for rate output
- [ ] Modify `homogeneous_freezing()` for rate output
- [ ] Provide temperature for fractionation
- [ ] Test rate calculations
- [ ] Document changes

**Acceptance Criteria**:
- All freezing rates correctly output
- Temperatures provided
- Mass conserved

---

#### Issue #5.2: Modify All Melting Functions
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-5`, `p3-modification`, `melting`  
**Assignees**: TBD  
**Dependencies**: #5.1  
**Estimated Effort**: 2-3 days  

**Description**:
Modify P3 melting functions to output rates.

**Files to Modify**: `micro_p3.F90`
- `ice_melting()` (~line 2500)
- `ice_complete_melting()` (~line 2550)

**Tasks**:
- [ ] Modify `ice_melting()` for rate output
- [ ] Modify `ice_complete_melting()` for rate output
- [ ] Test rate calculations
- [ ] Document changes

**Acceptance Criteria**:
- Melting rates correctly output
- Mass conserved

---

#### Issue #5.3: Implement Freezing Fractionation
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-5`, `fractionation`, `freezing`  
**Assignees**: TBD  
**Dependencies**: #5.2  
**Estimated Effort**: 3-4 days  

**Description**:
Implement isotope fractionation during freezing processes. Generally small effects except at warm temperatures.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Add immersion freezing fractionation
- [ ] Add rain freezing fractionation
- [ ] Add homogeneous freezing (no fractionation)
- [ ] Implement temperature-dependent kinetics
- [ ] Test freezing fractionation
- [ ] Document implementation

**Acceptance Criteria**:
- Freezing fractionation implemented
- Temperature dependence correct
- Homogeneous freezing conservative
- Tests pass

**Fractionation Approach**:
- Immersion freezing: Small kinetic effect (α_k ≈ 1.001-1.005)
- Rain freezing: Minimal (fast process)
- Homogeneous freezing: None (α = 1.0, instantaneous)

---

#### Issue #5.4: Implement Melting (Conservative Transfer)
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-5`, `melting`, `conservative`  
**Assignees**: TBD  
**Dependencies**: #5.3  
**Estimated Effort**: 2-3 days  

**Description**:
Implement melting as conservative mass transfer (no fractionation at 0°C).

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Add ice melting to liquid
- [ ] Add complete melting
- [ ] Add snow melting to rain
- [ ] Verify conservative transfer
- [ ] Test melting
- [ ] Document implementation

**Acceptance Criteria**:
- Melting transfers isotopes conservatively
- Mass and isotope ratios conserved
- Tests pass

---

#### Issue #5.5: Extend wtrc_apply_rates() for Freeze/Melt
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-5`, `core`, `apply-rates`  
**Assignees**: TBD  
**Dependencies**: #5.4  
**Estimated Effort**: 4-5 days  

**Description**:
Fully integrate freezing and melting into `wtrc_apply_rates()`.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Integrate all freezing processes
- [ ] Integrate all melting processes
- [ ] Test all three phase types together
- [ ] Verify mass conservation
- [ ] Test mixed-phase cloud evolution
- [ ] Document implementation

**Acceptance Criteria**:
- All phase transitions handled
- Vapor-liquid, vapor-ice, liquid-ice all working
- Mass conserved across all transitions
- Tests pass

**Complete Process List** (Phases 2-5):
- Vapor ↔ liquid (Phase 3)
- Vapor ↔ ice (Phase 4)  
- Liquid ↔ ice (Phase 5)
- Liquid → liquid (Phase 3)

---

#### Issue #5.6: Phase 5 Testing and Validation
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-5`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #5.5  
**Estimated Effort**: 4-5 days  

**Description**:
Test and validate freezing and melting isotope effects.

**Tasks**:
- [ ] Create freezing rain test case
- [ ] Create melting snow test case
- [ ] Test mixed-phase cloud evolution
- [ ] Verify phase transition isotope effects
- [ ] Test mass conservation through phase changes
- [ ] Compare to EAMv2-wiso
- [ ] Create validation report

**Acceptance Criteria**:
- Freezing/melting isotope effects correct
- Mass conserved through transitions
- Mixed-phase clouds behave physically
- Results comparable to EAMv2-wiso

**Test Cases**:
1. Supercooled liquid freezing
2. Ice melting at 0°C
3. Mixed-phase cloud lifecycle
4. Freezing rain event

---

### Phase 6: Collection and Riming

#### Issue #6.1: Modify Ice-Liquid Collection Functions
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-6`, `p3-modification`, `collection`  
**Assignees**: TBD  
**Dependencies**: #5.6  
**Estimated Effort**: 3-4 days  

**Description**:
Modify P3 collection functions to output rates for riming and wet growth.

**Files to Modify**: `micro_p3.F90`
- `ice_cldliq_collection()` (~line 2000)
- `ice_cldliq_wet_growth()` (~line 2100)
- `ice_rain_collection()` (~line 2150)

**Tasks**:
- [ ] Modify `ice_cldliq_collection()` for riming rate
- [ ] Modify `ice_cldliq_wet_growth()` for wet growth/shedding
- [ ] Modify `ice_rain_collection()` for rate
- [ ] Provide rime mass fraction
- [ ] Test rate calculations
- [ ] Document changes

**Acceptance Criteria**:
- Collection rates correctly output
- Wet growth split into freezing and shedding
- Rime mass tracked
- Mass conserved

---

#### Issue #6.2: Implement Collection Fractionation
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-6`, `fractionation`, `riming`  
**Assignees**: TBD  
**Dependencies**: #6.1  
**Estimated Effort**: 4-5 days  

**Description**:
Implement isotope fractionation during collection processes (riming, wet growth).

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Add riming fractionation (kinetic effects)
- [ ] Add wet growth fractionation (partial equilibration)
- [ ] Add rain collection by ice
- [ ] Handle shedding in wet growth
- [ ] Test collection fractionation
- [ ] Document implementation

**Acceptance Criteria**:
- Collection fractionation implemented
- Wet growth handled correctly
- Rime mass isotope composition tracked
- Tests pass

**Fractionation Approach**:
- Dry riming: Small kinetic effect during impact/freezing
- Wet growth: Partial equilibration, fraction depends on T
- Shedding: Conservative (liquid shed as rain)

---

#### Issue #6.3: Handle Bergeron Process
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-6`, `bergeron`, `vapor-transfer`  
**Assignees**: TBD  
**Dependencies**: #6.2  
**Estimated Effort**: 3-4 days  

**Description**:
Handle Bergeron process where vapor transfers from liquid to ice via vapor phase (two-step fractionation).

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Identify Bergeron process in P3 (implicit)
- [ ] Implement two-step fractionation (liq→vap→ice)
- [ ] Test Bergeron fractionation
- [ ] Document approach

**Acceptance Criteria**:
- Bergeron process isotope effects correct
- Two-step fractionation works
- Tests pass

**Note**: In P3 this is implicit in supersaturation adjustment, not explicit like MG2 `bergo`/`bergso` rates.

---

#### Issue #6.4: Extend wtrc_apply_rates() for Collection
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-6`, `core`, `apply-rates`  
**Assignees**: TBD  
**Dependencies**: #6.3  
**Estimated Effort**: 4-5 days  

**Description**:
Integrate collection processes into `wtrc_apply_rates()`.

**File to Modify**: `water_tracers.F90`

**Tasks**:
- [ ] Integrate riming
- [ ] Integrate wet growth
- [ ] Integrate rain collection
- [ ] Integrate Bergeron process
- [ ] Test all collection processes
- [ ] Verify mass conservation
- [ ] Document implementation

**Acceptance Criteria**:
- All collection processes handled
- Coexistence with other processes works
- Mass conserved
- Tests pass

---

#### Issue #6.5: Phase 6 Testing and Validation
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-6`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #6.4  
**Estimated Effort**: 3-4 days  

**Description**:
Test and validate collection process isotope effects.

**Tasks**:
- [ ] Create riming test case
- [ ] Test wet growth shedding
- [ ] Test Bergeron process
- [ ] Verify rime isotope composition
- [ ] Compare to observations
- [ ] Create validation report

**Acceptance Criteria**:
- Collection isotope effects physically reasonable
- Rime mass composition tracked
- Mass conserved
- No spurious signals

**Test Cases**:
1. Ice riming cloud droplets
2. Wet growth with shedding
3. Bergeron process in mixed-phase cloud

---

### Phase 7: Sedimentation with Fractionation

#### Issue #7.1: Modify rain_sedimentation()
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-7`, `p3-modification`, `sedimentation`  
**Assignees**: TBD  
**Dependencies**: #6.5  
**Estimated Effort**: 3-4 days  

**Description**:
Modify rain sedimentation to provide layer-by-layer evaporation rates and fall velocities for isotope fractionation.

**File to Modify**: `micro_p3.F90` (~line 3200)

**Tasks**:
- [ ] Extract rain flux per layer
- [ ] Extract evaporation rate per layer
- [ ] Provide fall velocity
- [ ] Provide layer thickness
- [ ] Test rate extraction
- [ ] Document changes

**Acceptance Criteria**:
- Sedimentation and evaporation rates separated per layer
- Data format compatible with isotope interface
- Mass conserved

---

#### Issue #7.2: Modify ice_sedimentation()
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-7`, `p3-modification`, `sedimentation`  
**Assignees**: TBD  
**Dependencies**: #7.1  
**Estimated Effort**: 3-4 days  

**Description**:
Modify ice sedimentation to provide layer-by-layer sublimation rates.

**File to Modify**: `micro_p3.F90` (~line 3400)

**Tasks**:
- [ ] Extract ice flux per layer
- [ ] Extract sublimation rate per layer
- [ ] Provide fall velocity
- [ ] Test rate extraction
- [ ] Document changes

**Acceptance Criteria**:
- Sedimentation and sublimation separated per layer
- Data format compatible
- Mass conserved

---

#### Issue #7.3: Port wtrc_sediment() Function
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-7`, `sedimentation`, `stewart-model`  
**Assignees**: TBD  
**Dependencies**: #7.2  
**Estimated Effort**: 5-7 days  

**Description**:
Port complete sedimentation isotope handler from EAMv2-wiso. This handles layer-by-layer fractionation during precipitation fall.

**Source**: `water_tracers.F90:3315` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_sediment()` function
- [ ] Adapt for P3 sedimentation structure
- [ ] Implement CFL sub-stepping
- [ ] Implement layer-by-layer Stewart model
- [ ] Handle rain partial equilibration
- [ ] Handle ice sublimation
- [ ] Test sedimentation isotopes
- [ ] Document function

**Acceptance Criteria**:
- Function compiles and runs
- CFL condition stable
- Layer-by-layer fractionation works
- Mass conserved
- Tests pass

**Key Features**:
- CFL-limited sub-stepping
- Stewart model applied per layer
- Raindrop size evolution during fall
- Clear vs cloudy sky separation
- Partial equilibration calculation

---

#### Issue #7.4: Integrate Sedimentation with P3
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-7`, `integration`, `sedimentation`  
**Assignees**: TBD  
**Dependencies**: #7.3  
**Estimated Effort**: 4-5 days  

**Description**:
Integrate sedimentation isotope handler with P3 interface.

**File to Modify**: `water_tracers_p3.F90`

**Tasks**:
- [ ] Call `wtrc_sediment()` from P3 interface
- [ ] Pass sedimentation fluxes
- [ ] Pass evaporation/sublimation rates
- [ ] Handle below-cloud fractionation
- [ ] Update surface precipitation isotopes
- [ ] Test integration
- [ ] Document implementation

**Acceptance Criteria**:
- Sedimentation isotopes work correctly
- Below-cloud fractionation applied
- Surface precipitation isotopes updated
- Mass conserved

---

#### Issue #7.5: Phase 7 Testing and Validation
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-7`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #7.4  
**Estimated Effort**: 5-7 days  

**Description**:
Test and validate sedimentation isotope fractionation. This is critical for surface precipitation isotope composition.

**Tasks**:
- [ ] Create precipitation through dry layer test
- [ ] Test humidity sensitivity
- [ ] Test vertical profiles
- [ ] Verify δD evolution during fall
- [ ] Compare to observations (GNIP)
- [ ] Compare to EAMv2-wiso
- [ ] Test mass conservation
- [ ] Create validation report

**Acceptance Criteria**:
- Below-cloud evaporation enriches precipitation
- Magnitude consistent with observations
- Vertical profiles physically reasonable
- Mass conserved
- Results comparable to EAMv2-wiso

**Test Cases**:
1. Rain falling through RH gradient
2. Snow sublimation
3. Comparison to GNIP station data
4. Vertical isotope profile validation

---

### Phase 8: Convection Interface

#### Issue #8.1: Port wtrc_q1q2_pjr() for ZM Deep Convection
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-8`, `convection`, `deep`  
**Assignees**: TBD  
**Dependencies**: #7.5  
**Estimated Effort**: 5-7 days  

**Description**:
Port ZM deep convection isotope interface from EAMv2-wiso.

**Source**: `water_tracers.F90:4960` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_q1q2_pjr()` function
- [ ] Update for EAMv3 ZM interface
- [ ] Handle updraft condensation fractionation
- [ ] Handle downdraft evaporation fractionation
- [ ] Handle detrainment to stratiform
- [ ] Track convective precipitation isotopes
- [ ] Test convection isotopes
- [ ] Document function

**Acceptance Criteria**:
- Function compiles and runs with EAMv3 ZM
- Convective fractionation works
- Detrainment handled correctly
- Tests pass

**Processes Handled**:
- Updraft: Rayleigh condensation
- Downdraft: Partial re-evaporation
- Convective precip: Tracked separately from stratiform
- Detrainment: Cloud water and ice to environment

---

#### Issue #8.2: Port wtrc_shallow() for Shallow Convection
**Status**: To Do  
**Priority**: Medium  
**Labels**: `phase-8`, `convection`, `shallow`  
**Assignees**: TBD  
**Dependencies**: #8.1  
**Estimated Effort**: 3-4 days  

**Description**:
Port shallow convection isotope interface from EAMv2-wiso.

**Source**: `water_tracers.F90:6234` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Copy `wtrc_shallow()` function
- [ ] Update for EAMv3 shallow convection (UW or Park)
- [ ] Handle shallow condensation
- [ ] Handle shallow detrainment
- [ ] Test shallow convection isotopes
- [ ] Document function

**Acceptance Criteria**:
- Function works with EAMv3 shallow convection
- Isotope effects physically reasonable
- Tests pass

---

#### Issue #8.3: Phase 8 Testing and Validation
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-8`, `testing`, `validation`  
**Assignees**: TBD  
**Dependencies**: #8.2  
**Estimated Effort**: 4-5 days  

**Description**:
Test and validate convective precipitation isotope effects.

**Tasks**:
- [ ] Create tropical deep convection test
- [ ] Test shallow convection case
- [ ] Compare convective vs stratiform isotope signatures
- [ ] Verify total precipitation (convective + stratiform)
- [ ] Compare to observations
- [ ] Compare to EAMv2-wiso
- [ ] Create validation report

**Acceptance Criteria**:
- Convective precipitation more depleted than stratiform
- Convective vs stratiform separation works
- Total water budget closed
- Results physically reasonable

**Test Cases**:
1. Tropical deep convection
2. Trade cumulus (shallow)
3. Convective vs stratiform partitioning
4. Tropical isotope observations comparison

---

### Phase 9: Conservation and Diagnostics

#### Issue #9.1: Port Conservation Functions
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-9`, `conservation`, `validation`  
**Assignees**: TBD  
**Dependencies**: #8.3  
**Estimated Effort**: 3-4 days  

**Description**:
Port mass conservation checking and fixing functions from EAMv2-wiso.

**Source**: `water_tracers.F90` (multiple functions) in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Port `wtrc_check_h2o()` - Mass conservation check
- [ ] Port `wtrc_adjust_h2o()` - Fix small negatives
- [ ] Port `wtrc_mass_fixer()` - Reset tracers if needed
- [ ] Extend P3 conservation functions for isotopes
- [ ] Test conservation checks
- [ ] Document functions

**Acceptance Criteria**:
- Conservation checks work correctly
- Small numerical errors fixed
- Mass conservation enforced
- Tests pass

---

#### Issue #9.2: Port Diagnostic Functions and Output
**Status**: To Do  
**Priority**: High  
**Labels**: `phase-9`, `diagnostics`, `output`  
**Assignees**: TBD  
**Dependencies**: #9.1  
**Estimated Effort**: 4-5 days  

**Description**:
Port diagnostic output functions and add isotope fields to history files.

**Source**: `water_tracers.F90` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Port `wtrc_output_precip()` - Precipitation isotope output
- [ ] Port `wtrc_setup_diag()` - Diagnostic setup
- [ ] Add δD, δ18O, δ17O calculation routines
- [ ] Add history fields for all phases
- [ ] Add precipitation isotope output
- [ ] Add budget diagnostics
- [ ] Test output
- [ ] Document output fields

**Acceptance Criteria**:
- All isotope fields in history files
- δD, δ18O, δ17O calculated correctly
- Precipitation isotopes output
- Output compatible with EAMv2-wiso format

**Output Fields** (minimum):
- `HDO_vapor`, `H218O_vapor` - Atmospheric vapor
- `HDO_cld`, `H218O_cld` - Cloud liquid
- `HDO_ice`, `H218O_ice` - Cloud ice
- `PRECIP_HDO`, `PRECIP_H218O` - Total precipitation
- `PRECRL_HDO`, `PRECRL_H218O` - Stratiform rain
- `PRECRS_HDO`, `PRECRS_H218O` - Stratiform snow
- `PRECRC_HDO`, `PRECRC_H218O` - Convective rain
- `PRECSC_HDO`, `PRECSC_H218O` - Convective snow
- δ notation fields: `dD_vapor`, `d18O_vapor`, etc.

---

#### Issue #9.3: Phase 9 Testing and Global Validation
**Status**: To Do  
**Priority**: Critical  
**Labels**: `phase-9`, `testing`, `validation`, `global`  
**Assignees**: TBD  
**Dependencies**: #9.2  
**Estimated Effort**: 7-10 days  

**Description**:
Comprehensive global testing and validation of complete isotope implementation.

**Tasks**:
- [ ] Run global simulation (multi-month)
- [ ] Check global mass conservation
- [ ] Validate against GNIP precipitation observations
- [ ] Validate against atmospheric vapor observations
- [ ] Compare to EAMv2-wiso global results
- [ ] Verify δD-δ18O relationship (meteoric water line)
- [ ] Check seasonal cycle
- [ ] Check spatial patterns
- [ ] Performance testing
- [ ] Create comprehensive validation report

**Acceptance Criteria**:
- Global mass conserved to < 0.1% error
- Spatial patterns match observations
- Seasonal cycle correct
- Results comparable to EAMv2-wiso
- Performance acceptable (< 20% slowdown)
- No numerical instabilities

**Validation Datasets**:
1. GNIP precipitation isotopes (global network)
2. TES/AIRS satellite vapor isotopes
3. Aircraft campaigns (vertical profiles)
4. Ice core records (for snow)
5. EAMv2-wiso benchmark simulations

**Metrics**:
- Global mean δD and δ18O
- Latitude gradients
- Altitude effects
- Temperature correlation
- Amount effect
- δD-δ18O slope (should be ~8 for precipitation)

---

### Phase 10: Additional Physics (Optional/Future)

#### Issue #10.1: Port Methane Oxidation Source
**Status**: To Do  
**Priority**: Low  
**Labels**: `phase-10`, `chemistry`, `optional`  
**Assignees**: TBD  
**Dependencies**: #9.3  
**Estimated Effort**: 2-3 days  

**Description**:
Port methane oxidation water vapor source for stratospheric isotopes.

**Source**: `water_tracers.F90:6751` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Port `wtrc_chem_ch4ox_tend()`
- [ ] Integrate with chemistry module
- [ ] Test stratospheric HDO source
- [ ] Document function

**Acceptance Criteria**:
- Methane oxidation adds vapor with correct isotope ratio
- Stratospheric vapor isotopes improved
- Tests pass

---

#### Issue #10.2: Port Radioactive Decay (HTO)
**Status**: To Do  
**Priority**: Low  
**Labels**: `phase-10`, `tritium`, `optional`  
**Assignees**: TBD  
**Dependencies**: #10.1  
**Estimated Effort**: 1-2 days  

**Description**:
Port tritium radioactive decay for HTO tracer.

**Source**: `water_tracers.F90:6832` in EAMv2-wiso  
**Destination**: `water_tracers.F90` in EAMv3-wiso

**Tasks**:
- [ ] Port `wtrc_rad_decay()`
- [ ] Implement HTO decay (half-life ~12.3 years)
- [ ] Test decay rates
- [ ] Document function

**Acceptance Criteria**:
- HTO decay correct
- Other isotopes unaffected
- Tests pass

---

#### Issue #10.3: CLUBB Turbulence Isotope Interface
**Status**: To Do  
**Priority**: Low  
**Labels**: `phase-10`, `clubb`, `optional`  
**Assignees**: TBD  
**Dependencies**: #10.2  
**Estimated Effort**: 7-10 days  

**Description**:
Port CLUBB turbulence isotope interface if using CLUBB parameterization.

**Source**: Multiple CLUBB files in EAMv2-wiso `components/eam/src/physics/clubb/`  
**Destination**: Same in EAMv3-wiso

**Tasks**:
- [ ] Assess CLUBB usage in EAMv3
- [ ] Port CLUBB isotope modules if needed
- [ ] Integrate with water tracers
- [ ] Test CLUBB isotope effects
- [ ] Document implementation

**Acceptance Criteria**:
- CLUBB works with isotopes if enabled
- Turbulent mixing isotope effects correct
- Tests pass

---

#### Issue #10.4: Surface Flux Isotope Interface
**Status**: To Do  
**Priority**: Low  
**Labels**: `phase-10`, `surface`, `optional`  
**Assignees**: TBD  
**Dependencies**: #10.3  
**Estimated Effort**: 5-7 days  

**Description**:
Interface isotopes with surface models (ocean, land) for isotopic surface fluxes.

**Tasks**:
- [ ] Assess surface model coupling needs
- [ ] Interface with ocean evaporation isotopes
- [ ] Interface with land evapotranspiration isotopes
- [ ] Test surface flux isotopes
- [ ] Document interface

**Acceptance Criteria**:
- Surface fluxes provide isotope ratios
- Ocean evaporation fractionation correct
- Land fluxes physically reasonable
- Tests pass

---

## Dependency Graph

```
Phase 0 (Planning)
  └─> Phase 1 (Infrastructure)
       └─> Phase 2 (Rain Evaporation)
            └─> Phase 3 (Vapor-Liquid)
                 └─> Phase 4 (Vapor-Ice)
                      └─> Phase 5 (Freeze/Melt)
                           └─> Phase 6 (Collection)
                                └─> Phase 7 (Sedimentation)
                                     ├─> Phase 8 (Convection)
                                     │    └─> Phase 9 (Conservation & Diagnostics)
                                     │         └─> Phase 10 (Optional)
                                     └─> Phase 9 (can start earlier)
```

**Critical Path**: 0 → 1 → 2 → 3 → 4 → 7 → 9  
**Optional Path**: 9 → 10

---

## Priority Labels

- **Critical**: Must be completed, on critical path
- **High**: Important for completeness, should be done
- **Medium**: Adds value, can be deferred if needed
- **Low**: Optional, nice-to-have

---

## Effort Estimates

- **Total Estimated Effort**: ~190-250 person-days
- **Phase 1**: ~15-20 days
- **Phase 2**: ~20-25 days
- **Phase 3**: ~15-20 days
- **Phase 4**: ~20-25 days
- **Phase 5**: ~15-20 days
- **Phase 6**: ~15-20 days
- **Phase 7**: ~20-25 days
- **Phase 8**: ~12-16 days
- **Phase 9**: ~14-19 days
- **Phase 10**: ~15-22 days (optional)

**Core Implementation** (Phases 1-9): ~155-205 days  
**With Optional Features** (Phase 10): ~170-227 days

---

**End of Issue Templates Document**
