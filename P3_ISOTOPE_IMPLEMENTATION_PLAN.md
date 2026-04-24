# P3 Water Isotope Implementation - Concrete Implementation Plan

**Project**: Port water isotope infrastructure from EAMv2-wiso to EAMv3-wiso  
**Created**: 2026-04-23  
**Status**: Phase 0 Complete → Phase 1 Ready to Start  
**Total Effort**: 190-250 person-days (6-12 months with 1-2 FTEs)

---

## Executive Summary

**Goal**: Enable water isotope simulations (HDO, H₂¹⁸O, H₂¹⁷O, HTO) in E3SM Atmosphere Model v3 using P3 microphysics.

**Approach**: 10-phase incremental implementation with scientific validation at each phase.

**Critical Path** (120-160 days):
```
Phase 1 (Infrastructure) → Phase 2 (Rain Evaporation) → Phase 3 (Vapor-Liquid) 
  → Phase 4 (Vapor-Ice) → Phase 7 (Sedimentation) → Phase 9 (Conservation)
```

**Highest Scientific Priority**: Phase 2 (Rain Evaporation) - controls precipitation isotope composition.

---

## Implementation Roadmap

### **Phase 1: Core Infrastructure Port (15-20 days) - NEXT**

**Priority**: 🔴 CRITICAL PATH  
**Start**: After Phase 0 review approval  
**Scientific Impact**: None (infrastructure only)  
**Testing**: Compilation and unit tests

#### Week 1-2: Module Porting

**Day 1-3**: Issue #1.1 - Port `water_isotopes.F90`
```bash
# Copy fractionation physics module
cp $EAMv2_ROOT/share/util/water_isotopes.F90 \
   $EAMv3_ROOT/share/util/

# Key tasks:
# - Update module dependencies for EAMv3
# - Verify all fractionation functions compile
# - Create unit tests for wiso_alpl(), wiso_alpi(), wiso_akel(), wiso_akci()
```

**Deliverable**: Core fractionation physics available in EAMv3-wiso

---

**Day 4-5**: Issue #1.2 - Port `water_types.F90`
```bash
# Copy water phase type definitions
cp $EAMv2_ROOT/share/util/water_types.F90 \
   $EAMv3_ROOT/share/util/

# Key tasks:
# - Verify water type constants (iwtvap, iwtliq, iwtice, iwtstrain, iwtstsnow)
# - Test master fractionation interface wtype_get_alpha()
```

**Deliverable**: Water type system available

---

**Day 6-7**: Issue #1.3 - Port `water_tracer_vars.F90`
```bash
# Copy configuration variables
cp $EAMv2_ROOT/components/eam/src/physics/cam/water_tracer_vars.F90 \
   $EAMv3_ROOT/components/eam/src/physics/cam/

# Key tasks:
# - Set up namelist parameters
# - Define tracer indexing arrays
```

**Deliverable**: Configuration infrastructure ready

---

#### Week 2-3: Core Framework

**Day 8-14**: Issue #1.4 - Port core `water_tracers.F90` functions

**CRITICAL**: Only port Phase 1 subset (~20% of full module)

```fortran
! Functions to port in Phase 1:
! - Initialization: wtrc_init(), wtrc_register(), wtrc_init_cnst()
! - Rate handling: wtrc_init_rates(), wtrc_add_rates(), wtrc_add_rate()
! - Utilities: wtrc_get_alpha(), wtrc_ratio(), wtrc_ratio_all()
! - Queries: wtrc_is_vap(), wtrc_is_liq(), wtrc_is_ice()

! DEFER to later phases:
! - wtrc_apply_rates() → Phase 2
! - stewart_isoevap() → Phase 2
! - wtrc_sediment() → Phase 7
```

**Deliverable**: Core tracer framework operational

---

**Day 15-17**: Issue #1.5 - Extend `physics_state` structure

```fortran
! File: components/eam/src/physics/cam/physics_types.F90
! Add to physics_state:
!   real(r8), pointer :: qtrc(:,:,:)  ! (ncol, pver, num_tracers)
! Add to physics_ptend:
!   real(r8), pointer :: qtrc_tend(:,:,:)
```

**Deliverable**: Physics state can hold isotope tracers

---

**Day 18-19**: Issue #1.6 - Register constituents

```fortran
! File: components/eam/src/physics/cam/physpkg.F90
! Call during initialization:
call wtrc_register()
! This registers 6 isotopes × 3-5 phases = 18-30 constituents
```

**Deliverable**: Isotope tracers registered in model

---

**Day 20**: Issue #1.7 - Add namelist configuration

```xml
<!-- File: components/eam/bld/namelist_files/namelist_definition.xml -->
<entry id="wisotope" type="logical" category="water_tracers">
  Enable water isotope tracers
</entry>
```

**Deliverable**: Model can be configured via namelist

---

#### Phase 1 Testing & Validation

```bash
# Build with isotopes disabled (backward compatibility)
./case.build

# Build with isotopes enabled
./xmlchange --append CAM_CONFIG_OPTS="-wiso"
./case.build

# Run short test
./case.submit

# Check that model:
# 1. Compiles successfully
# 2. Initializes without errors
# 3. Runs to completion
# 4. Isotope tracers are registered but not yet modified by physics
```

**Success Criteria**:
- ✅ All modules compile without errors
- ✅ Model runs with `wisotope=.true.` in namelist
- ✅ Isotope tracers initialized to prescribed values
- ✅ No impact on standard physics (isotopes not yet active)
- ✅ Unit tests pass for all ported functions

---

### **Phase 2: Rain Evaporation (20-25 days) - HIGHEST PRIORITY**

**Priority**: 🔴 CRITICAL PATH + 🔴 HIGHEST SCIENTIFIC PRIORITY  
**Scientific Impact**: Below-cloud evaporation controls surface precipitation δD/δ18O  
**Testing**: Single-column model + precipitation isotope validation

#### Why Phase 2 is Critical

Rain evaporation below cloud base is the **dominant process** controlling precipitation isotope ratios. Getting this right enables:
- Validation against GNIP (Global Network of Isotopes in Precipitation) observations
- Meaningful comparison to ice core records
- Understanding of moisture recycling

#### Week 4: P3 Modifications

**Day 21-23**: Issue #2.1 - Modify `evaporate_rain()`

```fortran
! File: components/eam/src/physics/p3/eam/micro_p3.F90 (~line 2800)
! Add output parameter:
subroutine evaporate_rain(..., qr_evap_rate)
  real(rtype), intent(out) :: qr_evap_rate(ncol, nlev)
  
  ! Calculate evaporation rate
  qr_evap_rate(i,k) = qr_evap_tend(i,k)
end subroutine
```

**Deliverable**: Rain evaporation rate accessible to isotope interface

---

**Day 24-27**: Issue #2.2 - Create `wtrc_p3_inter_phase2()`

```fortran
! New file: components/eam/src/physics/p3/eam/water_tracers_p3.F90
subroutine wtrc_p3_inter_phase2(state, ptend, pbuf, &
                                qr_evap_rate, T, p, qv, &
                                liqcldf, icecldf, dtime)
  ! Initialize rate arrays
  call wtrc_init_rates(pre_rates_grid, ncol, top_lev)
  
  ! Add rain evaporation rate
  call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                      iwtstrain, iwtvap, iwtstrain, qr_evap_rate)
  
  ! Apply fractionation (simplified for Phase 2)
  call wtrc_apply_rates_phase2(state, ptend, pre_rates_grid, &
                               T, p, qv, liqcldf, dtime)
end subroutine
```

**Deliverable**: Simplified P3 isotope interface for Phase 2

---

#### Week 5-6: Stewart Model Implementation

**Day 28-31**: Issue #2.3 - Port `stewart_isoevap()`

This is the **core scientific function** for rain evaporation fractionation.

```fortran
! Source: EAMv2-wiso water_tracers.F90:3678
! Physics: Stewart (1975) rain evaporation model
!
! Key equations:
!   R_final = ((R_initial - γ·R_vapor)·f^β) + γ·R_vapor
!
! Where:
!   γ = (α_e·h)/(1 - α_e·α_k·(1-h))     ! Effective fractionation
!   β = (1 - α_e·α_k·(1-h))/(α_e·α_k·(1-h))  ! Kinetic exponent
!   α_e = equilibrium fractionation factor
!   α_k = kinetic fractionation factor
!   h = relative humidity
!   f = fraction remaining
```

**Testing**:
```fortran
! Test case: Rain falling through dry layer (RH=70%)
! Expected: δD becomes more negative (depletion in heavy isotope)
! Typical: Δ(δD) ~ -5 to -15 per mil
```

**Deliverable**: Stewart model ported and tested

---

**Day 32-34**: Issue #2.4 - Port `wtrc_equil_time()`

```fortran
! Calculate equilibration fraction based on:
!   - Raindrop size (Marshall-Palmer distribution)
!   - Fall distance (layer thickness)
!   - Diffusion time scale
!   - Temperature and pressure

! Output: fequil ∈ [0,1]
!   0 = no equilibration (large drops, short fall)
!   1 = complete equilibration (small drops, long fall)
```

**Deliverable**: Equilibration calculator ported

---

**Day 35-37**: Issue #2.5 - Port `wtrc_liqvap_equil()`

```fortran
! Iterative solver for liquid-vapor equilibration
! Conserves total water while exchanging isotopes
! Iterates until convergence (typically 3-5 iterations)
```

**Deliverable**: Equilibration solver ported

---

**Day 38-40**: Issue #2.6 - Port `wtrc_apply_rates()` (Phase 2 version)

```fortran
! Simplified version for Phase 2 - only rain evaporation
subroutine wtrc_apply_rates_phase2(...)
  ! Loop over grid points
  do i = 1, ncol
    do k = top_lev, pver
      ! Get rain evaporation rate
      if (qr_evap_rate(i,k) > 0) then
        ! Apply Stewart model
        call stewart_isoevap(ispec, ...)
      endif
    end do
  end do
end subroutine
```

**Deliverable**: Rate application framework (simplified)

---

#### Week 7: Integration & Testing

**Day 41-43**: Issue #2.7 - Integrate with P3

```fortran
! File: components/eam/src/physics/p3/eam/micro_p3_interface.F90
! After P3 microphysics call:
if (do_water_isotopes) then
  call wtrc_p3_inter_phase2(state, ptend, pbuf, &
                            qr_evap_rate, T, p, qv, &
                            liqcldf, icecldf, dtime)
endif
```

**Deliverable**: Isotopes integrated into physics sequence

---

**Day 44-45**: Issue #2.8 - Phase 2 Testing

**Single-Column Model Tests**:

```bash
# Test 1: Tropical convection with dry layer
# Location: Tropical Pacific (0°N, 180°E)
# Conditions: Convective rain falling through RH=60% layer
# Expected: δD ~ -50‰ at cloud base → -65‰ at surface

# Test 2: Arctic stratiform
# Location: Barrow, AK
# Conditions: Light rain, cold temperatures
# Expected: δD ~ -200‰ (very depleted)

# Test 3: Mid-latitude frontal system
# Location: Oklahoma ARM site
# Expected: δD ~ -30 to -80‰ depending on season
```

**Validation Against EAMv2-wiso**:
```bash
# Run identical case in both models
# Compare:
#   - δD_precip at surface
#   - δ18O_precip at surface
#   - Vertical profiles of δD_vapor
```

**Success Criteria**:
- ✅ Rain evaporation decreases δD (expected behavior)
- ✅ Stronger effect at lower relative humidity
- ✅ Results qualitatively match EAMv2-wiso MG2
- ✅ Mass conservation within 0.1%

---

### **Phase 3: Vapor-Liquid Phase Changes (15-20 days)**

**Priority**: 🟠 HIGH  
**Scientific Impact**: Rayleigh distillation during condensation  
**Dependencies**: Phase 2 complete

#### Week 8-9: Implementation

**Day 46-48**: Issues #3.1-3.2 - Modify P3 functions

```fortran
! Add rate outputs to:
!   - cloud_water_autoconversion()
!   - cloud_rain_accretion()
! Both are liquid→liquid (no fractionation, just mass tracking)
```

---

**Day 49-51**: Issue #3.3 - Extract condensation rate

```fortran
! P3 handles condensation implicitly (supersaturation removal)
! Need to extract rate for isotope interface
! Location: within P3 main scheme
```

---

**Day 52-55**: Issue #3.4 - Port `wtrc_vap_distil()`

```fortran
! Rayleigh distillation model for condensation
! Physics:
!   R_vapor_final = R_vapor_initial · f^(α-1)
! where f = fraction remaining after condensation
```

---

**Day 56-60**: Issue #3.5 - Extend `wtrc_apply_rates()`

```fortran
! Add vapor-liquid processes to rate application
! Handle both directions:
!   - Condensation: vapor → liquid (with fractionation)
!   - Evaporation: liquid → vapor (with fractionation)
```

---

**Day 61-65**: Issue #3.6 - Phase 3 Testing

**Test Cases**:
- Cloud formation in rising parcel
- Verify condensation depletes vapor in heavy isotopes
- Check cloud liquid is enriched relative to vapor

---

### **Phase 4: Vapor-Ice Phase Changes (20-25 days)**

**Priority**: 🟠 HIGH  
**Scientific Impact**: Ice fractionation (stronger than liquid)  
**Dependencies**: Phase 3 complete

**Key Physics**: Ice-vapor fractionation is **stronger** than liquid-vapor:
- At -20°C: α(HDO, ice-vapor) ≈ 1.18 vs α(HDO, liq-vapor) ≈ 1.09
- Supersaturation affects kinetic fractionation significantly

#### Week 10-12: Implementation

**Issues #4.1-4.3**: Modify P3 ice functions
- `ice_nucleation()`
- `ice_deposition_sublimation()` - **Bidirectional!**
- `ice_supersat_conservation()`

**Issue #4.4**: Implement ice-vapor fractionation
```fortran
! Use wiso_alpi() for equilibrium
! Use wiso_akci() for kinetic (supersaturation-dependent)
alpha_total = wiso_alpi(ispec, T) * wiso_akci(ispec, Si)
```

**Issue #4.5**: Extend `wtrc_apply_rates()` for ice processes

**Issue #4.6**: Test cirrus clouds, mixed-phase clouds

---

### **Phase 5-6: Freezing, Melting, Collection (30-40 days)**

**Priority**: 🟡 MEDIUM  
**Scientific Impact**: Phase transitions, rime formation  
**Dependencies**: Phase 4 complete

These phases can proceed in parallel once Phase 4 is done.

**Phase 5** (15-20 days): Freezing and melting processes
- Homogeneous freezing (no fractionation)
- Immersion freezing (small kinetic effect)
- Melting (no fractionation at 0°C)

**Phase 6** (15-20 days): Collection and riming
- Ice collecting liquid (riming with kinetic fractionation)
- Wet growth with shedding (complex: partial freezing)
- Self-collection (conservative, no fractionation)

---

### **Phase 7: Sedimentation (20-25 days) - CRITICAL PATH**

**Priority**: 🔴 CRITICAL PATH  
**Scientific Impact**: Layer-by-layer fractionation during fall  
**Dependencies**: Phases 2-6 complete

#### Why Phase 7 is Critical

Sedimentation is critical because precipitation must fall through multiple atmospheric layers, potentially evaporating/sublimating in each layer with cumulative fractionation effects.

#### Week 18-19: Implementation

**Issue #7.1-7.2**: Modify P3 sedimentation functions
- `rain_sedimentation()` - with below-cloud evaporation
- `ice_sedimentation()` - with below-cloud sublimation

**Issue #7.3**: Port `wtrc_sediment()`

```fortran
! Complex implementation:
!   1. CFL-limited sub-stepping for stability
!   2. Layer-by-layer equilibration
!   3. Raindrop size evolution during fall
!   4. Stewart model applied per layer
!   5. Separate clear-sky vs cloudy-sky fractionation

! Critical: Below-cloud evaporation in dry air produces
! strongest isotope signal in precipitation
```

**Issue #7.4**: Integrate with P3 sedimentation scheme

**Issue #7.5**: Test precipitation falling through dry layers

**Success Criteria**:
- ✅ Precipitation isotopes match GNIP observations qualitatively
- ✅ Correct altitude effect (more depleted at higher elevation)
- ✅ Correct temperature effect (more depleted at colder temperatures)

---

### **Phase 8: Convection Interface (12-16 days)**

**Priority**: 🟠 HIGH  
**Scientific Impact**: Convective vs stratiform precipitation distinction  
**Dependencies**: Phase 7 complete (can start in parallel after Phase 1)

#### Convection Isotope Physics

Convective precipitation is isotopically distinct from stratiform:
- **Convective**: Rapid updraft condensation, less equilibration, typically less depleted
- **Stratiform**: Slow condensation, more equilibration, more depleted
- **Amount effect**: Heavier rainfall → more depleted (observed in tropics)

#### Implementation

**Issue #8.1**: Port `wtrc_q1q2_pjr()` for ZM deep convection
```fortran
! File: components/eam/src/physics/cam/zm_conv_intr.F90
! Physics:
!   - Updraft: Rayleigh distillation
!   - Downdraft: Partial evaporation
!   - Convective precipitation: separate tracking
```

**Issue #8.2**: Port `wtrc_shallow()` for shallow convection

**Issue #8.3**: Test tropical convection, midlatitude systems

---

### **Phase 9: Conservation & Diagnostics (14-19 days) - CRITICAL PATH**

**Priority**: 🔴 CRITICAL PATH  
**Scientific Impact**: Ensures physical consistency, enables validation  
**Dependencies**: Phase 8 complete

#### Conservation Checks

```fortran
! Check that for all isotopes:
!   Σ(q_vapor + q_liquid + q_ice + q_rain + q_snow) = constant
!
! Tolerance: < 0.1% error over multi-month simulation
```

#### Diagnostic Output

```fortran
! Add to history output:
!   - δD_vapor(z), δ18O_vapor(z) - vertical profiles
!   - δD_cld, δ18O_cld - cloud water isotopes
!   - δD_ice, δ18O_ice - ice isotopes
!   - δD_precip, δ18O_precip - total precipitation
!   - δD_precip_conv, δ18O_precip_conv - convective only
!   - δD_precip_strat, δ18O_precip_strat - stratiform only
```

#### Global Validation

```bash
# Multi-month global simulation at ne30pg2 resolution
# Compare to:
#   1. GNIP precipitation data (surface stations)
#   2. TES/AIRS satellite vapor isotopes (upper troposphere)
#   3. Aircraft campaigns (e.g., HIPPO, ATom)
#   4. EAMv2-wiso results (consistency check)

# Key metrics:
#   - δD-δ18O meteoric water line (slope should be ~8)
#   - Temperature gradient (δ18O ~ 0.5-0.7‰/°C)
#   - Altitude effect (δ18O ~ -0.1 to -0.3‰ per 100m)
#   - Seasonal cycle amplitude
```

---

### **Phase 10: Optional Physics (15-22 days)**

**Priority**: 🟢 LOW (Optional)  
**Scientific Impact**: Stratospheric isotopes, radioactive tracer  
**Dependencies**: Phase 9 complete

- Methane oxidation (adds H₂O to stratosphere)
- HTO radioactive decay
- CLUBB turbulence interface
- Surface flux isotopes (requires land/ocean coupling)

---

## Resource Requirements

### Personnel

**Minimum Team** (1 FTE, 10-12 months):
- Lead developer with Fortran + E3SM experience

**Recommended Team** (2 FTEs, 6-8 months):
- **Developer 1**: Core implementation (Phases 1-4, 7-9)
- **Developer 2**: P3 modifications (Phases 2-7) + Convection (Phase 8)

**Optional**:
- **Isotope physicist**: Validate fractionation physics
- **Tester**: Design and run validation cases

### Computational Resources

**Development/Testing**:
- HPC allocation: ~10,000 node-hours
- Single-column tests: ~10 node-hours each
- Short global tests: ~100 node-hours each

**Validation** (Phase 9):
- Multi-month global run: ~5,000 node-hours
- Sensitivity tests: ~2,000 node-hours

**Total**: ~15,000-20,000 node-hours

### Data Requirements

**Reference Code**:
- EAMv2-wiso repository access: `/Users/rfiorella/code/E3SM/EAMv2-wiso`

**Validation Data**:
- GNIP precipitation isotope database (free, online)
- TES/AIRS satellite data (NASA archive)
- Aircraft campaign data (varies by campaign)

---

## Risk Management

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| P3 process rates incompatible with isotope framework | Medium | High | Early testing in Phase 2; may need custom rate extraction |
| Iterative solver convergence issues with P3 | Medium | Medium | Test with multiple time steps, may need to adjust convergence criteria |
| Performance degradation > 20% | Low | Medium | Profile code, optimize critical loops, consider selective output |
| Conservation errors > 0.1% | Low | High | Implement strict conservation checks at each phase |

### Scientific Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Results differ significantly from EAMv2-wiso | Medium | High | Detailed process-level comparison, may indicate MG2-P3 physics differences |
| Isotopes don't match observations | Medium | High | Adjust tuning parameters (φ, dkfac), may need process improvements |
| Convective-stratiform partitioning incorrect | Low | Medium | Compare to radar-based precipitation datasets |

### Schedule Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Phase 2 takes longer than estimated | Medium | High | Phase 2 is critical; allocate buffer time, seek expert help if needed |
| Testing reveals fundamental issues | Low | High | Thorough unit testing at each phase to catch issues early |
| Personnel turnover | Low | High | Document everything, code reviews, pair programming |

---

## Testing Strategy (All Phases)

### Unit Tests (Per Phase)

```fortran
! For each function:
program test_wiso_alpl
  ! Test equilibrium fractionation
  real(r8) :: T, alpha, expected
  T = 273.15  ! 0°C
  alpha = wiso_alpl(isphdo, T)
  expected = 1.084  ! Known value for HDO at 0°C
  assert(abs(alpha - expected) < 0.001)
end program
```

### Single-Column Tests (Phases 2+)

```bash
# Tropical Pacific
create_newcase --case test_tropical_scm --compset F2000climo \
               --res T42_T42 --machine pm-cpu

# Arctic
create_newcase --case test_arctic_scm --compset F2000climo \
               --res T42_T42 --machine pm-cpu
```

### Integration Tests (Phases 2+)

```bash
# Short global test (5-10 days)
create_newcase --case test_global_short --compset WCYCL1850 \
               --res ne30pg2_r05_IcoswISC30E3r5 --machine pm-cpu

# Check: no crashes, mass conserved, isotopes reasonable
```

### Validation Tests (Phase 9)

```bash
# Multi-month global run
create_newcase --case validate_global_1year --compset WCYCL1850 \
               --res ne30pg2_r05_IcoswISC30E3r5 --machine pm-cpu

# Post-process:
python scripts/compare_to_gnip.py    # vs precipitation observations
python scripts/compare_to_tes.py     # vs satellite vapor
python scripts/compare_to_eamv2.py   # vs MG2 results
```

---

## Deliverables by Phase

| Phase | Code Deliverables | Documentation | Tests | Timeline |
|-------|------------------|---------------|-------|----------|
| **1** | Core modules ported, compiled | Module documentation | Unit tests | Weeks 1-3 |
| **2** | Rain evaporation working | Stewart model doc | SCM tests | Weeks 4-7 |
| **3** | Condensation fractionation | Process doc | Cloud tests | Weeks 8-9 |
| **4** | Ice processes working | Ice physics doc | Cirrus tests | Weeks 10-12 |
| **5** | Freezing/melting working | Phase change doc | Mixed-phase tests | Weeks 13-14 |
| **6** | Collection working | Riming doc | Collection tests | Weeks 15-16 |
| **7** | Sedimentation working | Sedimentation doc | Precipitation tests | Weeks 18-19 |
| **8** | Convection working | Convection doc | Tropical tests | Weeks 20-22 |
| **9** | Full validation | Validation report | Global comparison | Weeks 23-25 |
| **10** | Optional features | Feature docs | Feature tests | Weeks 26-27 |

---

## Success Metrics

### Phase-Level Success Criteria

| Phase | Metric | Target | Test Method |
|-------|--------|--------|-------------|
| 1 | Compilation | No errors with isotopes on/off | Build test |
| 2 | Rain evaporation δD change | -5 to -15‰ in dry air | SCM test |
| 3 | Condensation fractionation | Cloud δD < vapor δD | SCM test |
| 4 | Ice deposition | Ice δD << vapor δD | Cirrus test |
| 5 | Phase transitions | Mass conserved | Unit test |
| 6 | Collection | No spurious signals | Mixed-phase test |
| 7 | Surface precipitation | Matches GNIP ±20‰ | Global test |
| 8 | Convective/stratiform | Distinct isotope signals | Tropical test |
| 9 | Mass conservation | Error < 0.1% | Global 1-year |

### Final Project Success Criteria

**Scientific**:
- ✅ δD-δ18O meteoric water line slope: 7.5-8.5 (observed: ~8)
- ✅ Temperature gradient: 0.5-0.7‰/°C for δ18O
- ✅ Altitude effect: -0.1 to -0.3‰ per 100m for δ18O
- ✅ Seasonal cycle captured at GNIP sites

**Technical**:
- ✅ Mass conservation < 0.1% error globally
- ✅ Performance overhead < 20%
- ✅ Bit-for-bit reproducible when isotopes disabled
- ✅ Works with all E3SM compsets

**Operational**:
- ✅ Documented in E3SM user guide
- ✅ Example cases created
- ✅ Training materials available
- ✅ Code merged to master branch

---

## Next Steps (Immediate)

### Before Starting Phase 1

1. **Team Review Meeting** (1 day)
   - Review this implementation plan
   - Assign personnel
   - Confirm computational resources
   - Set milestone dates

2. **Development Environment Setup** (2 days)
   - Create feature branch: `feature/p3-isotopes-phase1`
   - Set up testing framework
   - Configure CI/CD if available
   - Clone EAMv2-wiso for reference

3. **Kick-off Phase 1** (Day 1)
   - Begin Issue #1.1: Port water_isotopes.F90
   - Set up weekly status meetings
   - Create shared documentation space

### Weekly Progress Tracking

```bash
# Create project board
git clone https://github.com/E3SM-Project/E3SM.git
cd E3SM
gh project create --title "P3 Water Isotopes"

# Import issues from p3_isotope_kanban.json
python scripts/import_issues.py p3_isotope_kanban.json
```

---

## Communication Plan

**Weekly Status Updates**: Every Friday
- Progress on current phase
- Blockers/issues
- Next week's goals

**Monthly Science Reviews**: End of each month
- Demo current functionality
- Show preliminary results
- Get science team feedback

**Phase Completion Reviews**: After each phase
- Formal review of deliverables
- Decision to proceed to next phase
- Update timeline if needed

---

## Conclusion

This implementation plan provides a concrete, phased approach to porting water isotope infrastructure from EAMv2-wiso (MG2) to EAMv3-wiso (P3). The critical path focuses on rain evaporation first (Phase 2), as this is the dominant process for precipitation isotopes.

**Estimated Timeline**: 6-8 months with 2 FTEs, 10-12 months with 1 FTE

**Key to Success**:
1. Thorough testing at each phase before proceeding
2. Early validation against observations (Phase 2)
3. Maintain communication with science team
4. Document everything as you go

**Questions?** Review with team and adjust plan as needed before beginning Phase 1.

---

**Document Version**: 1.0  
**Last Updated**: 2026-04-23  
**Next Review**: After Phase 1 completion

