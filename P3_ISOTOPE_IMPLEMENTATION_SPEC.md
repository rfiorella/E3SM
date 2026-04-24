# **Water Isotope Infrastructure Port to P3 Microphysics**
## **Implementation Specification**

**Document Version**: 1.0  
**Date**: 2026-04-21  
**Status**: Awaiting Phase 0 review

---

## **Project Overview**

Port water isotope/tracer infrastructure from EAMv2-wiso (Morrison-Gettelman microphysics) to EAMv3-wiso (P3 microphysics) to enable simulation of stable water isotopes (HDO, H218O, H217O) and tritium (HTO) in the E3SM atmosphere model.

---

## **Isotope Species to Track**

All water isotope species from EAMv2-wiso:

1. **H2O** - Bulk water (standard tracer)
2. **H216O** - Standard oxygen-16 water (nearly identical to bulk)
3. **HDO** - Heavy hydrogen (deuterium)
4. **H218O** - Heavy oxygen-18
5. **H217O** - Heavy oxygen-17 (with mass-dependent scaling)
6. **HTO** - Tritium (with radioactive decay)

---

## **Water Phase Types**

Must track isotopes in all water phases:

1. **Water vapor** (iwtvap)
2. **Cloud liquid** (iwtliq)
3. **Cloud ice** (iwtice)
4. **Stratiform rain** (iwtstrain)
5. **Stratiform snow** (iwtstsnow)
6. **Convective rain** (iwtcvrain) - *requires convection interface*
7. **Convective snow** (iwtcvsnow) - *requires convection interface*

---

## **Phased Implementation Plan**

### **Phase 0: Infrastructure Setup** ✓

**Objective**: Document existing systems and create implementation roadmap

**Tasks**:
- [x] Catalog P3 microphysics functions requiring isotope support
- [x] Document EAMv2-wiso isotope infrastructure
- [x] Create detailed specification document
- [ ] Review specification with science team
- [ ] Finalize design decisions

**Deliverables**:
- This specification document
- Function mapping between MG2 and P3 processes

---

### **Phase 1: Core Infrastructure Port**

**Objective**: Port essential isotope infrastructure modules from EAMv2-wiso to EAMv3-wiso

#### **Files to Port** (from `/Users/rfiorella/code/E3SM/EAMv2-wiso`):

1. **`share/util/water_isotopes.F90`** → EAMv3-wiso
   - Core fractionation calculations
   - Equilibrium fractionation factors: `wiso_alpl()`, `wiso_alpi()`
   - Kinetic fractionation factors: `wiso_akel()`, `wiso_akci()`
   - Isotope species definitions and constants
   - Diffusivity ratios, standard ratios

2. **`share/util/water_types.F90`** → EAMv3-wiso
   - Water phase type definitions
   - Master fractionation interface: `wtype_get_alpha()`

3. **`components/eam/src/physics/cam/water_tracer_vars.F90`** → EAMv3-wiso
   - Configuration variables
   - Namelist parameters
   - Tracer indexing arrays

4. **`components/eam/src/physics/cam/water_tracers.F90`** → EAMv3-wiso
   - Core water tracer module (~7,800 lines)
   - Process rate framework
   - Key functions needed immediately:
     - `wtrc_init()`, `wtrc_register()`, `wtrc_init_cnst()`
     - `wtrc_init_rates()`, `wtrc_add_rates()`, `wtrc_add_rate()`
     - `wtrc_get_alpha()` - fractionation factor wrapper
     - `wtrc_ratio()`, `wtrc_ratio_all()`
     - Conservation checking functions

#### **Integration Tasks**:

1. Add isotope tracers to constituent registry
2. Extend `physics_state` and `physics_ptend` for isotopes
3. Update `physpkg.F90` to initialize water tracers
4. Create namelist entries for isotope configuration

#### **Testing**:

- Compile test with isotopes disabled (backward compatibility)
- Compile test with isotopes enabled
- Verify tracer registration and initialization

#### **Deliverables**:

- Water isotope infrastructure compiles in EAMv3-wiso
- Tracers initialized but not yet modified by physics

---

### **Phase 2: Stratiform Rain Evaporation (Critical Path)**

**Objective**: Implement isotope fractionation for rain evaporation - the most important process for precipitation isotope composition

**Scientific Priority**: **HIGHEST** - Rain evaporation below cloud base is the dominant control on precipitation isotope ratios reaching the surface

#### **Functions to Modify**:

1. **`evaporate_rain`** (micro_p3.F90)
   - Extract evaporation rate (kg/kg/s)
   - Provide to isotope interface

#### **New Functions to Create**:

2. **`wtrc_p3_inter_phase2()`** (new file: `water_tracers_p3.F90`)
   - Simplified P3 interface for Phase 2
   - Collects only rain evaporation rate
   - Calls isotope tendency calculator
   - Based on `wtrc_mg_inter()` structure

#### **Functions to Port from EAMv2-wiso**:

3. **`stewart_isoevap()`** (water_tracers.F90:3678)
   - Stewart (1975) rain evaporation model
   - Includes:
     - Rayleigh distillation during evaporation
     - Kinetic fractionation (humidity-dependent)
     - Equilibrium fractionation (temperature-dependent)
     - Marshall-Palmer raindrop size distribution
     - Partial equilibration calculation

4. **`wtrc_equil_time()`** (water_tracers.F90)
   - Calculates fraction of rain experiencing equilibration
   - Function of droplet size, fall distance, temperature

5. **`wtrc_liqvap_equil()`** (water_tracers.F90)
   - Iterative liquid-vapor equilibration solver
   - Exchanges isotopes between rain and vapor

6. **`wtrc_apply_rates()`** (water_tracers.F90:2068)
   - Core isotope tendency calculator
   - Iterative solver for process rates
   - Applies fractionation factors
   - **Initially simplified to handle only rain evaporation**

#### **Process Flow**:

```
P3 micro_p3_tend()
  ↓
P3 evaporate_rain() → outputs evaporation rate
  ↓
wtrc_p3_inter_phase2() → collects rate
  ↓
wtrc_apply_rates() → applies fractionation
  ↓
stewart_isoevap() → calculates isotope partitioning
  ↓
Updates vapor and rain isotope tracers
```

#### **Key Scientific Parameters**:

- Equilibrium fractionation: Horita & Wesolowski (1994)
- Kinetic fractionation: Diffusivity ratios (difrm)
- Rain droplet size: Marshall-Palmer distribution
- Tuning parameter φ = 0.9 (from Bony et al. 2008)

#### **Testing**:

- Single-column test with prescribed rain evaporation
- Compare against MG2 isotope behavior
- Verify mass conservation
- Check δD and δ18O values against observations

#### **Deliverables**:

- Rain evaporation isotope fractionation working
- Precipitation isotope output to history files
- Validation against EAMv2-wiso MG2 results

---

### **Phase 3: Vapor-Liquid Phase Changes**

**Objective**: Add fractionation for condensation and cloud water processes

#### **Functions to Modify/Extend**:

1. **`cloud_water_autoconversion`** (micro_p3.F90)
   - Extract autoconversion rate (cloud → rain)
   - No fractionation (liquid → liquid) but mass tracking

2. **`cloud_rain_accretion`** (micro_p3.F90)
   - Extract accretion rate (cloud → rain)
   - No fractionation but mass tracking

#### **Functions to Port/Extend**:

3. **`wtrc_apply_rates()`** - Extend Phase 2 version
   - Add condensation/evaporation with equilibrium fractionation
   - Handle cloud supersaturation removal
   - Rayleigh distillation for condensation

4. **`wtrc_vap_distil()`** (water_tracers.F90)
   - Rayleigh distillation routine
   - Used for condensation fractionation

#### **P3 Processes to Interface**:

- Cloud supersaturation adjustment (implicit condensation)
- Cloud evaporation
- Autoconversion (mass tracking only)
- Accretion (mass tracking only)

#### **Update `wtrc_p3_inter`**:

- Expand from Phase 2 version
- Add liquid process rates
- Pre-sedimentation rate collection

#### **Testing**:

- Cloud formation test case
- Verify condensation fractionation
- Check cloud water isotope ratios

#### **Deliverables**:

- Condensation/evaporation fractionation working
- Cloud water isotope composition tracked

---

### **Phase 4: Vapor-Ice Phase Changes**

**Objective**: Add fractionation for ice nucleation, deposition, and sublimation

**Scientific Complexity**: **HIGH** - Ice fractionation has different physics than liquid, includes supersaturation dependence

#### **Functions to Modify/Extend**:

1. **`ice_nucleation`** (micro_p3.F90)
   - Extract nucleation rate (vapor → ice)
   - Provide temperature and supersaturation

2. **`ice_deposition_sublimation`** (micro_p3.F90)
   - Extract bidirectional rates
   - Critical: distinguish deposition vs sublimation

3. **`ice_supersat_conservation`** (micro_p3.F90)
   - Extract supersaturation removal rate
   - Similar to deposition

#### **Fractionation Physics to Implement**:

4. **Ice-vapor equilibrium fractionation**
   - `wiso_alpi()` from water_isotopes.F90
   - Based on Merlivat & Nief (1967), Majoube (1971)
   - Temperature-dependent

5. **Kinetic fractionation for ice**
   - `wiso_akci()` from water_isotopes.F90
   - Supersaturation-dependent (important!)
   - Different from liquid kinetic effects

#### **Update `wtrc_p3_inter`**:

- Add ice process rates
- Handle bidirectional vapor-ice exchange
- Apply ice-specific fractionation

#### **Testing**:

- Cirrus cloud case (pure ice)
- Mixed-phase cloud case
- Verify depositional fractionation
- Check ice supersaturation effects

#### **Deliverables**:

- Ice deposition/sublimation fractionation working
- Ice nucleation fractionation implemented
- Ice water isotope composition tracked

---

### **Phase 5: Freezing and Melting**

**Objective**: Add fractionation for liquid-ice transitions

#### **Functions to Modify/Extend**:

1. **`cldliq_immersion_freezing`** (micro_p3.F90)
   - Extract freezing rate
   - Kinetic fractionation during freezing

2. **`rain_immersion_freezing`** (micro_p3.F90)
   - Extract rain freezing rate
   - Similar to cloud liquid freezing

3. **`homogeneous_freezing`** (micro_p3.F90)
   - Extract homogeneous freezing rate
   - Minimal fractionation (rapid process)

4. **`ice_melting`** (micro_p3.F90)
   - Extract melting rate
   - Minimal/no fractionation at 0°C

5. **`ice_complete_melting`** (micro_p3.F90)
   - Complete phase transition
   - Simple mass transfer

#### **Fractionation Physics**:

- Freezing: Small kinetic effects (temperature-dependent)
- Melting: Minimal fractionation (equilibrium at Tm)
- Fast processes: Conservative (no fractionation)

#### **Update `wtrc_p3_inter`**:

- Add freezing/melting rates
- Distinguish fast vs slow processes
- Apply appropriate fractionation factors

#### **Testing**:

- Freezing rain case
- Mixed-phase cloud evolution
- Verify phase transition isotope effects

#### **Deliverables**:

- Freezing/melting fractionation working
- Correct isotope transfer between phases

---

### **Phase 6: Collection and Riming**

**Objective**: Handle isotope mass conservation during collection processes

#### **Functions to Modify/Extend**:

1. **`ice_cldliq_collection`** (micro_p3.F90)
   - Riming of cloud liquid onto ice
   - Kinetic fractionation during impact

2. **`ice_cldliq_wet_growth`** (micro_p3.F90)
   - Partial freezing with shedding
   - Complex: some liquid shed, some frozen

3. **`ice_rain_collection`** (micro_p3.F90)
   - Rain collected by ice
   - Similar to riming

4. **`ice_self_collection`** (micro_p3.F90)
   - Aggregation (ice + ice)
   - Conservative (no fractionation)

5. **`rain_self_collection`** (micro_p3.F90)
   - Coalescence (rain + rain)
   - Conservative

6. **`droplet_self_collection`** (micro_p3.F90)
   - Cloud droplet collection
   - Conservative

#### **Fractionation Physics**:

- Riming: Kinetic effects during liquid-ice contact
- Wet growth: Partial equilibration
- Self-collection: Conservative (same phase)

#### **Update `wtrc_p3_inter`**:

- Add collection rates
- Handle rime mass tracking
- Bergeron process (vapor transfer from liquid to ice)

#### **Testing**:

- Riming case study
- Verify mass conservation in collection
- Check rime isotope composition

#### **Deliverables**:

- Collection process isotopes working
- Rime mass isotope tracking
- Mass conservation verified

---

### **Phase 7: Sedimentation with Fractionation**

**Objective**: Implement vertical transport with below-cloud fractionation

#### **Functions to Modify/Extend**:

1. **`rain_sedimentation`** (micro_p3.F90)
   - Vertical transport
   - **Critical**: Below-cloud evaporation with Stewart model
   - Layer-by-layer fractionation

2. **`ice_sedimentation`** (micro_p3.F90)
   - Ice/snow vertical transport
   - Below-cloud sublimation fractionation

3. **`cloud_sedimentation`** (micro_p3.F90)
   - Cloud droplet fall (usually small)
   - Evaporation equilibration

#### **Functions to Port/Extend**:

4. **`wtrc_sediment()`** (water_tracers.F90:3315)
   - Full sedimentation isotope handler
   - CFL-limited sub-stepping
   - Layer-by-layer equilibration
   - Stewart model for rain evaporation during fall

#### **Critical Implementation Details**:

- Must handle partial equilibration during fall
- Raindrop size evolves as evaporation occurs
- Calculate equilibration fraction per layer
- Apply fractionation only to evaporating fraction

#### **Process Flow for Rain Sedimentation**:

```
For each layer k (top to bottom):
  1. Calculate rain flux into layer
  2. Calculate evaporation in layer
  3. Determine raindrop size from flux
  4. Calculate equilibration fraction (wtrc_equil_time)
  5. Apply Stewart model to evaporating fraction
  6. Calculate rain flux out of layer
  7. Update vapor and rain isotopes
```

#### **Update `wtrc_p3_inter`**:

- Interface sedimentation rates
- Call wtrc_sediment() for each precipitating species
- Handle cloud fraction effects (clear vs cloudy)

#### **Testing**:

- Precipitation falling through dry layer
- Verify δD decrease with evaporation
- Check vertical isotope profiles
- Compare to observations (e.g., GNIP data)

#### **Deliverables**:

- Sedimentation with fractionation working
- Below-cloud isotope evolution correct
- Surface precipitation isotopes validated

---

### **Phase 8: Convection Interface**

**Objective**: Track convective rain and snow isotopes separately from stratiform

**Requirement**: P3 is stratiform-only; convection is handled by separate schemes (ZM deep, UW shallow)

#### **New Files/Functions Required**:

1. **Extend existing convection interfaces**:
   - `zm_conv_intr.F90` - Zhang-McFarlane deep convection
   - `convect_shallow.F90` or `uwshcu.F90` - Shallow convection

2. **Port from EAMv2-wiso**:
   - `wtrc_q1q2_pjr()` (water_tracers.F90:4960) - ZM convection isotopes
     - Updraft condensation with fractionation
     - Downdraft evaporation with fractionation
     - Detrainment of condensate
     - Convective precipitation
   
   - `wtrc_shallow()` (water_tracers.F90:6234) - Shallow convection
     - Similar structure to deep convection
     - Handles shallow cumulus isotopes

3. **Add convective precipitation types**:
   - `iwtcvrain` - Convective rain isotope tracker
   - `iwtcvsnow` - Convective snow isotope tracker

#### **Integration Points**:

- Convection runs before stratiform microphysics
- Convective detrainment provides condensate to P3
- Need to track stratiform vs convective precipitation separately

#### **Output Requirements**:

- `PRECRC_H218O`, `PRECRC_HDO` - Convective rain isotopes
- `PRECSC_H218O`, `PRECSC_HDO` - Convective snow isotopes
- `PRECSL_H218O`, `PRECSL_HDO` - Stratiform rain isotopes (from P3)
- `PRECSS_H218O`, `PRECSS_HDO` - Stratiform snow isotopes (from P3)

#### **Testing**:

- Tropical convection case
- Compare convective vs stratiform isotope signatures
- Verify total precipitation = convective + stratiform

#### **Deliverables**:

- Convective precipitation isotopes tracked
- Separate output for convective vs stratiform
- Total water isotope budget closed

---

### **Phase 9: Conservation and Diagnostics**

**Objective**: Ensure mass conservation and provide comprehensive isotope diagnostics

#### **Conservation Functions to Implement**:

1. **Extend P3 conservation functions** for isotopes:
   - `water_vapor_conservation` - vapor isotope mass
   - `cloud_water_conservation` - cloud liquid isotope mass
   - `rain_water_conservation` - rain isotope mass
   - `ice_water_conservation` - ice isotope mass

2. **Port from EAMv2-wiso**:
   - `wtrc_check_h2o()` (water_tracers.F90:1388) - Mass conservation checks
   - `wtrc_adjust_h2o()` (water_tracers.F90:1513) - Fix small negative values
   - `wtrc_mass_fixer()` (water_tracers.F90:6677) - Reset standard tracer to bulk

#### **Diagnostic Functions**:

3. **Isotope ratio calculations**:
   - Convert tracer mass to δD, δ18O, δ17O notation
   - δ notation: δ = (R_sample/R_standard - 1) × 1000 [per mil]

4. **Output to history files**:
   - Atmospheric columns: δD_vapor, δ18O_vapor, δ17O_vapor
   - Cloud phases: δD_cld, δ18O_cld, δD_ice, δ18O_ice
   - Precipitation: δD_precip, δ18O_precip (total, convective, stratiform)

5. **Budget diagnostics**:
   - Column-integrated isotope mass
   - Process-level isotope tendencies
   - Fractionation diagnostics

#### **Port from EAMv2-wiso**:

- `wtrc_output_precip()` (water_tracers.F90:4161) - Precipitation isotope output
- `wtrc_setup_diag()` (water_tracers.F90:1356) - Diagnostic setup
- Add history fields via `addfld()` and `add_default()`

#### **Testing**:

- Global mass conservation check (multi-month run)
- Compare isotope budgets to EAMv2-wiso
- Verify δD-δ18O relationship (meteoric water line)

#### **Deliverables**:

- Mass conservation enforced
- Comprehensive isotope diagnostics
- Output compatible with EAMv2-wiso format

---

### **Phase 10: Additional Physics (Optional/Future)**

**Objective**: Add remaining isotope physics processes

#### **Additional Processes**:

1. **Methane oxidation**:
   - `wtrc_chem_ch4ox_tend()` (water_tracers.F90:6751)
   - Adds water vapor from CH4 + O2 → H2O + CO2
   - Important for stratospheric water isotopes

2. **Radioactive decay**:
   - `wtrc_rad_decay()` (water_tracers.F90:6832)
   - HTO decay (half-life ~12.3 years)
   - Only affects tritium tracer

3. **CLUBB turbulence interface**:
   - Isotope-aware turbulent mixing
   - PDF-based condensation with fractionation
   - Port CLUBB isotope modules if using CLUBB

4. **Surface flux isotopes**:
   - Ocean evaporation with fractionation
   - Land evapotranspiration isotopes
   - Requires surface model coupling

#### **Testing**:

- Stratospheric HDO from methane oxidation
- HTO decay rates
- Surface flux validation

#### **Deliverables**:

- Complete isotope physics package
- Ready for scientific applications

---

## **Key Files and Their Locations**

### **Source Files** (EAMv2-wiso → EAMv3-wiso)

| EAMv2-wiso Path | Purpose | Phase |
|----------------|---------|-------|
| `share/util/water_isotopes.F90` | Core fractionation | 1 |
| `share/util/water_types.F90` | Phase type definitions | 1 |
| `components/eam/src/physics/cam/water_tracer_vars.F90` | Configuration | 1 |
| `components/eam/src/physics/cam/water_tracers.F90` | Main module | 1+ |
| `components/eam/src/physics/cam/micro_mg_cam.F90` | MG2 interface (reference) | 2 |
| `components/eam/src/physics/cam/zm_conv_intr.F90` | Deep convection | 8 |
| `components/eam/src/physics/cam/uwshcu.F90` | Shallow convection | 8 |

### **New Files to Create** (EAMv3-wiso)

| File | Purpose | Phase |
|------|---------|-------|
| `components/eam/src/physics/p3/eam/water_tracers_p3.F90` | P3 isotope interface | 2 |
| `components/eam/src/physics/p3/eam/micro_p3_isotopes.F90` | P3 isotope utilities | 3+ |

### **Files to Modify** (EAMv3-wiso)

| File | Modifications | Phase |
|------|---------------|-------|
| `components/eam/src/physics/p3/eam/micro_p3.F90` | Add process rate outputs | 2+ |
| `components/eam/src/physics/p3/eam/micro_p3_interface.F90` | Call isotope interface | 2 |
| `components/eam/src/physics/cam/physpkg.F90` | Initialize isotopes | 1 |
| `components/eam/src/physics/cam/physics_types.F90` | Extend state structure | 1 |

---

## **P3 Functions Requiring Isotope Infrastructure**

### **Priority 1: Core Phase Change Processes** (Require Fractionation)

#### **A. Vapor ↔ Liquid Processes**

1. **`cloud_water_autoconversion`** - Converts cloud water → rain
2. **`evaporate_rain`** - Rain evaporation → vapor (**CRITICAL - Stewart model**)
3. **`cloud_rain_accretion`** - Cloud water collection by rain

#### **B. Vapor ↔ Ice Processes**

4. **`ice_nucleation`** - Vapor → ice (**CRITICAL - kinetic fractionation**)
5. **`ice_deposition_sublimation`** - Vapor ↔ ice (**CRITICAL - bidirectional**)
6. **`ice_supersat_conservation`** - Removes supersaturation by deposition

#### **C. Liquid ↔ Ice Processes**

7. **`ice_melting`** - Ice → liquid
8. **`ice_complete_melting`** - Complete ice → liquid
9. **`cldliq_immersion_freezing`** - Cloud liquid → ice
10. **`rain_immersion_freezing`** - Rain → ice
11. **`homogeneous_freezing`** - Cloud liquid → ice (rapid)

#### **D. Collection/Riming Processes**

12. **`ice_cldliq_collection`** - Ice collects cloud liquid (riming)
13. **`ice_cldliq_wet_growth`** - Wet growth regime
14. **`ice_rain_collection`** - Ice collects rain

#### **E. Ice-Ice and Rain-Rain Processes**

15. **`ice_self_collection`** - Ice aggregation (conservative)
16. **`rain_self_collection`** - Rain coalescence (conservative)
17. **`droplet_self_collection`** - Cloud droplet collection (conservative)

### **Priority 2: Sedimentation Processes**

18. **`rain_sedimentation`** - Vertical transport with evaporation
19. **`ice_sedimentation`** - Ice transport with sublimation
20. **`cloud_sedimentation`** - Cloud droplet fall
21. **`generalized_sedimentation`** - Framework for sedimentation

### **Priority 3: Supporting Functions**

22. **`get_cloud_dsd2`** - Cloud droplet size distribution
23. **`get_rain_dsd2`** - Rain size distribution
24. **`update_prognostic_ice`** - Updates ice state
25. **`update_prognostic_liquid`** - Updates liquid state
26-32. **Conservation functions** - Mass conservation for all phases

### **Priority 4: Interface and Driver Functions**

33. **`p3_main`** - Main P3 driver
34. **`p3_main_part1`** - First calculation segment
35. **`p3_main_part2`** - Second calculation segment
36. **`p3_main_part3`** - Sedimentation segment
37. **`micro_p3_tend`** - Interface between EAM and P3

---

## **Critical New Infrastructure Required**

### **1. Main Interface Function**

**`wtrc_p3_inter()`** - New function modeled after `wtrc_mg_inter()` (water_tracers.F90:1812)
- Collects all P3 process rates
- Maps P3 processes to water types (vapor, liquid, ice, rain, snow)
- Calls `wtrc_apply_rates()` with fractionation
- Handles pre-sedimentation, sedimentation, and post-sedimentation processes

### **2. Process Rate Collection**

Extend P3 to output process rates (kg/kg/s) for:
- Each vapor source/sink
- Each liquid source/sink  
- Each ice source/sink
- Each rain source/sink
- Temperature and cloud fraction for fractionation calculations

### **3. Isotope State Variables**

Add prognostic variables for each P3 water species:
- `qtrcr_vap` - isotope tracer in vapor
- `qtrcr_cld` - isotope tracer in cloud liquid
- `qtrcr_ice` - isotope tracer in cloud ice
- `qtrcr_rain` - isotope tracer in rain
- Multiple isotopes (HDO, H218O, H217O, HTO)

### **4. Precipitation Isotopes**

Track isotopic composition of:
- Surface rain rate and isotope ratio
- Surface snow rate and isotope ratio

---

## **Implementation Strategy Summary**

Based on the MG2 implementation in EAMv2-wiso, the recommended approach is:

1. **Phase 1**: Port core infrastructure (water_isotopes, water_types, water_tracers modules)
2. **Phase 2**: Implement rain evaporation (Stewart model) - highest priority for precipitation isotopes
3. **Phase 3**: Add vapor-liquid condensation/evaporation
4. **Phase 4**: Add vapor-ice deposition/sublimation
5. **Phase 5**: Add freezing/melting transitions
6. **Phase 6**: Add collection and riming processes
7. **Phase 7**: Implement sedimentation with fractionation
8. **Phase 8**: Add convection interface (ZM deep, UW shallow)
9. **Phase 9**: Add conservation checks and diagnostics
10. **Phase 10**: Add additional physics (methane oxidation, HTO decay, CLUBB, surface fluxes)

---

## **Key Dependencies from EAMv2-wiso**

Required modules/functions to port:
- `water_isotopes.F90` - Core fractionation physics
- `water_types.F90` - Water phase type definitions  
- `water_tracers.F90` - Process rate framework
- `water_tracer_vars.F90` - Configuration variables
- `wtrc_apply_rates()` - Main isotope tendency calculator
- `wtrc_liqvap_equil()` - Equilibration solver
- `stewart_isoevap()` - Rain evaporation model
- `wtrc_equil_time()` - Partial equilibration calculator
- `wtrc_get_alpha()` - Fractionation factor lookup

---

## **Success Criteria by Phase**

| Phase | Success Metric |
|-------|----------------|
| 1 | Code compiles with isotopes registered |
| 2 | Rain evaporation δD matches Stewart theory |
| 3 | Condensation fractionation produces reasonable cloud δD |
| 4 | Ice clouds show expected depletion patterns |
| 5 | Phase transitions conserve mass correctly |
| 6 | Collection processes don't create spurious signals |
| 7 | Precipitation isotopes match observations qualitatively |
| 8 | Convective vs stratiform separation works |
| 9 | Global mass conservation < 0.1% error |
| 10 | Full physics matches EAMv2-wiso behavior |

---

## **Dependencies and Risks**

### **Technical Dependencies**:

- EAMv2-wiso isotope code must be portable to EAMv3
- P3 process rates must be accessible for isotope calculations
- Physics timestep coordination between P3 and isotope solver

### **Scientific Risks**:

- P3 process representation may differ from MG2, affecting isotope behavior
- Iterative solver may not converge with P3 process rates
- Below-cloud fractionation may be more complex in P3 due to ice-rain interactions

### **Resource Risks**:

- Computational cost: isotopes add ~5-6 tracers × 3-4 phases = 15-24 additional prognostic variables
- Memory: Large rate arrays for process tracking
- I/O: Additional history output for diagnostics

---

## **Testing Strategy** (To Be Developed)

*Detailed testing plan to be developed after Phase 0 review*

Key test categories identified:
- Unit tests (individual functions)
- Single-column model tests
- Idealized cases (tropical, Arctic, stratosphere)
- Validation against observations (GNIP, aircraft)
- Comparison with EAMv2-wiso MG2 results

---

## **Questions for Phase 0 Review**

Before proceeding to Phase 1, please confirm:

1. **Scope agreement**: Is the phased approach with rain evaporation first acceptable?

2. **Convection timing**: Should Phase 8 (convection) be moved earlier or kept late?

3. **Simplifications**: Any physics processes we can skip initially?

4. **Code structure**: Create new files (`water_tracers_p3.F90`) or extend existing?

5. **Validation data**: What observational datasets should we target for validation?

6. **Timeline**: What is the target completion date for each phase?

7. **Resources**: Who will be involved in coding, testing, and validation?

---

## **Next Steps**

After Phase 0 approval:

1. **Phase 1 kickoff**: Begin porting core infrastructure modules
2. **Code review process**: Establish review workflow
3. **Testing framework**: Set up unit tests and test cases
4. **Documentation**: Create developer guide for isotope implementation

---

## **References**

### **Scientific References**:

- **Stewart (1975)**: Rain evaporation isotope model
- **Horita & Wesolowski (1994)**: Liquid-vapor equilibrium fractionation
- **Merlivat & Nief (1967)**: Ice-vapor equilibrium fractionation
- **Majoube (1971)**: Temperature dependence of fractionation
- **Bony et al. (2008)**: Below-cloud relative humidity effects

### **Code References**:

- **EAMv2-wiso**: `/Users/rfiorella/code/E3SM/EAMv2-wiso`
- **EAMv3-wiso**: `/Users/rfiorella/code/E3SM/EAMv3-wiso`
- **P3 scheme**: Morrison & Milbrandt (2015), Milbrandt & Morrison (2016)

---

**End of Specification Document**
