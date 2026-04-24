# **P3 Microphysics to Water Isotope Function Mapping**

**Document Version**: 1.0  
**Date**: 2026-04-21  
**Related**: P3_ISOTOPE_IMPLEMENTATION_SPEC.md

---

## **Table of Contents**

1. [Overview](#overview)
2. [MG2 to P3 Process Mapping](#mg2-to-p3-process-mapping)
3. [P3 Functions Detailed Mapping](#p3-functions-detailed-mapping)
4. [Water Type Transition Matrix](#water-type-transition-matrix)
5. [Fractionation Factor Requirements](#fractionation-factor-requirements)
6. [EAMv2-wiso Function Dependencies](#eamv2-wiso-function-dependencies)

---

## **Overview**

This document provides detailed mappings between:
- MG2 microphysics processes (EAMv2-wiso) and P3 microphysics processes (EAMv3-wiso)
- P3 functions and required isotope infrastructure
- Water phase transitions and fractionation requirements
- EAMv2-wiso functions and their P3 equivalents

---

## **MG2 to P3 Process Mapping**

### **Vapor-Liquid Processes**

| MG2 Process | MG2 Rate Variable | P3 Equivalent Function | P3 Implementation Phase | Notes |
|-------------|-------------------|------------------------|-------------------------|-------|
| Condensation (supersaturation removal) | `qcreso` | Part of P3 saturation adjustment | Phase 3 | Implicit in P3 |
| Cloud evaporation | `qcreso` (negative) | Cloud evaporation in P3 | Phase 3 | Need to extract rate |
| Autoconversion | `prco` | `cloud_water_autoconversion` | Phase 3 | Cloud → rain |
| Accretion | `prao` | `cloud_rain_accretion` | Phase 3 | Cloud collected by rain |
| Rain evaporation | `preo` | `evaporate_rain` | Phase 2 | **CRITICAL - Stewart model** |

### **Vapor-Ice Processes**

| MG2 Process | MG2 Rate Variable | P3 Equivalent Function | P3 Implementation Phase | Notes |
|-------------|-------------------|------------------------|-------------------------|-------|
| Ice nucleation | `mnuccco`, `mnuccto` | `ice_nucleation` | Phase 4 | Heterogeneous + homogeneous |
| Ice deposition | `qireso` (positive) | `ice_deposition_sublimation` | Phase 4 | Vapor → ice |
| Ice sublimation | `qireso` (negative) | `ice_deposition_sublimation` | Phase 4 | Ice → vapor |
| Deposition (supersat removal) | `pcmei` (positive) | `ice_supersat_conservation` | Phase 4 | Sets Si = 100% |

### **Liquid-Ice Processes**

| MG2 Process | MG2 Rate Variable | P3 Equivalent Function | P3 Implementation Phase | Notes |
|-------------|-------------------|------------------------|-------------------------|-------|
| Immersion freezing (cloud) | `mnuccco` | `cldliq_immersion_freezing` | Phase 5 | Cloud liquid → ice |
| Contact freezing | `mnuccto` | Part of `cldliq_immersion_freezing` | Phase 5 | Combined in P3 |
| Homogeneous freezing | `homoo` | `homogeneous_freezing` | Phase 5 | Rapid freezing |
| Rain freezing | `mnuccro` | `rain_immersion_freezing` | Phase 5 | Rain → ice/snow |
| Ice melting | `melto` | `ice_melting`, `ice_complete_melting` | Phase 5 | Ice → liquid |
| Snow melting | `meltso` | Part of ice melting in P3 | Phase 5 | P3 doesn't separate snow type |

### **Collection Processes**

| MG2 Process | MG2 Rate Variable | P3 Equivalent Function | P3 Implementation Phase | Notes |
|-------------|-------------------|------------------------|-------------------------|-------|
| Ice-cloud collection | `msacwio` | `ice_cldliq_collection` | Phase 6 | Riming |
| Ice-cloud wet growth | Part of `msacwio` | `ice_cldliq_wet_growth` | Phase 6 | Shedding occurs |
| Ice-rain collection | `praio`, `prcio` | `ice_rain_collection` | Phase 6 | Rain collected by ice |
| Snow-cloud collection | `psacwso` | Part of ice-cloud in P3 | Phase 6 | P3 unified ice |
| Snow-rain collection | `pracso` | Part of ice-rain in P3 | Phase 6 | P3 unified ice |
| Bergeron process (liquid) | `bergo` | Implicit in P3 supersaturation | Phase 6 | Vapor transfer |
| Bergeron process (ice) | `bergso` | Implicit in P3 supersaturation | Phase 6 | Vapor transfer |

### **Sedimentation Processes**

| MG2 Process | MG2 Variable | P3 Equivalent Function | P3 Implementation Phase | Notes |
|-------------|--------------|------------------------|-------------------------|-------|
| Rain sedimentation | `fr` (fall velocity) | `rain_sedimentation` | Phase 7 | With evaporation |
| Snow sedimentation | `fs` (fall velocity) | `ice_sedimentation` | Phase 7 | With sublimation |
| Cloud sedimentation | Usually negligible | `cloud_sedimentation` | Phase 7 | Often small |

### **Special MG2 Processes**

| MG2 Process | MG2 Rate Variable | P3 Equivalent | P3 Implementation Phase | Notes |
|-------------|-------------------|---------------|-------------------------|-------|
| Ice-ice collision | N/A in MG2 | `ice_self_collection` | Phase 6 | Aggregation |
| Rain-rain collision | N/A in MG2 | `rain_self_collection` | Phase 6 | Coalescence |
| Droplet self-collection | N/A in MG2 | `droplet_self_collection` | Phase 6 | Number only |

---

## **P3 Functions Detailed Mapping**

### **Phase 2: Rain Evaporation (Critical Path)**

#### **Function 1: `evaporate_rain`**

**Location**: `components/eam/src/physics/p3/eam/micro_p3.F90` (~line 2800)

**Purpose**: Calculate rain evaporation rate

**Inputs**:
- `qr` - rain mixing ratio (kg/kg)
- `nr` - rain number concentration (#/kg)
- `qv` - water vapor mixing ratio (kg/kg)
- `T` - temperature (K)
- `p` - pressure (Pa)
- `dt` - timestep (s)

**Outputs**:
- `qr_evap_tend` - rain evaporation tendency (kg/kg/s)
- `nr_evap_tend` - rain number evaporation tendency (#/kg/s)

**Isotope Requirements**:
- **Extract**: Evaporation rate for wtrc_p3_inter
- **Water Types**: iwtstrain → iwtvap
- **Fractionation**: Stewart model (kinetic + equilibrium)
- **Dependencies**: 
  - `stewart_isoevap()` - EAMv2-wiso water_tracers.F90:3678
  - `wtrc_equil_time()` - Calculate equilibration fraction
  - `wtrc_get_alpha()` - Get fractionation factors

**Interface Function**:
```fortran
! In wtrc_p3_inter_phase2():
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtstrain, iwtvap, iwtstrain, qr_evap_rate)

! Later in wtrc_apply_rates():
call stewart_isoevap(ispec, Rr0, rmass0, rmass, irmass0, &
                     vmass0, vmass, ivmass0, airtemp, airpres, &
                     pdel, dz, dtime, irmass, ivmass)
```

---

### **Phase 3: Vapor-Liquid Phase Changes**

#### **Function 2: `cloud_water_autoconversion`**

**Location**: `micro_p3.F90` (~line 1500)

**Purpose**: Convert cloud water to rain

**Inputs**:
- `qc` - cloud water mixing ratio (kg/kg)
- `nc` - cloud droplet number (#/kg)
- `rho` - air density (kg/m³)

**Outputs**:
- `qc2qr_autoconv_tend` - autoconversion tendency (kg/kg/s)
- `nc2nr_autoconv_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtliq → iwtstrain
- **Fractionation**: None (liquid → liquid)
- **Action**: Mass tracking only

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtstrain, iwtliq, qc2qr_autoconv_rate)
```

---

#### **Function 3: `cloud_rain_accretion`**

**Location**: `micro_p3.F90` (~line 1600)

**Purpose**: Cloud water collected by rain

**Inputs**:
- `qc`, `nc` - cloud water and number
- `qr`, `nr` - rain mixing ratio and number
- `rho` - air density

**Outputs**:
- `qc2qr_accret_tend` - accretion tendency (kg/kg/s)
- `nc2nr_accret_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtliq → iwtstrain
- **Fractionation**: None (liquid → liquid)
- **Action**: Mass tracking only

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtstrain, iwtliq, qc2qr_accret_rate)
```

---

#### **Function 4: Cloud Supersaturation Adjustment (Implicit)**

**Location**: Part of P3 main scheme

**Purpose**: Remove supersaturation by condensation

**Isotope Requirements**:
- **Water Types**: iwtvap → iwtliq
- **Fractionation**: Equilibrium + kinetic (Rayleigh distillation)
- **Dependencies**:
  - `wtrc_vap_distil()` - Rayleigh distillation
  - `wtrc_get_alpha()` - Equilibrium fractionation

**Interface Function**:
```fortran
! Condensation removes supersaturation
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtvap, iwtliq, iwtvap, qv_cond_rate)
```

---

### **Phase 4: Vapor-Ice Phase Changes**

#### **Function 5: `ice_nucleation`**

**Location**: `micro_p3.F90` (~line 1200)

**Purpose**: Ice nucleation from vapor

**Inputs**:
- `qv` - vapor mixing ratio
- `T` - temperature
- `Si` - ice supersaturation
- `nc` - cloud droplet number

**Outputs**:
- `qv2qi_nucleat_tend` - nucleation tendency (kg/kg/s)
- `ni_nucleat_tend` - ice number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtvap → iwtice (heterogeneous), iwtliq → iwtice (homogeneous)
- **Fractionation**: 
  - Equilibrium: `wiso_alpi()` for vapor → ice
  - Kinetic: `wiso_akci()` - supersaturation dependent
- **Critical**: Supersaturation affects kinetic fractionation strongly

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtvap, iwtice, iwtvap, qv2qi_nucleat_rate)
! Apply ice-vapor fractionation in wtrc_apply_rates
```

---

#### **Function 6: `ice_deposition_sublimation`**

**Location**: `micro_p3.F90` (~line 2200)

**Purpose**: Bidirectional vapor-ice exchange

**Inputs**:
- `qi` - ice mixing ratio
- `ni` - ice number
- `qv` - vapor mixing ratio
- `T` - temperature
- `Si` - ice supersaturation

**Outputs**:
- `qi_dep_subl_tend` - net deposition/sublimation (kg/kg/s)
- Sign: positive = deposition, negative = sublimation

**Isotope Requirements**:
- **Water Types**: 
  - Deposition: iwtvap → iwtice
  - Sublimation: iwtice → iwtvap
- **Fractionation**: Bidirectional
  - **Deposition**: Equilibrium (`wiso_alpi`) + kinetic (`wiso_akci`)
  - **Sublimation**: Equilibrium only (kinetic fractionation reversed)
- **Critical**: Must distinguish direction for correct fractionation sign

**Interface Function**:
```fortran
! Split into positive (deposition) and negative (sublimation)
if (qi_dep_subl_rate > 0.0) then
  call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                      iwtvap, iwtice, iwtvap, qi_dep_subl_rate)
else
  call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                      iwtice, iwtvap, iwtice, abs(qi_dep_subl_rate))
endif
```

---

#### **Function 7: `ice_supersat_conservation`**

**Location**: `micro_p3.F90` (~line 2300)

**Purpose**: Remove ice supersaturation

**Inputs**:
- `qv` - vapor mixing ratio
- `T` - temperature
- `Si` - ice supersaturation (> 100%)

**Outputs**:
- `qv2qi_dep_tend` - deposition to remove supersaturation (kg/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtvap → iwtice
- **Fractionation**: Same as deposition (equilibrium + kinetic)
- **Note**: Similar to `ice_deposition_sublimation` but only positive direction

**Interface Function**:
```fortran
call wtrc_add_rates(post_rates_grid, ncol, top_lev, &
                    iwtvap, iwtice, iwtvap, qv2qi_supersat_rate)
```

---

### **Phase 5: Freezing and Melting**

#### **Function 8: `cldliq_immersion_freezing`**

**Location**: `micro_p3.F90` (~line 1800)

**Purpose**: Cloud liquid freezing (immersion + contact)

**Inputs**:
- `qc` - cloud liquid mixing ratio
- `nc` - cloud droplet number
- `T` - temperature
- `INP` - ice nucleating particle concentration

**Outputs**:
- `qc2qi_immfrz_tend` - immersion freezing tendency (kg/kg/s)
- `nc2ni_immfrz_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtliq → iwtice
- **Fractionation**: Small kinetic effect (temperature dependent)
- **Note**: Fast process, minimal fractionation at cold temperatures

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtice, iwtliq, qc2qi_immfrz_rate)
```

---

#### **Function 9: `rain_immersion_freezing`**

**Location**: `micro_p3.F90` (~line 1900)

**Purpose**: Rain freezing to ice/snow

**Inputs**:
- `qr` - rain mixing ratio
- `nr` - rain number
- `T` - temperature

**Outputs**:
- `qr2qi_immfrz_tend` - rain freezing tendency (kg/kg/s)
- `nr2ni_immfrz_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtstrain → iwtstsnow (or iwtice)
- **Fractionation**: Minimal (fast freezing)
- **Note**: In MG2 this was `mnuccro` with special handling

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtstrain, iwtstsnow, iwtstrain, qr2qi_immfrz_rate, &
                    wtfri=wtfri_rain)  ! Frozen fraction
```

---

#### **Function 10: `homogeneous_freezing`**

**Location**: `micro_p3.F90` (~line 1750)

**Purpose**: Homogeneous freezing of cloud droplets (T < -37°C)

**Inputs**:
- `qc` - cloud liquid mixing ratio
- `nc` - cloud droplet number
- `T` - temperature (must be < T_homogfrz)

**Outputs**:
- `qc2qi_homfrz_tend` - homogeneous freezing tendency (kg/kg/s)
- `nc2ni_homfrz_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtliq → iwtice
- **Fractionation**: None (instantaneous, complete freezing)
- **Action**: Conservative mass transfer

**Interface Function**:
```fortran
call wtrc_add_rates(post_rates_grid, ncol, top_lev, &
                    iwtliq, iwtice, iwtliq, qc2qi_homfrz_rate)
```

---

#### **Function 11: `ice_melting`**

**Location**: `micro_p3.F90` (~line 2500)

**Purpose**: Ice melting to liquid

**Inputs**:
- `qi` - ice mixing ratio
- `ni` - ice number
- `T` - temperature (must be > T_zerodegc)

**Outputs**:
- `qi2qr_melt_tend` - melting tendency (kg/kg/s)
- `ni2nr_melt_tend` - number tendency (#/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtice → iwtliq (or iwtstrain)
- **Fractionation**: None (equilibrium at 0°C)
- **Action**: Conservative mass transfer

**Interface Function**:
```fortran
call wtrc_add_rates(post_rates_grid, ncol, top_lev, &
                    iwtice, iwtliq, iwtice, qi2qr_melt_rate)
```

---

#### **Function 12: `ice_complete_melting`**

**Location**: `micro_p3.F90` (~line 2550)

**Purpose**: Complete melting of remaining ice

**Inputs**:
- `qi` - ice mixing ratio
- `T` - temperature

**Outputs**:
- Complete transfer of ice to liquid

**Isotope Requirements**:
- **Water Types**: iwtice → iwtliq
- **Fractionation**: None
- **Action**: Conservative complete transfer

---

### **Phase 6: Collection and Riming**

#### **Function 13: `ice_cldliq_collection`**

**Location**: `micro_p3.F90` (~line 2000)

**Purpose**: Ice collection of cloud liquid (riming)

**Inputs**:
- `qi`, `ni` - ice mixing ratio and number
- `qc`, `nc` - cloud liquid and number
- `T` - temperature
- `rho` - air density

**Outputs**:
- `qc2qi_riming_tend` - riming tendency (kg/kg/s)
- `qi_rime_mass_tend` - rime mass accumulation

**Isotope Requirements**:
- **Water Types**: iwtliq → iwtice (rime mass)
- **Fractionation**: Kinetic effects during impact/freezing
- **Note**: Rime mass tracking in P3 (important for isotopes)

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtice, iwtliq, qc2qi_riming_rate)
```

---

#### **Function 14: `ice_cldliq_wet_growth`**

**Location**: `micro_p3.F90` (~line 2100)

**Purpose**: Wet growth regime with shedding

**Inputs**:
- `qi`, `ni` - ice
- `qc`, `nc` - cloud liquid
- `T` - temperature (near 0°C)
- Determines fraction frozen vs shed

**Outputs**:
- `qc2qi_wetgrowth_tend` - wet growth tendency
- `qi2qr_shed_tend` - shedding tendency (liquid shed as rain)

**Isotope Requirements**:
- **Water Types**: 
  - Freezing: iwtliq → iwtice
  - Shedding: iwtliq → iwtstrain
- **Fractionation**: Partial equilibration
- **Complexity**: HIGH - need to track both pathways

**Interface Function**:
```fortran
! Frozen fraction
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtice, iwtliq, qc2qi_wetgrowth_frz)
! Shed fraction
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtliq, iwtstrain, iwtliq, qc2qr_wetgrowth_shed)
```

---

#### **Function 15: `ice_rain_collection`**

**Location**: `micro_p3.F90` (~line 2150)

**Purpose**: Ice collection of rain

**Inputs**:
- `qi`, `ni` - ice
- `qr`, `nr` - rain
- `T` - temperature

**Outputs**:
- `qr2qi_collect_tend` - collection tendency (kg/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtstrain → iwtice (or iwtstsnow)
- **Fractionation**: Minimal (liquid already fractionated)
- **Action**: Mass transfer with possible partial freezing

**Interface Function**:
```fortran
call wtrc_add_rates(pre_rates_grid, ncol, top_lev, &
                    iwtstrain, iwtice, iwtstrain, qr2qi_collect_rate)
```

---

#### **Function 16-18: Self-Collection Processes**

**Functions**:
- `ice_self_collection` (micro_p3.F90 ~line 2250)
- `rain_self_collection` (micro_p3.F90 ~line 1650)
- `droplet_self_collection` (micro_p3.F90 ~line 1450)

**Purpose**: Aggregation/coalescence within same phase

**Isotope Requirements**:
- **Fractionation**: None (same phase)
- **Action**: Number concentration changes only, mass conserved
- **Note**: No isotope interface needed (conservative)

---

### **Phase 7: Sedimentation**

#### **Function 19: `rain_sedimentation`**

**Location**: `micro_p3.F90` (~line 3200)

**Purpose**: Rain vertical transport with evaporation

**Inputs**:
- `qr`, `nr` - rain mixing ratio and number
- `qv` - vapor mixing ratio
- `T`, `p` - thermodynamic state
- `dz` - layer thickness
- `dt` - timestep

**Outputs**:
- `qr_sed_tend` - sedimentation tendency (kg/kg/s)
- `qr_evap_sed_tend` - evaporation during fall (kg/kg/s)
- `precip_flux` - surface precipitation flux (kg/m²/s)

**Isotope Requirements**:
- **Water Types**: 
  - Transport: iwtstrain (vertical redistribution)
  - Evaporation: iwtstrain → iwtvap (Stewart model)
- **Fractionation**: Stewart model applied layer-by-layer
- **Complexity**: HIGH - iterative, CFL-limited
- **Dependencies**: 
  - `wtrc_sediment()` - Full sedimentation handler
  - `stewart_isoevap()` - Per-layer evaporation
  - `wtrc_equil_time()` - Equilibration fraction

**Interface Function**:
```fortran
! Call dedicated sedimentation routine
call wtrc_sediment(state, ptend, pbuf, precr, preci, top_lev, dtime, &
                   qr_sed_flux, qr_evap_rate, liqcldf, icecldf, &
                   fr, fs, T, p, qv, zi)
```

**Critical Implementation**:
- CFL sub-stepping for stability
- Layer-by-layer equilibration
- Raindrop size evolution during fall
- Clear vs cloudy sky separation

---

#### **Function 20: `ice_sedimentation`**

**Location**: `micro_p3.F90` (~line 3400)

**Purpose**: Ice/snow vertical transport with sublimation

**Inputs**:
- `qi`, `ni` - ice mixing ratio and number
- `qv` - vapor mixing ratio
- `T`, `p` - thermodynamic state
- `dz`, `dt` - layer thickness and timestep

**Outputs**:
- `qi_sed_tend` - sedimentation tendency (kg/kg/s)
- `qi_subl_sed_tend` - sublimation during fall (kg/kg/s)
- `snow_flux` - surface snow flux (kg/m²/s)

**Isotope Requirements**:
- **Water Types**: 
  - Transport: iwtstsnow (or iwtice)
  - Sublimation: iwtstsnow → iwtvap
- **Fractionation**: Ice-vapor equilibrium during sublimation
- **Dependencies**: Similar to rain sedimentation

**Interface Function**:
```fortran
! Part of wtrc_sediment call
! Snow/ice handled separately from rain
```

---

#### **Function 21: `cloud_sedimentation`**

**Location**: `micro_p3.F90` (~line 3100)

**Purpose**: Cloud droplet sedimentation (usually small)

**Inputs**:
- `qc`, `nc` - cloud liquid
- `dz`, `dt` - layer thickness and timestep

**Outputs**:
- `qc_sed_tend` - cloud sedimentation tendency (kg/kg/s)

**Isotope Requirements**:
- **Water Types**: iwtliq (vertical redistribution)
- **Fractionation**: Minimal (stays in cloud)
- **Action**: Simple mass transport

---

### **Phase 8: Convection Interface**

#### **Function 22: Zhang-McFarlane Deep Convection**

**Location**: `components/eam/src/physics/cam/zm_conv_intr.F90`

**Purpose**: Deep convective transport and precipitation

**Isotope Function**: `wtrc_q1q2_pjr()` (from EAMv2-wiso water_tracers.F90:4960)

**Processes**:
1. Updraft condensation (with fractionation)
2. Downdraft evaporation (with fractionation)
3. Convective precipitation
4. Detrainment to stratiform

**Water Types**:
- Updraft: iwtvap → iwtliq (condensation)
- Downdraft: iwtcvrain → iwtvap (evaporation)
- Detrainment: iwtliq → environment

**Fractionation**:
- Condensation: Rayleigh distillation
- Evaporation: Partial equilibration

---

#### **Function 23: Shallow Convection**

**Location**: `components/eam/src/physics/cam/uwshcu.F90`

**Purpose**: Shallow cumulus transport

**Isotope Function**: `wtrc_shallow()` (from EAMv2-wiso water_tracers.F90:6234)

**Processes**:
1. Shallow convective mixing
2. Condensation in updrafts
3. Minimal precipitation

**Water Types**: Similar to deep convection but smaller magnitude

---

### **Phase 9: Conservation and Diagnostics**

#### **Function 24: Mass Conservation Checks**

**Functions to Port**:
- `wtrc_check_h2o()` - Check total water mass conservation
- `wtrc_adjust_h2o()` - Fix numerical issues
- `wtrc_mass_fixer()` - Reset tracers if needed

**Location**: EAMv2-wiso water_tracers.F90:1388, 1513, 6677

---

#### **Function 25: Diagnostic Output**

**Functions to Port**:
- `wtrc_output_precip()` - Precipitation isotope output
- `wtrc_setup_diag()` - Diagnostic field registration

**Output Fields**:
- δD, δ18O, δ17O for all phases
- Precipitation isotope ratios
- Column-integrated budgets

---

## **Water Type Transition Matrix**

### **Process Rate Matrix Structure**

In the process rate framework, transitions are tracked as:

```
rates(i,k, source_type, destination_type, species_type)
```

Where:
- `source_type` = water type being consumed
- `destination_type` = water type being produced  
- `species_type` = which tracer (H2O, HDO, H218O, etc.)

### **Transition Matrix for P3 Processes**

| Source → Destination | Process | P3 Function | Fractionation | Phase |
|---------------------|---------|-------------|---------------|-------|
| **Vapor → Liquid** |
| iwtvap → iwtliq | Condensation | Supersaturation adjustment | Equilibrium + kinetic | 3 |
| iwtvap → iwtliq | Homogeneous nucleation | Part of `ice_nucleation` | Equilibrium | 4 |
| **Liquid → Vapor** |
| iwtliq → iwtvap | Cloud evaporation | Part of saturation adjustment | Equilibrium + kinetic | 3 |
| iwtstrain → iwtvap | Rain evaporation | `evaporate_rain` | **Stewart model** | 2 |
| **Vapor → Ice** |
| iwtvap → iwtice | Ice nucleation | `ice_nucleation` | Equilibrium + kinetic | 4 |
| iwtvap → iwtice | Deposition | `ice_deposition_sublimation` | Equilibrium + kinetic | 4 |
| iwtvap → iwtice | Supersat removal | `ice_supersat_conservation` | Equilibrium + kinetic | 4 |
| **Ice → Vapor** |
| iwtice → iwtvap | Sublimation | `ice_deposition_sublimation` | Equilibrium | 4 |
| iwtstsnow → iwtvap | Snow sublimation | `ice_sedimentation` | Equilibrium | 7 |
| **Liquid → Ice** |
| iwtliq → iwtice | Immersion freezing | `cldliq_immersion_freezing` | Kinetic (small) | 5 |
| iwtliq → iwtice | Homogeneous freezing | `homogeneous_freezing` | None | 5 |
| iwtliq → iwtice | Riming | `ice_cldliq_collection` | Kinetic | 6 |
| iwtliq → iwtice | Wet growth (frozen) | `ice_cldliq_wet_growth` | Kinetic | 6 |
| iwtstrain → iwtice | Rain freezing | `rain_immersion_freezing` | Minimal | 5 |
| iwtstrain → iwtice | Rain collection | `ice_rain_collection` | Minimal | 6 |
| **Ice → Liquid** |
| iwtice → iwtliq | Ice melting | `ice_melting` | None | 5 |
| iwtice → iwtliq | Complete melting | `ice_complete_melting` | None | 5 |
| iwtstsnow → iwtstrain | Snow melting | Part of `ice_melting` | None | 5 |
| **Liquid → Liquid** |
| iwtliq → iwtstrain | Autoconversion | `cloud_water_autoconversion` | None | 3 |
| iwtliq → iwtstrain | Accretion | `cloud_rain_accretion` | None | 3 |
| iwtliq → iwtstrain | Wet growth (shed) | `ice_cldliq_wet_growth` | None | 6 |
| **Ice → Ice** |
| iwtice → iwtice | Self-collection | `ice_self_collection` | None | 6 |
| **Rain → Rain** |
| iwtstrain → iwtstrain | Self-collection | `rain_self_collection` | None | 6 |

### **Bergeron Process (Special Case)**

The Bergeron process transfers vapor from liquid to ice (vapor supersaturated over liquid, subsaturated over ice):

| Source → Destination | Actual Path | Fractionation |
|---------------------|-------------|---------------|
| iwtliq → iwtice | iwtliq → iwtvap → iwtice | Two-step: liq-vap equilibrium, then vap-ice equilibrium + kinetic |

In MG2 this was handled with special `bergo` and `bergso` rates. In P3 it's implicit in the supersaturation adjustment.

---

## **Fractionation Factor Requirements**

### **Equilibrium Fractionation Factors**

| Transition | Function | Reference | Temperature Dependence | Phase |
|-----------|----------|-----------|----------------------|-------|
| Liquid-Vapor | `wiso_alpl()` | Horita & Wesolowski (1994) | α = exp(A/T² + B/T + C) | 2,3 |
| Ice-Vapor | `wiso_alpi()` | Merlivat & Nief (1967), Majoube (1971) | α = exp(A/T² + B/T + C) | 4 |

**Implementation**: `water_isotopes.F90` functions

**Typical Values** (at 0°C):
- HDO liquid-vapor: α ≈ 1.084
- H218O liquid-vapor: α ≈ 1.0098
- HDO ice-vapor: α ≈ 1.161
- H218O ice-vapor: α ≈ 1.0135

### **Kinetic Fractionation Factors**

| Process | Function | Formula | Dependencies | Phase |
|---------|----------|---------|--------------|-------|
| Liquid evaporation | `wiso_akel()` | α_k = (D_H2O/D_iso)^n | Diffusivity ratio, turbulence | 2,3 |
| Ice condensation | `wiso_akci()` | α_k = f(Si) | Supersaturation | 4 |

**Diffusivity Ratios** (difrm):
- HDO: 0.9755
- H218O: 0.9723
- H217O: 0.9847

**Exponent** (dkfac):
- Typical value: n ≈ 0.5 (partially turbulent)
- Range: 0 (fully turbulent) to 1 (molecular diffusion)

### **Combined Fractionation**

For most processes, total fractionation is:

```
α_total = α_equilibrium × α_kinetic
```

**Special Cases**:
1. **Stewart Rain Evaporation**: Complex formula involving humidity, droplet size, and fall time
2. **Wet Growth**: Partial equilibration based on freezing efficiency
3. **Homogeneous Freezing**: α = 1 (no fractionation)

---

## **EAMv2-wiso Function Dependencies**

### **Core Functions Required for All Phases**

| Function | Location | Purpose | Used In Phase |
|----------|----------|---------|---------------|
| `wiso_alpl()` | water_isotopes.F90:145 | Liquid-vapor equilibrium | 2,3,7 |
| `wiso_alpi()` | water_isotopes.F90:235 | Ice-vapor equilibrium | 4,7 |
| `wiso_akel()` | water_isotopes.F90:325 | Kinetic fractionation (liquid) | 2,3,7 |
| `wiso_akci()` | water_isotopes.F90:385 | Kinetic fractionation (ice) | 4,7 |
| `wtrc_init_rates()` | water_tracers.F90:1086 | Initialize rate arrays | 2+ |
| `wtrc_add_rates()` | water_tracers.F90:1138 | Add process rate | 2+ |
| `wtrc_get_alpha()` | water_tracers.F90:1298 | Get fractionation factor | 2+ |
| `wtrc_ratio()` | water_tracers.F90:1262 | Calculate isotope ratio | 2+ |

### **Phase-Specific Function Dependencies**

#### **Phase 2: Rain Evaporation**

| Function | Location | Purpose |
|----------|----------|---------|
| `stewart_isoevap()` | water_tracers.F90:3678 | Stewart rain evaporation model |
| `wtrc_equil_time()` | water_tracers.F90:6580 | Calculate equilibration fraction |
| `wtrc_liqvap_equil()` | water_tracers.F90:6381 | Iterative equilibration solver |
| `wtrc_apply_rates()` | water_tracers.F90:2068 | Apply rates with fractionation |
| `wtrc_eff_sat()` | water_tracers.F90:6885 | Effective saturation for kinetic effects |

#### **Phase 3: Condensation**

| Function | Location | Purpose |
|----------|----------|---------|
| `wtrc_vap_distil()` | water_tracers.F90:2776 | Rayleigh distillation |
| `wtrc_dicm()` | water_tracers.F90:5715 | Iterative equilibration (alternative) |

#### **Phase 7: Sedimentation**

| Function | Location | Purpose |
|----------|----------|---------|
| `wtrc_sediment()` | water_tracers.F90:3315 | Full sedimentation handler |
| (includes calls to stewart_isoevap, wtrc_liqvap_equil, wtrc_equil_time) |

#### **Phase 8: Convection**

| Function | Location | Purpose |
|----------|----------|---------|
| `wtrc_q1q2_pjr()` | water_tracers.F90:4960 | ZM deep convection isotopes |
| `wtrc_shallow()` | water_tracers.F90:6234 | Shallow convection isotopes |

#### **Phase 9: Diagnostics**

| Function | Location | Purpose |
|----------|----------|---------|
| `wtrc_check_h2o()` | water_tracers.F90:1388 | Mass conservation check |
| `wtrc_adjust_h2o()` | water_tracers.F90:1513 | Fix negative values |
| `wtrc_mass_fixer()` | water_tracers.F90:6677 | Reset to bulk water |
| `wtrc_output_precip()` | water_tracers.F90:4161 | Precipitation output |
| `wtrc_setup_diag()` | water_tracers.F90:1356 | Diagnostic setup |

---

## **Implementation Checklist by Phase**

### **Phase 1: Core Infrastructure**
- [ ] Port water_isotopes.F90
- [ ] Port water_types.F90  
- [ ] Port water_tracer_vars.F90
- [ ] Port water_tracers.F90 (core functions only)
- [ ] Extend physics_state for isotope tracers
- [ ] Register isotope constituents
- [ ] Add namelist parameters

### **Phase 2: Rain Evaporation**
- [ ] Modify evaporate_rain() to output rate
- [ ] Create wtrc_p3_inter_phase2()
- [ ] Port stewart_isoevap()
- [ ] Port wtrc_equil_time()
- [ ] Port wtrc_liqvap_equil()
- [ ] Port wtrc_apply_rates() (simplified)
- [ ] Test rain evaporation fractionation

### **Phase 3: Vapor-Liquid**
- [ ] Modify cloud_water_autoconversion()
- [ ] Modify cloud_rain_accretion()
- [ ] Handle saturation adjustment
- [ ] Port wtrc_vap_distil()
- [ ] Extend wtrc_apply_rates()
- [ ] Test condensation fractionation

### **Phase 4: Vapor-Ice**
- [ ] Modify ice_nucleation()
- [ ] Modify ice_deposition_sublimation()
- [ ] Modify ice_supersat_conservation()
- [ ] Implement ice-vapor fractionation
- [ ] Handle supersaturation-dependent kinetics
- [ ] Test ice fractionation

### **Phase 5: Freezing/Melting**
- [ ] Modify all freezing functions
- [ ] Modify all melting functions
- [ ] Implement freezing fractionation
- [ ] Test phase transitions

### **Phase 6: Collection**
- [ ] Modify ice_cldliq_collection()
- [ ] Modify ice_cldliq_wet_growth()
- [ ] Modify ice_rain_collection()
- [ ] Implement riming fractionation
- [ ] Test collection processes

### **Phase 7: Sedimentation**
- [ ] Modify rain_sedimentation()
- [ ] Modify ice_sedimentation()
- [ ] Port wtrc_sediment()
- [ ] Implement layer-by-layer fractionation
- [ ] Test precipitation isotopes

### **Phase 8: Convection**
- [ ] Port wtrc_q1q2_pjr()
- [ ] Port wtrc_shallow()
- [ ] Integrate with ZM and UW schemes
- [ ] Test convective isotopes

### **Phase 9: Conservation**
- [ ] Port conservation functions
- [ ] Port diagnostic functions
- [ ] Add history output
- [ ] Test global conservation

---

**End of Function Mapping Document**
