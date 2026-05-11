# Implementation Roadmap for Remaining High-Priority Functions

**Date:** February 17, 2026  
**Status:** Ready for implementation with EAMxx infrastructure

---

## Overview

The remaining 8 high-priority functions (~2,750 lines) require EAMxx infrastructure that is not yet available in the skeleton code. This document provides detailed implementation guidance for each function.

## Required EAMxx Infrastructure

Before implementing the remaining functions, the following EAMxx components are needed:

### 1. Physics State Structure
```cpp
struct PhysicsState {
  int ncol;                          // Number of columns
  int lchnk;                         // Chunk ID
  View2D<Real> q;                    // Constituent mixing ratios [ncol, nlev, ncnst]
  View2D<Real> t;                    // Temperature [ncol, nlev]
  View2D<Real> pmid;                 // Mid-level pressure [ncol, nlev]
  View2D<Real> pdel;                 // Layer thickness (pressure) [ncol, nlev]
  View2D<Real> zi;                   // Interface heights [ncol, nlev+1]
  // ... other state variables
};
```

### 2. Physics Tendency Structure
```cpp
struct PhysicsTendency {
  std::string name;                  // Process name
  View3D<Real> q;                    // Constituent tendencies [ncol, nlev, ncnst]
  View1D<bool> lq;                   // Flags for which constituents have tendencies
  // ... other tendency variables
};
```

### 3. Physics Buffer (or Field Manager)
- Kokkos::View-based replacement for Fortran physics_buffer
- Used for precipitation output and temporary fields
- Needs get_field/set_field interface

### 4. Saturation Functions
- `qsat_water(temp, press, es, qs)` - Saturation specific humidity over liquid
- `qsat_ice(temp, press, esi, qsi)` - Saturation specific humidity over ice

---

## Function 1: wtrc_mg_inter (252 lines)

**Purpose:** Interface between MG2 microphysics and water tracer system

**Complexity:** MEDIUM - Straightforward rate mapping

**Implementation Strategy:**

### Algorithm Outline
```
1. Get number of columns from state
2. Initialize pre_rates and post_rates arrays (call wtrc_init_rates)
3. Split deposition/sublimation rates by sign (pcmei/ncmei)
4. Split melting rates by sign (pmelts/nmelts)
5. Map MG2 process rates to pre_rates matrix:
   - Vapor processes: deposition, sublimation, rain/snow evaporation
   - Liquid processes: freezing, accretion, autoconversion, Bergeron
   - Ice processes: accretion to snow, autoconversion
   - Rain processes: accretion to snow, freezing
6. Map MG2 process rates to post_rates matrix:
   - Condensation/deposition (supersaturation removal)
   - Freezing, melting
   - Precipitation phase changes
7. Call wtrc_apply_rates with prepared matrices
8. Return precipitation output (precr, preci)
```

### Key MG2 Process Rate Mappings

#### Pre-Sedimentation Rates
| MG2 Variable | Process | Source → Dest | Rate Type |
|--------------|---------|---------------|-----------|
| pcmei | Deposition | Vapor → Ice | Vapor |
| ncmei | Sublimation | Ice → Vapor | Ice |
| preo | Rain evap | Rain → Vapor | Rain |
| prdso | Snow sublim | Snow → Vapor | Snow |
| mnuccco+mnuccto+msacwio | Freezing | Liquid → Ice | Liquid |
| prao+prco | Autoconv/accr | Liquid → Rain | Liquid |
| psacwso | Accretion | Liquid → Snow | Liquid |
| bergo | Bergeron (liq) | Liquid → Liquid | Liquid |
| bergso | Bergeron (ice) | Ice → Ice | Ice |
| praio+prcio | Ice to snow | Ice → Snow | Ice |
| pracso | Rain to snow | Rain → Snow | Rain |
| mnuccro | Rain freezing | Rain → Snow/Ice | Rain |

#### Post-Sedimentation Rates
| MG2 Variable | Process | Source → Dest | Rate Type |
|--------------|---------|---------------|-----------|
| qcreso | Condensation | Vapor → Liquid | Vapor |
| qireso | Deposition | Vapor → Ice | Vapor |
| homoo | Homog freezing | Liquid → Ice | Liquid |
| melto | Melting | Ice → Liquid | Ice |
| frzrpst | Rain freezing | Rain → Snow | Rain |
| meltspst | Snow melting | Snow → Rain | Snow |

### Pseudocode
```cpp
void wtrc_mg_inter(PhysicsState& state, PhysicsTendency& ptend, PhysicsBuffer& pbuf,
                   View2D<Real>& precr, View2D<Real>& preci, int top_lev, Real dtime,
                   const View2D<Real>& alst_mic, const View2D<Real>& aist_mic,
                   // ... all MG2 process rate inputs ...
                   const View2D<Real>& wtfc, const View2D<Real>& wtfi,
                   const View2D<Real>& wtfr, const View2D<Real>& wtfs,
                   const View2D<Real>& wtprelat, const View2D<Real>& wtsedlat,
                   const View2D<Real>& wtpostlat)
{
  int ncol = state.ncol;
  
  // Allocate rate matrices
  View5D<Real> pre_rates("pre_rates", ncol, PVER, PWTYPE, PWTYPE, PWTYPE);
  View5D<Real> post_rates("post_rates", ncol, PVER, PWTYPE, PWTYPE, PWTYPE);
  
  // Initialize to zero
  wtrc_init_rates(top_lev, pre_rates);
  wtrc_init_rates(top_lev, post_rates);
  
  // Split deposition/sublimation and melting by sign
  View2D<Real> pcmei("pcmei", ncol, PVER);
  View2D<Real> ncmei("ncmei", ncol, PVER);
  View2D<Real> pmelts("pmelts", ncol, PVER);
  View2D<Real> nmelts("nmelts", ncol, PVER);
  
  Kokkos::parallel_for("split_rates", MDRangePolicy<2>({0,top_lev},{ncol,PVER}),
    KOKKOS_LAMBDA(int i, int k) {
      pcmei(i,k) = (cmeiout(i,k) >= 0.0) ? cmeiout(i,k) : 0.0;
      ncmei(i,k) = (cmeiout(i,k) < 0.0) ? cmeiout(i,k) : 0.0;
      pmelts(i,k) = (meltso(i,k) >= 0.0) ? meltso(i,k) : 0.0;
      nmelts(i,k) = (meltso(i,k) < 0.0) ? meltso(i,k) : 0.0;
    });
  
  // Add pre-sedimentation rates
  using namespace water_types;
  int vap = static_cast<int>(WaterType::Vapor);
  int liq = static_cast<int>(WaterType::CloudLiquid);
  int ice = static_cast<int>(WaterType::CloudIce);
  int rain = static_cast<int>(WaterType::StratRain);
  int snow = static_cast<int>(WaterType::StratSnow);
  
  wtrc_add_rates(pre_rates, ncol, top_lev, vap, ice, vap, pcmei);   // deposition
  wtrc_add_rates(pre_rates, ncol, top_lev, vap, ice, ice, ncmei);   // sublimation
  wtrc_add_rates(pre_rates, ncol, top_lev, rain, vap, rain, preo);  // rain evap
  wtrc_add_rates(pre_rates, ncol, top_lev, snow, vap, snow, prdso); // snow sublim
  
  // Liquid processes
  auto liq_freeze = mnuccco + mnuccto + msacwio;  // Need element-wise add
  wtrc_add_rates(pre_rates, ncol, top_lev, liq, ice, liq, liq_freeze);
  
  auto liq_to_rain = prao + prco;
  wtrc_add_rates(pre_rates, ncol, top_lev, liq, rain, liq, liq_to_rain);
  
  wtrc_add_rates(pre_rates, ncol, top_lev, liq, snow, liq, psacwso);
  wtrc_add_rates(pre_rates, ncol, top_lev, liq, liq, liq, bergo);  // Bergeron
  wtrc_add_rates(pre_rates, ncol, top_lev, ice, ice, ice, bergso);  // Bergeron
  
  // Ice processes
  auto ice_to_snow = praio + prcio;
  wtrc_add_rates(pre_rates, ncol, top_lev, ice, snow, ice, ice_to_snow);
  
  // Rain processes
  wtrc_add_rates(pre_rates, ncol, top_lev, rain, snow, rain, pracso);
  wtrc_add_rates(pre_rates, ncol, top_lev, rain, snow, rain, mnuccro, true, &wtfri_pre);
  
  // Add post-sedimentation rates
  wtrc_add_rates(post_rates, ncol, top_lev, vap, liq, vap, qcreso);  // condensation
  wtrc_add_rates(post_rates, ncol, top_lev, vap, ice, vap, qireso);  // deposition
  wtrc_add_rates(post_rates, ncol, top_lev, liq, ice, liq, homoo);   // freezing
  wtrc_add_rates(post_rates, ncol, top_lev, ice, liq, ice, melto);   // melting
  wtrc_add_rates(post_rates, ncol, top_lev, rain, snow, rain, frzrpst, true, &wtfri_post);
  wtrc_add_rates(post_rates, ncol, top_lev, snow, rain, snow, meltspst);
  
  // Call wtrc_apply_rates
  wtrc_apply_rates(state, ptend, pbuf, top_lev, dtime, true, precr, preci,
                   pre_rates, post_rates, alst_mic, aist_mic,
                   wtfc, wtfi, wtfr, wtfs, wtprelat, wtsedlat, wtpostlat,
                   frzro, meltso, wtfri_pre);
}
```

---

## Function 2: wtrc_sediment_mg1 (286 lines)

**Purpose:** Calculate sedimentation of cloud liquid and ice for water tracers (MG1)

**Complexity:** HIGH - CFL sub-stepping, cloud fractions, equilibration

**Implementation Strategy:**

### Algorithm Outline
```
1. Save initial temperature
2. Calculate cloud fractions (liquid, ice) with minimum threshold
3. Initialize precipitation arrays
4. Loop over columns:
   a. Get initial fall velocities (fc, fi)
   b. Calculate CFL-stable sub-step: nstep = max(rgvm*dtime/pdel+1, 1)
   c. Loop over sub-steps:
      i. Calculate sedimentation fluxes: falout = fall_velocity * q
      ii. At model top: apply flux divergence
      iii. Below model top:
          - Calculate cloud fraction ratio
          - Apply instant evap for flux entering clear sky
          - Update states (vapor, liquid, ice)
          - Equilibrate vapor/liquid for isotopes
          - Update temperature from latent heating
      iv. Reset fall velocities for next sub-step
   d. Accumulate surface precipitation
5. Output to physics buffer
```

### Key Physics

#### CFL Stability
```
Maximum fall velocity: rgvm = max(fi, fc)
Sub-steps needed: nstep = ceiling(rgvm * dt / dz)
Ensures Courant number <= 1 for stable sedimentation
```

#### Cloud Fraction Treatment
```
dum = cldf(k) / cldf(k-1)  // Fraction ratio
dum = min(dum, 1.0)

Flux entering from above:
- Fraction dum enters cloudy part → stays condensed
- Fraction (1-dum) enters clear part → instant evaporation
```

#### Isotopic Equilibration
For each sub-step, equilibrate newly evaporated condensate with vapor:
```
alpha = wtrc_get_alpha(qvap, temp, ispec, vapor, liquid, kinetic=false, RH=1.0)
call wtrc_liqvap_equil(alpha, feq=1.0, qvap_tot, qliq_tot, qvap_iso, qliq_iso, dliqiso)
```

### Pseudocode
```cpp
void wtrc_sediment_mg1(int niter, int ncol, int lchnk, int top_lev, Real dtime,
                       const View2D<Real>& fc, const View2D<Real>& fi,
                       const View2D<Real>& liqcldf, const View2D<Real>& icecldf,
                       const View2D<Real>& pdel, PhysicsBuffer& pbuf,
                       View2D<Real>& tloc, View3D<Real>& qloc)
{
  constexpr Real mincld = 0.0001;
  constexpr Real latvap = 2.5e6;   // Latent heat of vaporization (J/kg)
  constexpr Real latice = 3.34e5;  // Latent heat of fusion (J/kg)
  constexpr Real cpair = 1004.0;   // Specific heat of air (J/kg/K)
  constexpr Real gravit = 9.80616; // Gravity (m/s^2)
  
  // Save initial temperature
  View2D<Real> tloc0("tloc0", ncol, PVER);
  Kokkos::deep_copy(tloc0, tloc);
  
  // Calculate cloud amounts with minimum
  View2D<Real> lcldm("lcldm", ncol, PVER);
  View2D<Real> icldm("icldm", ncol, PVER);
  Kokkos::parallel_for("calc_cldf", MDRangePolicy<2>({0,top_lev},{ncol,PVER}),
    KOKKOS_LAMBDA(int i, int k) {
      lcldm(i,k) = max(liqcldf(i,k), mincld);
      icldm(i,k) = max(icecldf(i,k), mincld);
    });
  
  // Initialize precipitation
  View2D<Real> precr("precr", ncol, wtrc_nwset);
  View2D<Real> preci("preci", ncol, wtrc_nwset);
  Kokkos::deep_copy(precr, 0.0);
  Kokkos::deep_copy(preci, 0.0);
  
  // Loop over columns (sequential for CFL calculation)
  for (int i = 0; i < ncol; ++i) {
    // Get fall velocities
    View1D<Real> fcol("fcol", PVER);
    View1D<Real> fice("fice", PVER);
    Kokkos::parallel_for("get_fall_vel", RangePolicy(top_lev, PVER),
      KOKKOS_LAMBDA(int k) {
        fcol(k) = fc(i,k);
        fice(k) = fi(i,k);
      });
    
    // Calculate CFL-stable sub-step
    Real rgvm = 0.0;
    Kokkos::parallel_reduce("max_fall_vel", RangePolicy(top_lev, PVER),
      KOKKOS_LAMBDA(int k, Real& max_val) {
        Real vel = max(fice(k), fcol(k));
        if (vel > max_val) max_val = vel;
      }, Kokkos::Max<Real>(rgvm));
    
    int nstep = 1;
    Kokkos::parallel_reduce("calc_nstep", RangePolicy(top_lev, PVER),
      KOKKOS_LAMBDA(int k, int& max_step) {
        int step = static_cast<int>(rgvm * dtime / pdel(i,k) + 1.0);
        if (step > max_step) max_step = step;
      }, Kokkos::Max<int>(nstep));
    
    // Loop over sub-steps
    for (int n = 0; n < nstep; ++n) {
      // Calculate sedimentation fluxes
      View2D<Real> falouti("falouti", PVER, wtrc_nwset);
      View2D<Real> faloutc("faloutc", PVER, wtrc_nwset);
      
      Kokkos::parallel_for("calc_fluxes", MDRangePolicy<2>({top_lev,0},{PVER,wtrc_nwset}),
        KOKKOS_LAMBDA(int k, int m) {
          int idx_ice = wtrc_iatype(m, iwtice);
          int idx_liq = wtrc_iatype(m, iwtliq);
          falouti(k,m) = fice(k) * qloc(i,k,idx_ice);
          faloutc(k,m) = fcol(k) * qloc(i,k,idx_liq);
        });
      
      // Apply at model top
      int k = top_lev;
      Kokkos::parallel_for("apply_top", RangePolicy(0, wtrc_nwset),
        KOKKOS_LAMBDA(int m) {
          Real faltndi = falouti(k,m) / pdel(i,k);
          Real faltndc = faloutc(k,m) / pdel(i,k);
          int idx_ice = wtrc_iatype(m, iwtice);
          int idx_liq = wtrc_iatype(m, iwtliq);
          qloc(i,k,idx_ice) -= faltndi * dtime / nstep;
          qloc(i,k,idx_liq) -= faltndc * dtime / nstep;
        });
      
      // Loop over rest of column
      for (int k = top_lev + 1; k < PVER; ++k) {
        // Calculate cloud fraction ratios
        Real dum = min(lcldm(i,k) / lcldm(i,k-1), 1.0);
        Real dum1 = min(icldm(i,k) / icldm(i,k-1), 1.0);
        
        // Apply flux divergence with evaporation
        for (int m = 0; m < wtrc_nwset; ++m) {
          Real faltndqie = (falouti(k,m) - falouti(k-1,m)) / pdel(i,k);
          Real faltndi = (falouti(k,m) - dum1 * falouti(k-1,m)) / pdel(i,k);
          Real faltndqce = (faloutc(k,m) - faloutc(k-1,m)) / pdel(i,k);
          Real faltndc = (faloutc(k,m) - dum * faloutc(k-1,m)) / pdel(i,k);
          
          int idx_vap = wtrc_iatype(m, iwtvap);
          int idx_ice = wtrc_iatype(m, iwtice);
          int idx_liq = wtrc_iatype(m, iwtliq);
          
          // Add evap/sublimation to vapor
          qloc(i,k,idx_vap) -= (faltndqie - faltndi) * dtime / nstep;
          qloc(i,k,idx_vap) -= (faltndqce - faltndc) * dtime / nstep;
          
          // Update condensate
          qloc(i,k,idx_ice) -= faltndi * dtime / nstep;
          qloc(i,k,idx_liq) -= faltndc * dtime / nstep;
          
          // Equilibrate if isotope
          if (wisotope && m != 0) {
            Real alpha = wtrc_get_alpha(qloc(i,k,wtrc_iatype(0,iwtvap)),
                                       tloc0(i,k), iwspec(idx_vap),
                                       iwtvap, iwtliq, false, 1.0, false);
            Real dliqiso;
            wtrc_liqvap_equil(alpha, 1.0,
                             qloc(i,k,wtrc_iatype(0,iwtvap)),
                             qloc(i,k,wtrc_iatype(0,iwtliq)),
                             qloc(i,k,idx_vap),
                             qloc(i,k,idx_liq),
                             dliqiso);
          }
          
          // Update temperature (H2O only)
          if (m == 0) {
            tloc(i,k) += (faltndqie - faltndi) * (latvap + latice) * dtime / cpair / nstep;
            tloc(i,k) += (faltndqce - faltndc) * latvap * dtime / cpair / nstep;
          }
        }
        
        // Update fractionation temperature
        tloc0(i,k) = tloc(i,k);
        
        // Reset fall velocities
        fice(k) = max(fice(k)/pdel(i,k), fice(k-1)/pdel(i,k-1)) * pdel(i,k);
        fcol(k) = max(fcol(k)/pdel(i,k), fcol(k-1)/pdel(i,k-1)) * pdel(i,k);
      }
      
      // Accumulate surface precipitation
      for (int m = 0; m < wtrc_nwset; ++m) {
        precr(i,m) += faloutc(PVER-1,m) / gravit / nstep / 1000.0;  // m/s
        preci(i,m) += falouti(PVER-1,m) / gravit / nstep / 1000.0;  // m/s
      }
    }
  }
  
  // Output precipitation to physics buffer
  for (int iwset = 0; iwset < wtrc_nwset; ++iwset) {
    auto precr_field = pbuf.get_field(wtrc_srfpcp_indices(iwtstrain, iwset));
    auto preci_field = pbuf.get_field(wtrc_srfpcp_indices(iwtstsnow, iwset));
    Kokkos::parallel_for("output_precip", ncol,
      KOKKOS_LAMBDA(int i) {
        precr_field(i) = precr(i,iwset);
        preci_field(i) = preci(i,iwset);
      });
  }
}
```

---

## Remaining Functions (Abbreviated)

Due to length constraints, here are implementation summaries for the other functions:

### Function 3: wtrc_apply_rates_mg1 (789 lines)
- **Most complex function**
- Iteratively applies pre-rates, calls wtrc_sediment_mg1, applies post-rates
- Handles Rayleigh distillation, Bergeron process, equilibration
- Two-pass numerical error correction at end
- See lines 1020-1809 in Fortran for full algorithm

### Function 4: wtrc_sediment (485 lines)  
- Extends wtrc_sediment_mg1 to handle 4 hydrometeor types (adds rain, snow)
- More complex CFL calculation
- Stewart isotope evaporation model for rain
- See lines 3138-3623 in Fortran

### Function 5: wtrc_apply_rates (629 lines)
- Similar to wtrc_apply_rates_mg1 but for MG2
- Calls wtrc_sediment instead of wtrc_sediment_mg1
- Different precipitation tracking
- See lines 2068-2697 in Fortran

### Function 6: wtrc_q1q2_pjr (311 lines)
- ZM deep convection interface
- Updraft/downdraft calculations with fractionation
- Compressed array indexing (ideep)
- See lines 6752-7063 in Fortran

### Function 7: wtrc_mass_fixer (267 lines)
- Reset H2O tracer to match bulk Q
- Proportional adjustment of isotopes
- See lines 4928-5195 in Fortran

### Function 8: wtrc_check_h2o (331 lines)
- Mass conservation checking
- Detailed error diagnostics
- See lines 4372-4703 in Fortran

---

## Testing Strategy

### Unit Tests
1. Test wtrc_mg_inter rate mapping with synthetic MG2 rates
2. Test wtrc_sediment_mg1 CFL calculation
3. Test equilibration in sedimentation

### Integration Tests
1. Single-column test with full MG microphysics
2. Compare outputs against Fortran for identical inputs
3. Check mass conservation

### GPU Tests
1. Verify Kokkos parallel patterns work on device
2. Check for race conditions
3. Performance profiling

---

## Next Steps

1. **Set up EAMxx infrastructure** - Get physics state/tendency/buffer working
2. **Implement wtrc_mg_inter** - Straightforward mapping, good first target
3. **Implement wtrc_sediment_mg1** - Core algorithm for sedimentation
4. **Implement wtrc_apply_rates_mg1** - Brings it all together
5. **Test and validate** - Compare against Fortran
6. **Complete remaining 5 functions** - MG2 support, convection, diagnostics

---

**Estimated Total Effort:** 6-8 weeks full-time with EAMxx infrastructure in place

**Status:** Ready to begin implementation once infrastructure is available
