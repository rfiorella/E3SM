# Water Tracer Conversion from CAM to EAMxx

## Overview

This document describes the conversion of water tracer modules from CAM (Fortran) to EAMxx (C++/Kokkos). The water tracer system enables tracking of water isotopes (HDO, H218O, H217O, HTO) and generic water tracers through the atmosphere model with proper isotopic fractionation.

## Converted Modules

### 1. water_isotopes.hpp
**Original**: `water_isotopes.F90` (775 lines)  
**Location**: `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_isotopes.hpp`

**Key Features**:
- Isotope species definitions (H2O, H216O, HDO, H218O, H217O, HTO)
- Equilibrium fractionation factors (liquid-vapor, ice-vapor)
- Kinetic fractionation effects (Merlivat & Jouzel, 1979)
- Ocean evaporation calculations
- Radioactive decay (for tritium)
- Device-compatible (GPU/CPU) with `KOKKOS_INLINE_FUNCTION`
- Works with EKat Pack types for vectorization

**Scientific References**:
- Horita & Wesolowski (1994) - Liquid/vapor fractionation
- Merlivat & Nief (1967) - Ice/vapor fractionation
- Jouzel & Merlivat (1984) - Kinetic effects
- Stewart (1975) - Diffusive evaporation
- Majoube (1971a,b) - Fractionation coefficients
- Schoenemann et al. (2014) - H217O adjustments

**Main Functions**:
```cpp
// Equilibrium fractionation
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk);  // Liquid-vapor

template<typename ScalarT>
ScalarT wiso_alpi(int isp, const ScalarT& tk);  // Ice-vapor

// Kinetic fractionation
template<typename ScalarT>
ScalarT wiso_kmol(int isp, const ScalarT& rbot, const ScalarT& zbot, const ScalarT& ustar);

template<typename ScalarT>
ScalarT wiso_akel(int isp, const ScalarT& tk, const ScalarT& hum0, const ScalarT& alpeq);

template<typename ScalarT>
ScalarT wiso_akci(int isp, const ScalarT& tk, const ScalarT& alpeq, const ScalarT& rh = -1.0);

// Ocean evaporation
template<typename ScalarT>
ScalarT wiso_flxoce(/* ... */);

// Utilities
template<typename ScalarT>
ScalarT wiso_ratio(int isp, const ScalarT& qiso, const ScalarT& qtot);

template<typename ScalarT>
ScalarT wiso_delta(int isp, const ScalarT& qiso, const ScalarT& qtot);

template<typename ScalarT>
ScalarT wiso_decay(int isp, const ScalarT& q, Real dtime);
```

### 2. water_types.hpp
**Original**: `water_types.F90` (170 lines)  
**Location**: `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_types.hpp`

**Key Features**:
- Water type enumerations (Vapor, CloudLiquid, CloudIce, StratRain, StratSnow, ConvRain, ConvSnow)
- Phase transition fractionation logic
- Device-compatible utility functions
- Type checking functions (is_liquid, is_ice, is_precip, etc.)

**Main Functions**:
```cpp
// Get fractionation factor for phase transitions
template<typename ScalarT>
ScalarT wtype_get_alpha(int ispec, int isrctype, int idsttype,
                        const ScalarT& tk, const ScalarT& rh,
                        bool do_kinetic, const ScalarT& rhi = -1.0);

// Convenience functions
template<typename ScalarT>
ScalarT wtype_get_alpha_condensation(/* ... */);

template<typename ScalarT>
ScalarT wtype_get_alpha_evaporation(/* ... */);

// Utility functions (device-compatible)
KOKKOS_INLINE_FUNCTION
bool wtype_is_liquid(int itype);

KOKKOS_INLINE_FUNCTION
bool wtype_is_ice(int itype);

KOKKOS_INLINE_FUNCTION
bool wtype_is_precip(int itype);
```

### 3. CMake Integration

**Files Modified**:
1. `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt`
2. `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/share/core/eamxx_config.h.in`
3. `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/CMakeLists.txt`

**CMake Options**:
```cmake
option(SCREAM_TRACE_WATER "Whether to enable generic water tracers" OFF)
option(SCREAM_WATER_ISOTOPES "Whether to enable water isotope tracers" OFF)
```

**Preprocessor Definitions**:
```cpp
#ifdef SCREAM_TRACE_WATER
  // Generic water tracer code
#endif

#ifdef SCREAM_WATER_ISOTOPES
  // Isotope-specific code
#endif
```

## Key Design Decisions

### 1. Device Portability
All physics functions use `KOKKOS_INLINE_FUNCTION` to enable GPU execution. Functions are templated on scalar type to work with both:
- Scalar values: `Real` (double or float)
- Pack types: `ekat::Pack<Real,N>` for vectorization

### 2. Conditional Compilation
Water tracer features are enabled via CMake options:
- `SCREAM_TRACE_WATER`: Adds tracer dimension to water fields
- `SCREAM_WATER_ISOTOPES`: Enables isotope-specific calculations (requires TRACE_WATER)

### 3. Namespace Organization
```
scream::
  ├── water_isotopes::   // Isotope fractionation functions
  └── water_types::      // Water phase definitions and transitions
```

### 4. Type Safety
- Enum classes for species (`IsoSpecies`) and water types (`WaterType`)
- Constexpr arrays for physical constants
- Strong typing prevents mixing of species and constituent indices

## Usage Example

```cpp
#include "water_isotopes.hpp"
#include "water_types.hpp"

using namespace scream;
using namespace scream::water_isotopes;
using namespace scream::water_types;

// Calculate liquid-vapor fractionation for HDO at 280K
constexpr int isp_hdo = static_cast<int>(IsoSpecies::HDO) - 1;  // 0-indexed
Real temp = 280.0;
Real alpha_eq = wiso_alpl(isp_hdo, temp);

// Apply kinetic effects
Real rh = 0.8;  // 80% relative humidity
Real alpha_kin = wiso_akel(isp_hdo, temp, rh, alpha_eq);

// Get fractionation for vapor -> liquid condensation
constexpr int src = static_cast<int>(WaterType::Vapor);
constexpr int dst = static_cast<int>(WaterType::CloudLiquid);
Real alpha = wtype_get_alpha(isp_hdo, src, dst, temp, rh, true);
```

## Integration with Physics Processes

Water tracers need to be integrated into physics processes that handle water:
1. **Microphysics** (P3, MG) - condensation, evaporation, sedimentation
2. **Convection** (ZM, SHOC) - convective transport and detrainment
3. **Turbulence** (SHOC, vertical diffusion) - turbulent mixing
4. **Surface fluxes** - ocean/land evaporation with fractionation

Each process should:
1. Check if `SCREAM_TRACE_WATER` is defined
2. Apply fractionation factors from `water_types::wtype_get_alpha()`
3. Update tracer concentrations alongside bulk water

## Testing Recommendations

1. **Unit Tests**: Test fractionation functions against known values
   - Equilibrium fractionation vs. published data
   - Kinetic fractionation in different regimes
   
2. **Single-Column Tests**: Verify isotope conservation
   - Total mass conservation
   - Isotope ratio maintenance
   
3. **GPU Compatibility**: Verify device execution
   - Pack-based operations
   - Kokkos parallel patterns

## Next Steps

1. **Main water_tracers Module**: Convert `water_tracers.F90` (7,813 lines)
   - Process rate tracking and application
   - Sedimentation with tracers
   - Precipitation tracking
   - Mass conservation checks
   
2. **Physics Integration**: Update physics processes
   - P3 microphysics
   - SHOC convection/turbulence
   - Surface flux calculations
   
3. **Field Management**: Add tracer dimension to water fields
   - Modify field definitions in AtmosphereProcess
   - Update field requests to include tracer dimension
   - Conditional dimension based on SCREAM_TRACE_WATER
   
4. **Initialization**: Add tracer initialization
   - Read initial conditions
   - Set up tracer metadata
   - Configure output variables

## Files Changed

### Modified:
- `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt`
  - Added `SCREAM_TRACE_WATER` and `SCREAM_WATER_ISOTOPES` options
  - Added dependency check (isotopes require trace_water)

- `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/share/core/eamxx_config.h.in`
  - Added preprocessor defines for water tracer options

- `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/CMakeLists.txt`
  - Updated to conditionally define preprocessor macros

### Created/Rewritten:
- `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_isotopes.hpp`
  - Complete rewrite with device compatibility
  - All fractionation functions implemented
  - Template-based for scalar/pack compatibility

- `/Users/rfiorella/code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_types.hpp`
  - Complete rewrite with fractionation logic
  - Phase transition handling
  - Utility functions for type checking

## Differences from CAM Implementation

1. **Language**: Fortran → C++/Kokkos
2. **Device Support**: Added GPU compatibility
3. **Vectorization**: Uses EKat Pack types instead of explicit SIMD
4. **Type Safety**: Enum classes instead of integer parameters
5. **Namespaces**: Better organization than Fortran modules
6. **Templates**: Single implementation works for scalars and packs
7. **Conditional Compilation**: CMake-based instead of CPP defines

## Authors

**Original CAM Code**:
- David Noone <dcn@colorado.edu> (2003)
- Jesse Nusbaumer <nusbaume@colorado.edu> (2011)
- Chuck Bardeen (2012)

**EAMxx Conversion**:
- Rich Fiorella (2026)

## References

1. Horita, J., & Wesolowski, D. J. (1994). Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature. *Geochimica et Cosmochimica Acta*, 58(16), 3425-3437.

2. Merlivat, L., & Nief, G. (1967). Fractionnement isotopique lors des changements d'état solide-vapeur et liquide-vapeur de l'eau à des températures inférieures à 0° C. *Tellus*, 19(1), 122-127.

3. Jouzel, J., & Merlivat, L. (1984). Deuterium and oxygen 18 in precipitation: Modeling of the isotopic effects during snow formation. *Journal of Geophysical Research: Atmospheres*, 89(D7), 11749-11757.

4. Majoube, M. (1971a). Fractionnement en oxygène 18 et en deutérium entre l'eau et sa vapeur. *Journal de Chimie Physique*, 68, 1423-1436.

5. Majoube, M. (1971b). Fractionnement en 18O entre la glace et la vapeur d'eau. *Journal de Chimie Physique*, 68, 625-636.

6. Schoenemann, S. W., Schauer, A. J., & Steig, E. J. (2013). Measurement of SLAP2 and GISP δ17O and proposed VSMOW-SLAP normalization for δ17O and 17Oexcess. *Rapid Communications in Mass Spectrometry*, 27(5), 582-590.

7. Stewart, M. K. (1975). Stable isotope fractionation due to evaporation and isotopic exchange of falling waterdrops: Applications to atmospheric processes and evaporation of lakes. *Journal of Geophysical Research*, 80(9), 1133-1146.

8. Merlivat, L., & Jouzel, J. (1979). Global climatic interpretation of the deuterium‐oxygen 18 relationship for precipitation. *Journal of Geophysical Research: Oceans*, 84(C8), 5029-5033.
