# Water Tracer Registry

## Overview

The water tracer registry provides a compile-time mapping between tracer indices (CMP dimension slots) and isotopologue catalog indices (physical parameters in `water_isotopes.hpp`). This enables type-safe, efficient lookup of isotopologue metadata from device kernels.

## Quick Start

### Configuring Tracers

Create a tracer configuration file:

```cmake
# my_tracers.cmake
add_water_tracer(NAME bulk ISOTOPOLOGUE H2O)
add_water_tracer(NAME hdo ISOTOPOLOGUE HDO)
add_water_tracer(NAME h218o ISOTOPOLOGUE H218O)
```

Build with:
```bash
cmake -S . -B build \
  -DSCREAM_TRACE_WATER=ON \
  -DSCREAM_WATER_TRACERS_FILE=/path/to/my_tracers.cmake
```

### Using the Registry API

```cpp
#include "physics/water_tracers/water_tracer_registry.hpp"
#include "physics/water_tracers/water_isotopes.hpp"

using namespace scream::WaterTracers;
using namespace scream::WaterIsotopes;

// Device kernel
KOKKOS_INLINE_FUNCTION
void apply_fractionation(int tracer_idx, Real temperature) {
  // Get catalog index for this tracer
  int cat_idx = tracer_isotopologue(tracer_idx);
  
  // Skip H2O (no self-fractionation)
  if (cat_idx == 0) return;
  
  // Skip passive tags
  if (tracer_is_tag(tracer_idx)) return;
  
  // Get isotopologue name from catalog
  auto iso_name = WaterIsotopologueNames[cat_idx];
  
  // Compute fractionation factor
  Real alpha = AlphaEqLiquidVapor(iso_name, temperature);
  
  // Apply to tendency...
}
```

## API Reference

### Query Functions

All query functions are `constexpr` and callable from device kernels (except `find_tracer_by_name`).

#### `int tracer_isotopologue(int i)`
Returns catalog index for tracer `i`. Maps to `WaterIsotopologues<Scalar>` arrays in `water_isotopes.hpp`.

**Example:**
```cpp
int cat_idx = tracer_isotopologue(2);  // e.g., returns 2 for HDO
```

#### `std::string_view tracer_name(int i)`
Returns per-tracer name (unique identifier from config file).

**Example:**
```cpp
auto name = tracer_name(0);  // "bulk"
```

#### `bool tracer_is_tag(int i)`
Returns true if tracer is a passive tag (no fractionation).

**Example:**
```cpp
if (tracer_is_tag(i)) {
  // Skip fractionation for tags
}
```

#### `std::optional<int> find_tracer_by_name(std::string_view name)` (HOST ONLY)
Returns tracer index for given name, or `std::nullopt` if not found.

**Example:**
```cpp
auto idx = find_tracer_by_name("hdo");
if (idx.has_value()) {
  // Use idx.value()
}
```

### Data Structures

#### `WaterTracerInfo`
```cpp
struct WaterTracerInfo {
  std::string_view name;     // Per-tracer name
  int catalog_idx;            // Index into WaterIsotopologues arrays
  bool is_tag;                // True for passive tags
};
```

#### `WaterTracerRegistry`
```cpp
constexpr std::array<WaterTracerInfo, WTRC_MAX_CNST> WaterTracerRegistry;
```

## Configuration Reference

### `add_water_tracer()`

```cmake
add_water_tracer(
  NAME <unique_name>
  ISOTOPOLOGUE <catalog_name>
  [TAG]  # Optional: marks as passive tag
)
```

**Parameters:**
- `NAME`: Unique identifier for this tracer slot
- `ISOTOPOLOGUE`: One of: `H2O`, `H216O`, `HDO`, `H218O`, `H217O`, `HTO`
- `TAG`: Optional flag for passive tracers (no fractionation)

**Validation (configure-time):**
- First tracer must be H2O or H216O
- All isotopologue names must be in catalog
- No duplicate tracer names
- At least one tracer required when `SCREAM_TRACE_WATER=ON`

### Example Configurations

#### Single bulk H2O:
```cmake
add_water_tracer(NAME bulk ISOTOPOLOGUE H2O)
```

#### Bulk + isotopologues:
```cmake
add_water_tracer(NAME bulk ISOTOPOLOGUE H2O)
add_water_tracer(NAME hdo ISOTOPOLOGUE HDO)
add_water_tracer(NAME h218o ISOTOPOLOGUE H218O)
```

#### Bulk + tag + isotopologue:
```cmake
add_water_tracer(NAME bulk ISOTOPOLOGUE H2O)
add_water_tracer(NAME evap_tag ISOTOPOLOGUE H2O TAG)  # Passive tag
add_water_tracer(NAME hdo ISOTOPOLOGUE HDO)
```

## Index Spaces

Two distinct index spaces are used:

### Tracer Index (CMP dimension)
- Range: `[0, WTRC_MAX_CNST)`
- Build-dependent (set by config file)
- Addresses slots in `(COL, CMP, LEV)` field layout

### Catalog Index
- Range: `[0, 5]` (fixed)
- Codebase-constant (defined in `water_isotopes.hpp`)
- Addresses `WaterIsotopologues<Scalar>` constant arrays

**The registry maps tracer index → catalog index.**

## Important Notes

1. **Catalog ordering is fixed.** Appending new isotopologues is allowed; reordering breaks existing configs.

2. **Multiple tracers can share a catalog index.** Example: bulk H2O and passive H2O tags both map to catalog index 0.

3. **Tracer 0 must be H2O.** This is enforced by validation because bulk-water consumers hard-code CMP=0.

4. **No runtime overhead.** All queries are constexpr; lookups compile to direct array indexing.

5. **Device-callable.** Use freely in Kokkos kernels - the registry is constexpr device data.

## Migration from SCREAM_NUM_WATER_TRACERS

**Old (deprecated):**
```bash
cmake -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=3
```

**New:**
```bash
cmake -DSCREAM_TRACE_WATER=ON -DSCREAM_WATER_TRACERS_FILE=path/to/config.cmake
```

The count (`WTRC_MAX_CNST`) is derived from the config file. Setting `SCREAM_NUM_WATER_TRACERS` now triggers a hard error.

## Files

- **Public API:** `water_tracer_registry.hpp`
- **Generated (do not edit):** `build/water_tracer_registry.gen.hpp`
- **Template:** `cmake/water_tracers/water_tracer_registry.gen.hpp.in`
- **Example configs:** `cmake/water_tracers/bulk_only.cmake`, `registry_n4.cmake`
- **Helper scripts:** `scripts/water-tracers-config-flags`

## See Also

- `water_isotopes.hpp` - Isotopologue catalog and fractionation parameters
- `water_tracers.hpp` - Field layout and subview accessors
- `REGISTRY_FOLLOWUP.md` - Next steps (fractionation physics, tags, etc.)
- Spec: `specs/2026-05-06-eamxx-water-tracer-isotopologue-registry.md`
