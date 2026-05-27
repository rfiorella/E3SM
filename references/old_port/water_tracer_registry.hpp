#ifndef WATER_TRACER_REGISTRY_HPP
#define WATER_TRACER_REGISTRY_HPP

#include "water_tracer_registry.gen.hpp"
#include <optional>
#include <Kokkos_Core.hpp>

namespace scream {
namespace WaterTracers {

// Query API: tracer index → isotopologue catalog index
// Constexpr, callable from device kernels and host code
KOKKOS_INLINE_FUNCTION
constexpr int tracer_isotopologue(int tracer_idx) {
  return WaterTracerRegistry[tracer_idx].catalog_idx;
}

// Query API: tracer index → per-tracer name
// Constexpr, callable from device kernels and host code
KOKKOS_INLINE_FUNCTION
constexpr std::string_view tracer_name(int tracer_idx) {
  return WaterTracerRegistry[tracer_idx].name;
}

// Query API: tracer index → is-tag flag
// Constexpr, callable from device kernels and host code
KOKKOS_INLINE_FUNCTION
constexpr bool tracer_is_tag(int tracer_idx) {
  return WaterTracerRegistry[tracer_idx].is_tag;
}

// Host-only convenience: name → tracer index
// Returns std::nullopt if name not found
inline std::optional<int> find_tracer_by_name(std::string_view name) {
  for (int i = 0; i < WTRC_MAX_CNST; ++i) {
    if (WaterTracerRegistry[i].name == name) {
      return i;
    }
  }
  return std::nullopt;
}

} // namespace WaterTracers
} // namespace scream

#endif // WATER_TRACER_REGISTRY_HPP
