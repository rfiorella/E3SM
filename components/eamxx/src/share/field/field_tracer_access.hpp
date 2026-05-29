#ifndef SCREAM_FIELD_TRACER_ACCESS_HPP
#define SCREAM_FIELD_TRACER_ACCESS_HPP

#include <Kokkos_Core.hpp>
#include "ekat_assert.hpp"

namespace scream {

/*
 * Accessor patterns for tracer-dimension fields
 *
 * Water species fields are extended from (col, lev) to (tracer, col, lev)
 * to support isotope tracers. Slot 0 always contains the bulk water species.
 *
 * Two accessor patterns are provided:
 *   1. EXPLICIT: Direct indexing into slot 0: qv(0, icol, ilev)
 *   2. SUBVIEW:  Kokkos subview of slot 0: auto qv_bulk = get_tracer_bulk_subview(qv)
 *
 * The choice is controlled by CMake flag SCREAM_TRACER_ACCESS=[EXPLICIT|SUBVIEW]
 * and is validated for BFB preservation in spec 2a tests.
 */

// --- EXPLICIT INDEXING PATTERN ---
// Returns a functor that accesses slot 0 via explicit indexing
// Usage: auto qv_accessor = get_tracer_bulk_explicit(qv_view);
//        Real value = qv_accessor(icol, ilev);  // internally does qv_view(0, icol, ilev)

template<typename ViewT>
struct ExplicitTracerAccessor {
  ViewT view;

  KOKKOS_INLINE_FUNCTION
  ExplicitTracerAccessor(const ViewT& v) : view(v) {}

  // 2D access (col, lev) - internally indexes slot 0
  KOKKOS_INLINE_FUNCTION
  typename ViewT::value_type& operator()(int icol, int ilev) const {
    return view(0, icol, ilev);
  }

  // Alternative with explicit tracer index (for generality)
  KOKKOS_INLINE_FUNCTION
  typename ViewT::value_type& operator()(int itracer, int icol, int ilev) const {
    EKAT_KERNEL_ASSERT_MSG(itracer == 0,
      "ExplicitTracerAccessor: bulk water is at slot 0");
    return view(itracer, icol, ilev);
  }
};

template<typename ViewT>
KOKKOS_INLINE_FUNCTION
ExplicitTracerAccessor<ViewT> get_tracer_bulk_explicit(const ViewT& view) {
  return ExplicitTracerAccessor<ViewT>(view);
}

// --- SUBVIEW PATTERN ---
// Returns a Kokkos subview of slot 0, collapsing the tracer dimension
// Usage: auto qv_bulk = get_tracer_bulk_subview(qv_view);
//        Real value = qv_bulk(icol, ilev);  // 2D view directly

template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk_subview(const ViewT& view)
  -> decltype(Kokkos::subview(view, 0, Kokkos::ALL(), Kokkos::ALL()))
{
  return Kokkos::subview(view, 0, Kokkos::ALL(), Kokkos::ALL());
}

// --- CONDITIONAL ACCESSOR ---
// Select pattern at compile time based on SCREAM_TRACER_ACCESS macro
// This is the recommended interface for physics code

#if defined(SCREAM_TRACER_ACCESS_EXPLICIT)

template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk(const ViewT& view) {
  return get_tracer_bulk_explicit(view);
}

#elif defined(SCREAM_TRACER_ACCESS_SUBVIEW)

template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk(const ViewT& view) {
  return get_tracer_bulk_subview(view);
}

#else
// Default to explicit indexing if not specified
template<typename ViewT>
KOKKOS_INLINE_FUNCTION
auto get_tracer_bulk(const ViewT& view) {
  return get_tracer_bulk_explicit(view);
}

#endif

} // namespace scream

#endif // SCREAM_FIELD_TRACER_ACCESS_HPP
