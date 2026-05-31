#ifndef SCREAM_WATER_TRACER_RATIO_HPP
#define SCREAM_WATER_TRACER_RATIO_HPP

#include "share/scream_types.hpp"
#include <Kokkos_Core.hpp>

namespace scream {
namespace water_tracers {

//
// Compute ratio of tracer N to tracer 0 (bulk water)
//
// This utility is used for validating passive tracer transport.
// For correctly implemented tracers with scaled surface fluxes,
// ratios should be preserved to machine precision.
//
// @param tracer_amount: Amount of tracer N at a grid point
// @param bulk_amount: Amount of tracer 0 (bulk) at the same point
// @param min_bulk: Minimum bulk amount threshold (avoids division by zero)
// @return: Ratio tracer_N / tracer_0, or 0.0 if bulk < min_bulk
//
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT compute_tracer_ratio(
  const ScalarT& tracer_amount,
  const ScalarT& bulk_amount,
  const ScalarT& min_bulk = 1e-20
) {
  if (bulk_amount < min_bulk) {
    return 0.0;
  }
  return tracer_amount / bulk_amount;
}

//
// Compute tracer ratio field (parallel over columns and levels)
//
// @param tracer_field: 3D field with shape (tracer, col, lev)
// @param tracer_idx: Index of tracer to compute ratio for (typically 1)
// @param ratio_field: Output 2D field with shape (col, lev)
// @param min_bulk: Minimum bulk amount threshold
//
template<typename ViewT3D, typename ViewT2D>
void compute_tracer_ratio_field(
  const ViewT3D& tracer_field,  // (tracer, col, lev)
  const int tracer_idx,
  const ViewT2D& ratio_field,   // (col, lev) - output
  const Real min_bulk = 1e-20
) {
  const int ncols = tracer_field.extent(1);
  const int nlevs = tracer_field.extent(2);

  Kokkos::parallel_for(
    "compute_tracer_ratio_field",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(const int icol, const int ilev) {
      const Real bulk = tracer_field(0, icol, ilev);
      const Real tracer = tracer_field(tracer_idx, icol, ilev);
      ratio_field(icol, ilev) = compute_tracer_ratio(tracer, bulk, min_bulk);
    }
  );
}

//
// Compute global mean ratio for a 3D field
//
// @param tracer_field: 3D field with shape (tracer, col, lev)
// @param tracer_idx: Index of tracer to compute ratio for
// @param min_bulk: Minimum bulk amount threshold
// @return: Global mean ratio (averaged over all col, lev where bulk > min_bulk)
//
template<typename ViewT3D>
Real compute_mean_tracer_ratio(
  const ViewT3D& tracer_field,
  const int tracer_idx,
  const Real min_bulk = 1e-20
) {
  const int ncols = tracer_field.extent(1);
  const int nlevs = tracer_field.extent(2);

  Real sum_ratio = 0.0;
  int count = 0;

  Kokkos::parallel_reduce(
    "compute_mean_tracer_ratio",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(const int icol, const int ilev, Real& lsum, int& lcount) {
      const Real bulk = tracer_field(0, icol, ilev);
      if (bulk >= min_bulk) {
        const Real tracer = tracer_field(tracer_idx, icol, ilev);
        lsum += compute_tracer_ratio(tracer, bulk, min_bulk);
        lcount += 1;
      }
    },
    sum_ratio, count
  );

  return (count > 0) ? (sum_ratio / count) : 0.0;
}

} // namespace water_tracers
} // namespace scream

#endif // SCREAM_WATER_TRACER_RATIO_HPP
