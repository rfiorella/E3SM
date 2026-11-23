#ifndef EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP

#include "eamxx_water_tracers_functions.hpp"  // tracer_ratio

namespace scream {
namespace water_tracers {

/*
 * Device-callable transfers of the water_tracer_transfers process, one
 * function per class of bulk water transfer, so that each class can get its
 * own fractionation factor later.
 *
 * Each function takes the donor phase's tracer and bulk water, and the bulk
 * mass moved out of the donor phase, and returns the tracer mass moved
 * [kg/kg]. Same ScalarT convention as eamxx_water_tracers_functions.hpp (an
 * ekat::Pack; tracer_ratio masks dry donors, so these are FPE-safe).
 */

// -----------------------------------------------------------------------
// SHOC's net condensation `amount`: > 0 moves vapor to cloud liquid (the
// donor is vapor), < 0 moves cloud liquid to vapor (the donor is cloud
// liquid, and the result has the sign of amount too).
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT shoc_condensation (const ScalarT& tr_donor, const ScalarT& q_donor, const ScalarT& amount)
{
  // alpha = 1: the tracer moves at the donor's ratio. This is where the
  // liquid/vapor fractionation factor of SHOC's phase change goes.
  return amount*tracer_ratio(tr_donor, q_donor);
}

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
