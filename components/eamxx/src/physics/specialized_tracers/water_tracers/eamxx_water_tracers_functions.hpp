#ifndef EAMXX_WATER_TRACERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACERS_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::abs for ekat::Pack

namespace scream {
namespace water_tracers {

/*
 * Device-callable functions for water tracer ratios and ocean evaporation.
 *
 * ScalarT must be an ekat::Pack<Real,N> for some N; N=1 stands in for the
 * scalar case (the same convention as eamxx_water_isotopes_fractionation.hpp).
 * Lanes are selected with ekat::Mask, and every division or comparison that
 * could see a zero, NaN or fill value is masked first, so these functions are
 * safe in FPE-trapping builds.
 *
 * Fluxes are upward positive [kg/m2/s]: the opposite sign of the coupler's evap.
 */

// Minimum |qtot| below which tracer_ratio() treats the denominator as unreliable
// and returns r_std instead of dividing. Guards the masked division against
// 0/0 on dry cells in FPE-trapping builds; not tied to any specific Fortran
// source (water_isotopes.F90::wiso_ratio uses a different, unmasked scheme
// with a six-orders-larger qtiny=1e-16, and its own comment flags that value
// as unreconciled with this one).
static constexpr Real water_tracers_qmin = 1.e-22;

// The coupler fills out-of-domain ocean states with spval = 1e30
// (share/util/shr_flux_mod.F90); anything this large is not a physical value.
static constexpr Real ocean_state_fill_threshold = 1.e20;

// How a tracer's evaporation flux is formed (see ocean_flux).
enum class OceanFlux { BulkProportional, BulkFormula };

// -----------------------------------------------------------------------
// Tracer ratio R = qtrc/qtot, or r_std where |qtot| < water_tracers_qmin.
//
// The denominator is replaced by 1 on the masked lanes BEFORE dividing, so
// a dry cell (qtot = qtrc = 0) never forms 0/0.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT tracer_ratio (const ScalarT& qtrc, const ScalarT& qtot,
                      const typename ekat::ScalarTraits<ScalarT>::scalar_type r_std = 1)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto tiny = ekat::abs(qtot) < RealT(water_tracers_qmin);

  ScalarT denom = qtot;
  denom.set(tiny, RealT(1));

  ScalarT ratio = qtrc / denom;
  ratio.set(tiny, r_std);
  return ratio;
}

// -----------------------------------------------------------------------
// Bulk ocean flux of a water species, upward positive:
//
//   F = f * rho * u * r * (s_eff - q_i)
//
// f = ocean fraction, rho = bottom-level air density, u = ocean friction
// velocity (So_ustar), r = moisture transfer coefficient (So_re, the square
// root of the exchange coefficient in the Large scheme), s_eff = the species'
// saturation specific humidity at the sea surface, q_i = the species at the
// bottom level. With s_eff = So_ssq and q_i = qv this is the coupler's
// evaporation (shr_flux_mod.F90: evap = rbot*ustar*re*(qbot - ssq)) with the
// opposite sign, times f because the imported states are not fraction
// weighted.
//
// Callers rely on this exact operation order: E_ocn is computed with it, so
// e.g. bulk_ocean_flux(f,rho,u,r,s*1,q) == E_ocn bit-for-bit.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT bulk_ocean_flux (const ScalarT& f, const ScalarT& rho, const ScalarT& u,
                         const ScalarT& r, const ScalarT& s_eff, const ScalarT& q_i)
{
  return f*rho*u*r*(s_eff - q_i);
}

// -----------------------------------------------------------------------
// Valid-ocean mask: f > 0, u > 0, and finite r, s with |r|, |s| < 1e20.
//
// NaN lanes are found with isnan and zeroed before any ordered comparison,
// because an ordered comparison with a NaN operand may raise FE_INVALID.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ekat::Mask<ScalarT::n> ocean_valid (const ScalarT& f, const ScalarT& u,
                                    const ScalarT& r, const ScalarT& s)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto nan = ekat::isnan(f) || ekat::isnan(u) || ekat::isnan(r) || ekat::isnan(s);

  ScalarT fz = f, uz = u, rz = r, sz = s;
  fz.set(nan, RealT(0));
  uz.set(nan, RealT(0));
  rz.set(nan, RealT(0));
  sz.set(nan, RealT(0));

  const RealT fill = ocean_state_fill_threshold;
  return !nan && fz > RealT(0) && uz > RealT(0) &&
         ekat::abs(rz) < fill && ekat::abs(sz) < fill;
}

// -----------------------------------------------------------------------
// Bulk ocean evaporation E_ocn = bulk_ocean_flux(f, rho, u, r, s, q).
//
// On return f, u, r and s are zeroed wherever the ocean state is not valid
// (see ocean_valid), so no later arithmetic touches coupler fill values, and
// E_ocn is exactly 0 there. Downstream flux functions use f > 0 as the
// valid-ocean test.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ocean_evaporation (ScalarT& f, ScalarT& u, ScalarT& r, ScalarT& s,
                           const ScalarT& rho, const ScalarT& q)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto invalid = !ocean_valid(f, u, r, s);
  f.set(invalid, RealT(0));
  u.set(invalid, RealT(0));
  r.set(invalid, RealT(0));
  s.set(invalid, RealT(0));

  ScalarT e_ocn(0);
  e_ocn.set(!invalid, bulk_ocean_flux(f, rho, u, r, s, q));
  return e_ocn;
}

// -----------------------------------------------------------------------
// True where f marks a valid ocean point and e_ocn is evaporation (not dew).
// Factored out so a caller that must compute something only on evaporating
// lanes (e.g. an isotope's temperature-dependent fractionation factor) can
// build the same mask ocean_flux uses internally, without a second copy of
// the boolean expression.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ekat::Mask<ScalarT::n> is_evaporating (const ScalarT& f, const ScalarT& e_ocn)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  return f > RealT(0) && e_ocn > RealT(0);
}

// -----------------------------------------------------------------------
// Ocean flux of a tracer with ocean ratio r_ocn and equilibrium fractionation
// factor alpha (default 1, i.e. no fractionation):
//
//   not valid ocean (f <= 0):  0
//   E_ocn > 0, BulkProportional: r_ocn * E_ocn
//   E_ocn > 0, BulkFormula:      bulk_ocean_flux(f, rho, u, r, s*r_ocn/alpha, q_i)
//   E_ocn <= 0 (dew):            tracer_ratio(q_i, q) * E_ocn
//
// alpha is used only in the BulkFormula evap branch; BulkProportional and dew
// are unaffected by it, so a caller that omits alpha (every generic tracer,
// and any isotope species or build where alpha == 1, e.g. H216O) gets exactly
// today's BulkFormula result bit-for-bit. Dew removes the tracer at its
// bottom-level ratio to qv.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ocean_flux (const OceanFlux kind,
                    const typename ekat::ScalarTraits<ScalarT>::scalar_type r_ocn,
                    const ScalarT& f, const ScalarT& rho, const ScalarT& u,
                    const ScalarT& r, const ScalarT& s, const ScalarT& q,
                    const ScalarT& q_i, const ScalarT& e_ocn,
                    const ScalarT& alpha = ScalarT(1))
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto ocean = f > RealT(0);
  const auto evap  = is_evaporating(f, e_ocn);
  const auto dew   = ocean && !(e_ocn > RealT(0));

  ScalarT flux(0);
  if (evap.any()) {
    if (kind == OceanFlux::BulkProportional) {
      flux.set(evap, r_ocn*e_ocn);
    } else {
      flux.set(evap, bulk_ocean_flux(f, rho, u, r, s*r_ocn/alpha, q_i));
    }
  }
  if (dew.any()) {
    flux.set(dew, tracer_ratio(q_i, q)*e_ocn);
  }
  return flux;
}

// -----------------------------------------------------------------------
// Total surface flux of a tracer, upward positive [kg/m2/s]: the ocean part
// (ocean_flux) plus the land and sea-ice part, where
//
//   s_evap = merged ocean+land+ice surface evaporation of bulk water
//            (the surf_evap field, upward positive), assumed finite
//   e_ocn  = the ocean-only bulk evaporation (ocean_evaporation), exactly 0
//            where the ocean state is not valid
//   L      = s_evap - e_ocn, the non-ocean bulk evaporation
//
// L > 0 evaporates the tracer at the land ratio r_land; L < 0 (dew) removes it
// at its bottom-level ratio tracer_ratio(q_i, q), like the ocean dew branch.
//
// Cases, for the sign of e_ocn (ocean) and L (land and ice), and the flux
// formed for a BulkProportional tracer (rho_d = tracer_ratio(q_i, q)):
//
//   e_ocn >= 0, L >= 0  both evaporate: r_land*s_evap + (r_ocn - r_land)*e_ocn
//                       (the same value as r_ocn*e_ocn + r_land*L, but with
//                       r_ocn = r_land = 1 it is exactly s_evap, so a tracer
//                       with R = 1 closes against bulk water to the bit)
//   e_ocn >= 0, L <  0  ocean evaporates, land dews: r_ocn*e_ocn + rho_d*L
//   e_ocn <  0, L >= 0  ocean dews, land evaporates: rho_d*e_ocn + r_land*L
//   e_ocn <  0, L <  0  both dew: rho_d*e_ocn + rho_d*L
//
// (A column with no valid ocean has e_ocn = 0, so it is the first or second
// row with r_land only.) A tracer with R = 1 everywhere closes to roundoff in
// the dew cases, because rho_d = 1 there.
//
// BulkFormula tracers (and isotopes) form the ocean part with the ocean bulk
// formula on the tracer's own bottom-level concentration q_i, so they do not
// close exactly against s_evap: only the non-ocean part is proportional. The
// ocean part and the non-ocean part are added as in the lower three rows.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT surface_tracer_flux (const OceanFlux kind,
                             const typename ekat::ScalarTraits<ScalarT>::scalar_type r_ocn,
                             const typename ekat::ScalarTraits<ScalarT>::scalar_type r_land,
                             const ScalarT& f, const ScalarT& rho, const ScalarT& u,
                             const ScalarT& r, const ScalarT& s, const ScalarT& q,
                             const ScalarT& q_i, const ScalarT& s_evap, const ScalarT& e_ocn,
                             const ScalarT& alpha = ScalarT(1))
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const ScalarT f_ocn = ocean_flux(kind, r_ocn, f, rho, u, r, s, q, q_i, e_ocn, alpha);

  const ScalarT l_evap = s_evap - e_ocn;
  const auto land_evap = l_evap > RealT(0);
  const auto land_dew  = l_evap < RealT(0);

  ScalarT flux = f_ocn;
  flux.set(land_evap, f_ocn + r_land*l_evap);
  if (land_dew.any()) {
    flux.set(land_dew, f_ocn + tracer_ratio(q_i, q)*l_evap);
  }
  if (kind == OceanFlux::BulkProportional) {
    const auto both_evap = !(e_ocn < RealT(0)) && !land_dew;
    flux.set(both_evap, r_land*s_evap + (r_ocn - r_land)*e_ocn);
  }
  return flux;
}

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACERS_FUNCTIONS_HPP
