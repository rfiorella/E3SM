#ifndef EAMXX_WATER_ISOTOPES_KINETIC_HPP
#define EAMXX_WATER_ISOTOPES_KINETIC_HPP

#include "share/core/eamxx_types.hpp"
#include "share/physics/physics_constants.hpp"
#include "eamxx_water_isotopes_parameters.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::log/ekat::pow overloads for ekat::Pack

namespace scream {
namespace water_isotopes {

/*
 * Kinetic isotope fractionation: of ocean evaporation (alpha_kinetic), and of
 * liquid evaporation into subsaturated air and deposition on ice from
 * supersaturated air (alpha_liquid_vapor_kinetic, alpha_ice_vapor_kinetic,
 * further down).
 *
 * Kinetic fractionation factor alpkn of ocean evaporation.
 *
 * Port of iCAM's wiso_kmol (share/util/water_isotopes.F90, original author
 * David Noone): Brutsaert (1975) turbulent/molecular resistance ratio with the
 * Merlivat & Jouzel (1979) kinetic modifier,
 *
 *   alpkn = 1 - kmol,  kmol = (difn - 1) / (difn + tmr),  difn = (D/D_i)^n,
 *
 * with a smooth (n = 2/3) or rough (n = 1/2) regime chosen by the surface
 * Reynolds number against recrit = 1. alpkn multiplies the bulk evaporation
 * flux of the isotope (see water_tracers::ocean_flux). It lies in (0,1] for
 * ustar*zbot >~ 1e-5 m2/s; below about 9e-6 m2/s the smooth-regime log term
 * takes it out of that range, as in the Fortran (which has no guard either).
 *
 * Each species uses its own diffusivity ratio D_i/D (iso_params.difrm, the
 * Fortran's difrm), H217O and HTO included, as wiso_kmol does. H216O has
 * alpkn = 1, as does every species under KineticFractionationFormulation::None.
 *
 * Same conventions as eamxx_water_isotopes_fractionation.hpp: ScalarT is an
 * ekat::Pack<Real,N> (N=1 for scalar use), and log()/pow() are called
 * unqualified so ADL picks the ekat::Pack overloads.
 */
struct WaterIsotopeKinetic
{
private:
  using PC = scream::physics::Constants<Real>;

  // Constants of wiso_kmol
  static constexpr double difair  = 2.36e-5;  // molecular diffusivity of air [m2/s]
  static constexpr double muair   = 1.7e-5;   // dynamic viscosity of air [kg/m/s]
  static constexpr double recrit  = 1.0;      // critical Reynolds number
  static constexpr double charnock_inv = 81.1;  // z0 = ustar^2/(81.1*g)

  // In-range state substituted into lanes outside the mask before any log or
  // pow, so that, e.g., ustar = 0 at invalid-ocean points cannot raise an FPE
  // (the same idea as WaterIsotopeFractionation's T_lane_fill).
  static constexpr double ustar_lane_fill = 0.3;   // [m/s]
  static constexpr double rbot_lane_fill  = 1.2;   // [kg/m3]
  static constexpr double zbot_lane_fill  = 10.0;  // [m]

  // kmol = 1 - alpkn of HDO (Hydrogen) or H218O (Oxygen), from wiso_kmol.
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT kmol(const ScalarT& rbot, const ScalarT& zbot, const ScalarT& ustar,
                      const typename ekat::ScalarTraits<ScalarT>::scalar_type difrmj)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    const ScalarT z0   = (ustar*ustar)/(RealT(charnock_inv)*RealT(PC::gravit.value));  // Charnock
    const ScalarT vmu  = RealT(muair)/rbot;   // kinematic viscosity
    const ScalarT Sc   = vmu/RealT(difair);   // Schmidt number
    const ScalarT reno = ustar*z0/vmu;        // Reynolds number

    const auto smooth = reno < RealT(recrit);
    const auto rough  = !smooth;

    ScalarT enn, tmr;
    if (smooth.any()) {
      enn.set(smooth, RealT(2)/RealT(3));
      tmr.set(smooth, ((RealT(1)/RealT(PC::Karman))*log(ustar*zbot/(RealT(30)*vmu)))
                      / (RealT(13.6)*pow(Sc, RealT(2)/RealT(3))));
    }
    if (rough.any()) {
      enn.set(rough, RealT(1)/RealT(2));
      tmr.set(rough, ((RealT(1)/RealT(PC::Karman))*log(zbot/z0) - RealT(5))
                     / (RealT(7.3)*pow(reno, RealT(1)/RealT(4))*pow(Sc, RealT(1)/RealT(2))));
    }

    const ScalarT difn = pow(ScalarT(RealT(1)/difrmj), enn);  // D/D_i, not D_i/D
    return (difn - RealT(1))/(difn + tmr);
  }

public:
  // -----------------------------------------------------------------------
  // alpkn for one species, given the bottom-level air density rbot [kg/m3],
  // bottom-level midpoint height zbot [m] and friction velocity ustar [m/s].
  //
  // mask selects the lanes holding real data (callers pass
  // water_tracers::is_evaporating(f, e_ocn)). Lanes outside it get in-range
  // placeholder inputs before any log/pow and return alpkn = 1.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_kinetic(
      const ScalarT& rbot,
      const ScalarT& zbot,
      const ScalarT& ustar,
      const WaterIsotopologues species,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& mask = ekat::Mask<ScalarT::n>(true))
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    if (iso_params.kinetic_formulation() == KineticFractionationFormulation::None ||
        species == WaterIsotopologues::H216O || !mask.any()) {
      return ScalarT(1);
    }

    ScalarT rbot_live{RealT(rbot_lane_fill)}, zbot_live{RealT(zbot_lane_fill)},
            ustar_live{RealT(ustar_lane_fill)};
    rbot_live.set(mask, rbot);
    zbot_live.set(mask, zbot);
    ustar_live.set(mask, ustar);

    const ScalarT k = kmol(rbot_live, zbot_live, ustar_live, iso_params.difrm(species));

    ScalarT alpkn(1);
    alpkn.set(mask, RealT(1) - k);
    return alpkn;
  }

  /*
   * Kinetic fractionation of liquid evaporation and of deposition on ice.
   *
   * Ports of wiso_heff, wiso_ssatf, wiso_akel and wiso_akci
   * (share/util/water_isotopes.F90, original author David Noone; RH over ice
   * in wiso_akci by Marina Dutsch). Each returns an effective fractionation
   * factor that replaces the equilibrium one, alpha_eq, in the same
   * CondensedOverVapor convention (R_condensed/R_vapor), with D_i/D the
   * molecular diffusivity ratio of the species (iso_params.difrm):
   *
   *  - alpha_liquid_vapor_kinetic (akel), liquid evaporating into air of
   *    relative humidity rh over liquid (Stewart 1975, J. Geophys. Res. 80,
   *    1133-1146):
   *      akel = alpha_eq*h/(alpha_eq*(D/D_i)^dkfac*(h - 1) + 1),
   *      h = effective_humidity(rh) = min(1, fkhum*min(1, rh) + 1 - fkhum).
   *    akel >= alpha_eq, = alpha_eq at rh >= 1, and grows as rh falls: the
   *    vapor the liquid gives off is depleted beyond equilibrium (it leaves
   *    with R_liquid/akel), since the heavy isotopologues diffuse more slowly.
   *  - alpha_ice_vapor_kinetic (akci), vapor depositing on ice at a
   *    supersaturation ratio S over ice (Jouzel & Merlivat 1984, J. Geophys.
   *    Res. 89, 11749-11757):
   *      akci = alpha_eq*S/(alpha_eq*(D/D_i)*(S - 1) + 1),  S in [1, ssatmx],
   *    applied only below tkini = 253.15 K, with the Fortran's hard step (=
   *    alpha_eq at and above tkini). akci <= alpha_eq, = alpha_eq at S = 1,
   *    and falls as S grows (below 1 near S = 2 for H218O): growth limited by
   *    diffusion takes the heavy isotopologues in less than equilibrium would.
   *    The caller gives S: the Fortran default is ice_supersaturation_empirical
   *    (T) (wiso_ssatf), the alternative the model's RH over ice.
   *
   * Every species uses its own alpha_eq and its own D_i/D, H217O and HTO
   * included, as wiso_akel and wiso_akci do (the equilibrium alpha_eq of
   * H217O and HTO keep their mass-scaling from H218O and HDO; the kinetic
   * part does not). H216O returns alpha_eq (= 1).
   *
   * Lanes outside mask return alpha_eq unchanged and get in-range inputs
   * before any division or pow, so padding cannot raise an FPE. Same ScalarT
   * convention as above. The constants are the Fortran's.
   */
  static constexpr double dkfac  = 0.58;    // exponent of D/D_i in akel (Stewart 1975)
  static constexpr double fkhum  = 0.25;    // weight of rh in the effective humidity
  static constexpr double tkini  = 253.15;  // [K] akci applies below it (Jouzel & Merlivat 1984)
  static constexpr double fsata  = 1.0;     // ssatf = fsata + fsatb*(T - tzero) ...
  static constexpr double fsatb  = -0.002;  // [1/K] ... (fsatb tuned to Antarctic d-excess)
  static constexpr double tzero  = 273.16;  // [K] SHR_CONST_TKTRIP
  static constexpr double ssatmx = 2.0;     // largest supersaturation ratio over ice

private:
  // akel's denominator alpha_eq*(D/D_i)^dkfac*(h - 1) + 1 is positive while
  // alpha_eq*(D/D_i)^dkfac < 1/(1 - h), which is >= 4 since h >= 0.75: for the
  // liquid fits, at any temperature above ~160 K (HTO, the largest alpha; ~130
  // K for HDO). Only a fit extrapolated far below any liquid water gets there;
  // this floor keeps akel finite and positive in that case (the replay bounds
  // the resulting beta anyway).
  static constexpr double akel_denominator_min = 1.e-6;

public:
  // -----------------------------------------------------------------------
  // wiso_heff: the effective humidity of evaporation into air of relative
  // humidity rh over liquid [-], min(1, fkhum*min(1, rh) + 1 - fkhum), in
  // [0.75, 1]. rh < 0 (only from a negative vapor) counts as 0; the Fortran
  // does not clip it.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT effective_humidity(const ScalarT& rh)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    // Local copies: device code may not bind a reference to a class-scope
    // constant (the Pack operators take their scalars by reference)
    constexpr RealT fk = fkhum, zero = 0, one = 1;
    const ScalarT h0 = ekat::min(ekat::max(rh, zero), one);
    return ekat::min(fk*h0 + one - fk, one);
  }

  // -----------------------------------------------------------------------
  // wiso_ssatf: the empirical supersaturation ratio over ice of temperature T
  // [K], min(max(fsata + fsatb*(T - tzero), 1), ssatmx) [-]. It rises from 1
  // at tzero by 0.002 per K of cooling (1.04 at 253.15 K, 1.15 at 200 K).
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT ice_supersaturation_empirical(const ScalarT& T)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    constexpr RealT a = fsata, b = fsatb, t0 = tzero, smax = ssatmx, one = 1;
    return ekat::min(ekat::max(a + b*(T - t0), one), smax);
  }

  // -----------------------------------------------------------------------
  // wiso_akel: the effective liquid/vapor factor of liquid evaporating into
  // air of relative humidity rh over liquid [-] (see above), given the
  // equilibrium factor alpha_eq (CondensedOverVapor) of the same species.
  // = alpha_eq exactly where the effective humidity is 1 (rh >= 1).
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_liquid_vapor_kinetic(
      const ScalarT& alpha_eq,
      const ScalarT& rh,
      const WaterIsotopologues species,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& mask = ekat::Mask<ScalarT::n>(true))
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    constexpr RealT one = 1, kfac = dkfac, den_min = akel_denominator_min;

    if (species == WaterIsotopologues::H216O || !mask.any()) {
      return alpha_eq;
    }
    ScalarT rh_live(one);
    rh_live.set(mask, rh);
    const ScalarT h = effective_humidity(rh_live);
    const auto kinetic = mask && h < one;
    if (!kinetic.any()) {
      return alpha_eq;
    }

    ScalarT a(one);
    a.set(kinetic, alpha_eq);

    // (D/D_i)^dkfac of the species ("dondi" in the Fortran)
    const RealT dondi = Kokkos::pow(one/iso_params.difrm(species), kfac);

    const ScalarT den = ekat::max(a*dondi*(h - one) + one, den_min);
    ScalarT akel = alpha_eq;
    akel.set(kinetic, a*h/den);
    return akel;
  }

  // -----------------------------------------------------------------------
  // wiso_akci: the effective ice/vapor factor of vapor depositing on ice at
  // temperature T [K] and supersaturation ratio S over ice [-] (clipped to
  // [1, ssatmx]; see above), given the equilibrium factor alpha_eq
  // (CondensedOverVapor) of the same species. = alpha_eq exactly at
  // T >= tkini and where S <= 1.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor_kinetic(
      const ScalarT& alpha_eq,
      const ScalarT& T,
      const ScalarT& S,
      const WaterIsotopologues species,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& mask = ekat::Mask<ScalarT::n>(true))
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    constexpr RealT one = 1, t_kin = tkini, smax = ssatmx;

    if (species == WaterIsotopologues::H216O || !mask.any()) {
      return alpha_eq;
    }
    // (T is compared on every lane, so padding gets a placeholder first: an
    // ordered comparison with NaN raises FE_INVALID)
    ScalarT T_live(t_kin), S_live(one);
    T_live.set(mask, T);
    S_live.set(mask && T_live < t_kin, S);
    const ScalarT sat = ekat::min(ekat::max(S_live, one), smax);
    const auto kinetic = mask && sat > one;
    if (!kinetic.any()) {
      return alpha_eq;
    }

    ScalarT a(one);
    a.set(kinetic, alpha_eq);

    // D/D_i of the species ("dondi" in the Fortran)
    const RealT dondi = one/iso_params.difrm(species);

    ScalarT akci = alpha_eq;
    akci.set(kinetic, a*sat/(a*dondi*(sat - one) + one));
    return akci;
  }
};

} // namespace water_isotopes
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_KINETIC_HPP
