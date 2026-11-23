#ifndef EAMXX_WATER_ISOTOPES_KINETIC_HPP
#define EAMXX_WATER_ISOTOPES_KINETIC_HPP

#include "share/core/eamxx_types.hpp"
#include "share/physics/physics_constants.hpp"
#include "eamxx_water_isotopes_parameters.hpp"
#include "eamxx_water_isotopes_fractionation.hpp"  // H217O/HTO mass-scaling exponents

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::log/ekat::pow overloads for ekat::Pack

namespace scream {
namespace water_isotopes {

/*
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
 * H217O and HTO are not in the Fortran. Their kmol is derived from H218O/HDO
 * with the same mass-scaling exponents as the equilibrium factor:
 * kmol_H217O = 1 - (1 - kmol_H218O)^0.529, kmol_HTO = 1 - (1 - kmol_HDO)^2.
 * H216O has alpkn = 1, as does every species under
 * KineticFractionationFormulation::None.
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
    using WIF   = WaterIsotopeFractionation;

    if (iso_params.kinetic_formulation() == KineticFractionationFormulation::None ||
        species == WaterIsotopologues::H216O || !mask.any()) {
      return ScalarT(1);
    }

    ScalarT rbot_live{RealT(rbot_lane_fill)}, zbot_live{RealT(zbot_lane_fill)},
            ustar_live{RealT(ustar_lane_fill)};
    rbot_live.set(mask, rbot);
    zbot_live.set(mask, zbot);
    ustar_live.set(mask, ustar);

    const IsoElement el = WaterIsotopeParameters<RealT>::element_of(species);
    const ScalarT k = kmol(rbot_live, zbot_live, ustar_live, iso_params.difrm(el));

    // alpkn of the parent species (HDO or H218O). For H217O/HTO,
    // kmol = 1 - (1 - kmol_parent)^exponent is alpkn = alpkn_parent^exponent.
    const ScalarT alpkn_parent = RealT(1) - k;

    ScalarT alpkn(1);
    switch (species) {
      case WaterIsotopologues::H217O:
        alpkn.set(mask, pow(alpkn_parent, RealT(WIF::H217O_exponent)));
        break;
      case WaterIsotopologues::HTO:
        alpkn.set(mask, pow(alpkn_parent, RealT(WIF::HTO_exponent)));
        break;
      default:  // HDO, H218O (H216O returned above)
        alpkn.set(mask, alpkn_parent);
        break;
    }
    return alpkn;
  }
};

} // namespace water_isotopes
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_KINETIC_HPP
