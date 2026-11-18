#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_fpe.hpp>
#include <ekat_pack.hpp>

#include <cfenv>
#include <cmath>
#include <limits>

namespace scream {
namespace {

namespace wt = water_tracers;
using Pack1 = ekat::Pack<Real,1>;

// An evaporating ocean column (q < s)
const Real f0 = 0.9, rho0 = 1.2, u0 = 0.3, r0 = 1.5e-3, s0 = 0.02, q0 = 0.01;

Pack1 E (const Real f, const Real q) {
  return wt::bulk_ocean_flux(Pack1(f), Pack1(rho0), Pack1(u0), Pack1(r0), Pack1(s0), Pack1(q));
}

Pack1 flux_of (const wt::OceanFlux kind, const Real r_ocn, const Real q, const Real q_i,
              const Real alpha = 1) {
  return wt::ocean_flux(kind, r_ocn, Pack1(f0), Pack1(rho0), Pack1(u0), Pack1(r0),
                        Pack1(s0), Pack1(q), Pack1(q_i), E(f0, q), Pack1(alpha));
}

TEST_CASE("ocean_flux") {
  using OF = wt::OceanFlux;

  SECTION("bulk_formula with R_ocn=1, q_i=q equals E_ocn bit-for-bit") {
    REQUIRE(flux_of(OF::BulkFormula, 1, q0, q0)[0] == E(f0, q0)[0]);
  }
  SECTION("bulk_proportional is R_ocn*E_ocn") {
    REQUIRE(flux_of(OF::BulkProportional, 0.5, q0, 3e-3)[0] == Real(0.5)*E(f0, q0)[0]);
    REQUIRE(flux_of(OF::BulkProportional, 1,   q0, 3e-3)[0] == E(f0, q0)[0]);
  }
  SECTION("bulk_formula uses s*R_ocn and q_i") {
    const Real r_ocn = 0.5, q_i = 2e-3;
    const Real ref = wt::bulk_ocean_flux(Pack1(f0), Pack1(rho0), Pack1(u0), Pack1(r0),
                                         Pack1(s0*r_ocn), Pack1(q_i))[0];
    REQUIRE(flux_of(OF::BulkFormula, r_ocn, q0, q_i)[0] == ref);
  }
  SECTION("dew removes the tracer at its ratio to qv") {
    const Real q = 0.03, q_i = 0.012;  // q > s: dew
    const Real e = E(f0, q)[0];
    REQUIRE(e < 0);
    const Real ref = (wt::tracer_ratio(Pack1(q_i), Pack1(q))*e)[0];
    REQUIRE(flux_of(OF::BulkProportional, 0.5, q, q_i)[0] == ref);
    REQUIRE(flux_of(OF::BulkFormula,      0.5, q, q_i)[0] == ref);
  }
  SECTION("bulk_formula with alpha != 1 scales s*r_ocn/alpha") {
    const Real r_ocn = 0.7, q_i = 2.5e-3, alpha = 1.02;
    const Real ref = wt::bulk_ocean_flux(Pack1(f0), Pack1(rho0), Pack1(u0), Pack1(r0),
                                         Pack1(s0*r_ocn/alpha), Pack1(q_i))[0];
    REQUIRE(flux_of(OF::BulkFormula, r_ocn, q0, q_i, alpha)[0] == ref);
  }
  SECTION("dew ignores alpha") {
    const Real q = 0.03, q_i = 0.012;  // q > s: dew
    const Real with_alpha    = flux_of(OF::BulkFormula, 0.5, q, q_i, Real(1.3))[0];
    const Real without_alpha = flux_of(OF::BulkFormula, 0.5, q, q_i)[0];
    REQUIRE(with_alpha == without_alpha);
  }
}

TEST_CASE("invalid_ocean") {
  // Trap explicitly: no step may compute with fill values or NaN.
  const int saved = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);

  const Real nan  = std::numeric_limits<Real>::quiet_NaN();
  const Real fill = 1e30;
  struct State { Real f, u, r, s; };
  const State states[] = {
    {0,  u0, r0,   s0},    // no ocean
    {f0, 0,  r0,   s0},    // no friction velocity
    {f0, u0, r0,   fill},  // coupler fill value
    {f0, u0, fill, s0},
    {f0, u0, nan,  s0},    // NaN
    {f0, u0, r0,   nan},
    {nan, u0, r0,  s0},
  };

  for (const auto& st : states) {
    Pack1 f(st.f), u(st.u), r(st.r), s(st.s);
    REQUIRE(not wt::ocean_valid(f, u, r, s)[0]);

    const Pack1 e = wt::ocean_evaporation(f, u, r, s, Pack1(rho0), Pack1(q0));
    REQUIRE(e[0] == 0);
    REQUIRE(f[0] == 0);
    REQUIRE(u[0] == 0);
    REQUIRE(r[0] == 0);
    REQUIRE(s[0] == 0);

    for (const auto kind : {wt::OceanFlux::BulkProportional, wt::OceanFlux::BulkFormula}) {
      REQUIRE(wt::ocean_flux(kind, Real(1), f, Pack1(rho0), u, r, s, Pack1(q0), Pack1(q0), e)[0] == 0);
      // A garbage alpha must never be evaluated off the evap branch.
      REQUIRE(wt::ocean_flux(kind, Real(1), f, Pack1(rho0), u, r, s, Pack1(q0), Pack1(q0), e,
                             Pack1(nan))[0] == 0);
    }
  }

  // A valid state is left untouched, and E_ocn is the bulk formula exactly
  Pack1 f(f0), u(u0), r(r0), s(s0);
  REQUIRE(wt::ocean_valid(f, u, r, s)[0]);
  const Pack1 e = wt::ocean_evaporation(f, u, r, s, Pack1(rho0), Pack1(q0));
  REQUIRE(f[0] == f0);
  REQUIRE(e[0] == E(f0, q0)[0]);

  ekat::disable_all_fpes();
  ekat::enable_fpes(saved);
}

// Total surface flux (ocean + land/ice part) for a given merged evaporation s_evap
Pack1 surface_flux_of (const wt::OceanFlux kind, const Real r_ocn, const Real r_land,
                       const Real q, const Real q_i, const Real s_evap, const Real f = f0,
                       const Real alpha = 1) {
  return wt::surface_tracer_flux(kind, r_ocn, r_land, Pack1(f), Pack1(rho0), Pack1(u0), Pack1(r0),
                                 Pack1(s0), Pack1(q), Pack1(q_i), Pack1(s_evap), E(f, q), Pack1(alpha));
}

TEST_CASE("surface_tracer_flux") {
  using OF = wt::OceanFlux;
  const Real rtol = 1e-14;

  // Trap explicitly: the masked lanes must never divide or compare garbage.
  const int saved = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);

  const Real e_ev = E(f0, q0)[0];                    // ocean evaporates: 4.86e-6
  const Real q_dew = 0.03, q_i_dew = 0.012;
  const Real e_dw = E(f0, q_dew)[0];                 // ocean dew: -4.86e-6
  REQUIRE(e_ev > 0);
  REQUIRE(e_dw < 0);

  SECTION("R_ocn = R_land = 1, everything evaporates: exactly s_evap") {
    for (const Real s_evap : {2e-5, e_ev, 1e-5, 4.9e-6}) {
      REQUIRE(surface_flux_of(OF::BulkProportional, 1, 1, q0, 3e-3, s_evap)[0] == s_evap);
    }
    // all land (f=0): the same
    REQUIRE(surface_flux_of(OF::BulkProportional, 1, 1, q0, 3e-3, Real(2e-5), Real(0))[0] == Real(2e-5));
  }
  SECTION("R_land != R_ocn, both evaporate") {
    const Real s_evap = 2e-5, r_ocn = 0.5, r_land = 2;
    const Real flux = surface_flux_of(OF::BulkProportional, r_ocn, r_land, q0, 3e-3, s_evap)[0];
    REQUIRE(flux == Approx(r_ocn*e_ev + r_land*(s_evap - e_ev)).epsilon(rtol));
    REQUIRE(flux == r_land*s_evap + (r_ocn - r_land)*e_ev);
  }
  SECTION("ocean evaporates, land dews") {
    const Real s_evap = 2e-6, r_ocn = 0.5, r_land = 2;  // s_evap < e_ev: L < 0
    const Real q_i = 3e-3;
    const Real flux = surface_flux_of(OF::BulkProportional, r_ocn, r_land, q0, q_i, s_evap)[0];
    REQUIRE(flux == Approx(r_ocn*e_ev + (q_i/q0)*(s_evap - e_ev)).epsilon(rtol));
    // r_land is not used for dew
    REQUIRE(flux == surface_flux_of(OF::BulkProportional, r_ocn, Real(7), q0, q_i, s_evap)[0]);
  }
  SECTION("ocean dews, land evaporates") {
    const Real s_evap = 1e-6, r_ocn = 0.5, r_land = 2;  // L = s_evap - e_dw > 0
    const Real flux = surface_flux_of(OF::BulkProportional, r_ocn, r_land, q_dew, q_i_dew, s_evap)[0];
    REQUIRE(flux == Approx((q_i_dew/q_dew)*e_dw + r_land*(s_evap - e_dw)).epsilon(rtol));
    // r_ocn is not used for dew
    REQUIRE(flux == surface_flux_of(OF::BulkProportional, Real(7), r_land, q_dew, q_i_dew, s_evap)[0]);
  }
  SECTION("both dew") {
    const Real s_evap = 2*e_dw;
    const Real flux = surface_flux_of(OF::BulkProportional, 0.5, 2, q_dew, q_i_dew, s_evap)[0];
    REQUIRE(flux == Approx((q_i_dew/q_dew)*s_evap).epsilon(rtol));
    // A tracer with the bulk ratio removes exactly the bulk dew
    REQUIRE(surface_flux_of(OF::BulkProportional, 1, 1, q_dew, q_dew, s_evap)[0] == Approx(s_evap).epsilon(rtol));
    REQUIRE(surface_flux_of(OF::BulkProportional, 1, 1, q_dew, q_dew, 1e-6)[0] == Approx(1e-6).epsilon(rtol));
    REQUIRE(surface_flux_of(OF::BulkProportional, 1, 1, q0, q0, 2e-6)[0] == Approx(2e-6).epsilon(rtol));
  }
  SECTION("f = 0 (all land and ice)") {
    const Real f = 0;  // E(0, q) = 0
    REQUIRE(surface_flux_of(OF::BulkProportional, 0.5, 2, q0, 3e-3, Real(2e-5), f)[0] == Real(2)*Real(2e-5));
    REQUIRE(surface_flux_of(OF::BulkFormula,      0.5, 2, q0, 3e-3, Real(2e-5), f)[0] == Real(2)*Real(2e-5));
    // dew over land
    REQUIRE(surface_flux_of(OF::BulkProportional, 0.5, 2, q0, 3e-3, Real(-2e-6), f)[0]
            == Approx((3e-3/q0)*Real(-2e-6)).epsilon(rtol));
    REQUIRE(surface_flux_of(OF::BulkProportional, 0.5, 2, q0, 3e-3, Real(0), f)[0] == 0);
  }
  SECTION("bulk_formula: ocean bulk formula plus the proportional land part") {
    const Real s_evap = 2e-5, r_ocn = 0.5, r_land = 2, q_i = 2e-3;
    const Real flux = surface_flux_of(OF::BulkFormula, r_ocn, r_land, q0, q_i, s_evap)[0];
    const Real f_ocn = flux_of(OF::BulkFormula, r_ocn, q0, q_i)[0];
    REQUIRE(flux == Approx(f_ocn + r_land*(s_evap - e_ev)).epsilon(rtol));
    // alpha acts on the ocean part only
    const Real a = surface_flux_of(OF::BulkFormula, r_ocn, r_land, q0, q_i, s_evap, f0, Real(1.02))[0];
    REQUIRE(a == Approx(flux_of(OF::BulkFormula, r_ocn, q0, q_i, Real(1.02))[0] +
                        r_land*(s_evap - e_ev)).epsilon(rtol));
    // dew on land and ocean evaporation
    const Real b = surface_flux_of(OF::BulkFormula, r_ocn, r_land, q0, q_i, Real(2e-6))[0];
    REQUIRE(b == Approx(f_ocn + (q_i/q0)*(Real(2e-6) - e_ev)).epsilon(rtol));
  }
  SECTION("invalid ocean state: only the land part, nothing computed with fill values") {
    const Real nan  = std::numeric_limits<Real>::quiet_NaN();
    const Real fill = 1e30;
    struct State { Real f, u, r, s; };
    const State states[] = {
      {0,  u0, r0,   s0}, {f0, 0, r0, s0}, {f0, u0, r0, fill}, {f0, u0, fill, s0},
      {f0, u0, nan, s0}, {f0, u0, r0, nan}, {nan, u0, r0, s0},
    };
    for (const auto& st : states) {
      Pack1 f(st.f), u(st.u), r(st.r), sv(st.s);
      const Pack1 e = wt::ocean_evaporation(f, u, r, sv, Pack1(rho0), Pack1(q0));
      REQUIRE(e[0] == 0);
      for (const auto kind : {OF::BulkProportional, OF::BulkFormula}) {
        for (const Real s_evap : {2e-5, -3e-6}) {
          const Real q_i = 3e-3;
          const Real flux = wt::surface_tracer_flux(kind, Real(0.5), Real(2), f, Pack1(rho0), u, r, sv,
                                                    Pack1(q0), Pack1(q_i), Pack1(s_evap), e,
                                                    Pack1(nan))[0];  // garbage alpha: never used here
          const Real ref = s_evap>0 ? Real(2)*s_evap : (q_i/q0)*s_evap;
          REQUIRE(flux == Approx(ref).epsilon(rtol));
        }
      }
    }
  }
  SECTION("dry bottom level (q = q_i = 0) does not form 0/0") {
    const Real flux = surface_flux_of(OF::BulkProportional, 1, 1, Real(0), Real(0), Real(-1e-6), Real(0))[0];
    REQUIRE(std::isfinite(flux));
  }

  ekat::disable_all_fpes();
  ekat::enable_fpes(saved);
}

TEST_CASE("is_evaporating") {
  // Evaporating: valid ocean, E_ocn > 0.
  REQUIRE(wt::is_evaporating(Pack1(f0), E(f0, q0))[0]);
  // Dew: valid ocean, E_ocn <= 0.
  REQUIRE(not wt::is_evaporating(Pack1(f0), E(f0, Real(0.03)))[0]);
  // Not a valid ocean point, even with a positive E_ocn value.
  REQUIRE(not wt::is_evaporating(Pack1(0), Pack1(1e-5))[0]);
}

} // anonymous namespace
} // namespace scream
