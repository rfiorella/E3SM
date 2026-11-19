#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracer_transfers_functions.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_kinetic.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_fpe.hpp>
#include <ekat_kokkos_types.hpp>
#include <ekat_pack.hpp>

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>

namespace scream {
namespace {

namespace wt = water_tracers;
using Pack1 = ekat::Pack<Real,1>;
using PackN = ekat::Pack<Real,SCREAM_PACK_SIZE>;

// Relative tolerance of the roundoff checks
constexpr Real tol = std::is_same<Real,double>::value ? 1e-13 : 1e-5;

// Traps FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW for its lifetime (as the fpe
// build does), so the degenerate-input tests prove that nothing raises them
struct FpeTrap {
  int saved;
  FpeTrap () : saved(ekat::get_enabled_fpes()) {
    ekat::disable_all_fpes();
    ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  }
  ~FpeTrap () {
    ekat::disable_all_fpes();
    ekat::enable_fpes(saved);
  }
};

bool rel_close (const Real a, const Real b, const Real rtol = tol) {
  return std::abs(a - b) <= rtol*std::max(std::abs(a), std::abs(b));
}

// donor_transfers on one donor, scalars in and out (a pack of one)
template <int N>
bool transfer (const Real tr, const Real q, const Real (&amount)[N], const Real (&beta)[N],
               Real (&moved)[N])
{
  Pack1 a[N], b[N], m[N];
  for (int k=0; k<N; ++k) {
    a[k] = amount[k];
    b[k] = beta[k];
  }
  const auto rayleigh = wt::donor_transfers(Pack1(tr), Pack1(q), a, b, m);
  for (int k=0; k<N; ++k) {
    moved[k] = m[k][0];
  }
  return rayleigh[0];
}

template <int N>
Real sum (const Real (&x)[N]) {
  Real s = 0;
  for (int k=0; k<N; ++k) s += x[k];
  return s;
}

} // anonymous namespace

TEST_CASE("donor_transfers_linear") {
  const Real tr = 3.7e-3, q = 7.3e-3;

  SECTION("every beta 1: the linear transfer, bit for bit") {
    const Real amount[3] = {1e-4, 2.5e-6, 0}, beta[3] = {1, 1, 1};
    Real moved[3];
    REQUIRE(not transfer(tr, q, amount, beta, moved));
    const Real ratio = wt::tracer_ratio(Pack1(tr), Pack1(q))[0];
    for (int k=0; k<3; ++k) {
      REQUIRE(moved[k] == amount[k]*ratio);
    }
    // including a dry donor (tracer_ratio's r_std) and an overdrawn one
    const Real big[1] = {2*q}, one[1] = {1};
    Real m1[1];
    REQUIRE(not transfer(tr, q, big, one, m1));
    REQUIRE(m1[0] == big[0]*ratio);
    REQUIRE(not transfer(0, 0, big, one, m1));
    REQUIRE(m1[0] == big[0]);
  }

  SECTION("the Rayleigh path at beta_eff = 1 is the linear transfer to roundoff") {
    // An idle sink with beta != 1 sends the donor through the Rayleigh path,
    // with beta_eff = 1
    for (const Real x : {1e-12, 1e-6, 1e-2, 0.5, 0.99}) {
      INFO("A/q = " << x);
      const Real amount[2] = {x*q, 0}, beta[2] = {1, 1.2};
      Real moved[2];
      REQUIRE(transfer(tr, q, amount, beta, moved));
      REQUIRE(rel_close(moved[0], amount[0]*(tr/q)));
      REQUIRE(moved[1] == 0);
    }
  }
}

TEST_CASE("donor_transfers_rayleigh") {
  const Real tr = 2.0e-3, q = 1.0e-2;

  SECTION("one sink: the closed-system Rayleigh law R/R0 = f^(beta-1)") {
    for (const Real beta : {1.1, 0.9, 1.02, 1/1.1, 1.5, 0.7}) {
      for (const Real f : {1 - 1e-9, 0.999, 0.9, 0.5, 0.1, 1e-3}) {
        INFO("beta = " << beta << ", f = " << f);
        const Real amount[1] = {(1 - f)*q}, b[1] = {beta};
        Real moved[1];
        REQUIRE(transfer(tr, q, amount, b, moved));
        // what is left, and its ratio
        const Real tr_left = tr - moved[0], q_left = q - amount[0];
        const Real R0 = tr/q;
        REQUIRE(rel_close(tr_left/q_left, R0*std::pow(q_left/q, beta - 1), 1e3*tol));
        // the first, small part removed has R0*beta
        if (f > 0.999) {
          REQUIRE(rel_close(moved[0]/amount[0], R0*beta, 1e-2*std::abs(beta - 1) + tol));
        }
      }
    }
  }

  SECTION("positivity as f -> 0, for beta above and below 1") {
    for (const Real beta : {0.5, 0.9, 1/1.1, 1.1, 2.0}) {
      for (const Real eps : {1e-3, 1e-8, 1e-14, 0.0}) {
        INFO("beta = " << beta << ", 1 - A/q = " << eps);
        const Real amount[1] = {q*(1 - eps)}, b[1] = {beta};
        Real moved[1];
        transfer(tr, q, amount, b, moved);
        REQUIRE(moved[0] >= 0);
        REQUIRE(moved[0] <= tr);
        REQUIRE(tr - moved[0] >= 0);
      }
      // the sinks take all of the bulk water, or more: all of the tracer leaves
      for (const Real over : {1.0, 1.5}) {
        const Real amount[2] = {0.5*over*q, 0.5*over*q}, b[2] = {beta, 1};
        Real moved[2];
        transfer(tr, q, amount, b, moved);
        REQUIRE(tr - (moved[0] + moved[1]) >= 0);
        REQUIRE(rel_close(moved[0] + moved[1], tr));
      }
    }
  }

  SECTION("split in proportion to amount*beta") {
    const Real amount[3] = {1e-4, 3e-4, 2e-5}, beta[3] = {1.1, 0.95, 1};
    Real moved[3];
    REQUIRE(transfer(tr, q, amount, beta, moved));
    const Real total = sum(moved);
    Real w = 0;
    for (int k=0; k<3; ++k) w += amount[k]*beta[k];
    for (int k=0; k<3; ++k) {
      REQUIRE(rel_close(moved[k], total*amount[k]*beta[k]/w));
    }
    // and the total follows beta_eff
    const Real A = sum(amount);
    REQUIRE(rel_close(total, tr*(1 - std::pow(1 - A/q, w/A))));
  }

  SECTION("independent of the order of the sinks") {
    Real amount[3] = {1e-4, 3e-4, 2e-5}, beta[3] = {1.1, 0.95, 1};
    Real ref[3];
    transfer(tr, q, amount, beta, ref);
    int perm[3] = {0, 1, 2};
    while (std::next_permutation(perm, perm + 3)) {
      Real a[3], b[3], m[3];
      for (int k=0; k<3; ++k) {
        a[k] = amount[perm[k]];
        b[k] = beta[perm[k]];
      }
      transfer(tr, q, a, b, m);
      for (int k=0; k<3; ++k) {
        REQUIRE(rel_close(m[k], ref[perm[k]]));
      }
    }
  }

  SECTION("invariant to splitting the step") {
    const Real amount[2] = {2e-3, 3e-3}, beta[2] = {1.1, 0.9};
    Real one_step[2];
    transfer(tr, q, amount, beta, one_step);
    for (const int nsub : {2, 5}) {
      INFO("substeps: " << nsub);
      Real t = tr, m = q, moved_sum[2] = {0, 0};
      for (int s=0; s<nsub; ++s) {
        const Real a[2] = {amount[0]/nsub, amount[1]/nsub};
        Real d[2];
        transfer(t, m, a, beta, d);
        t -= d[0] + d[1];
        m -= a[0] + a[1];
        moved_sum[0] += d[0];
        moved_sum[1] += d[1];
      }
      REQUIRE(rel_close(t, tr - sum(one_step), 10*tol));
      // per sink the split differs at second order in the ratio's change
      // over the step (the sinks see the same ratios either way only in sum)
      REQUIRE(rel_close(moved_sum[0] + moved_sum[1], sum(one_step), 10*tol));
    }
  }

  SECTION("beta is kept within [beta_min, beta_max]") {
    const Real amount[1] = {1e-3};
    Real moved[1];
    const Real zero[1] = {0};
    transfer(tr, q, amount, zero, moved);
    REQUIRE(moved[0] > 0);
    REQUIRE(rel_close(moved[0], tr*(1 - std::pow(1 - amount[0]/q, Real(wt::fractionation_beta_min))), 1e-6));
    const Real huge[1] = {1e30};
    transfer(tr, q, amount, huge, moved);
    REQUIRE(moved[0] == tr);  // f^1e6 underflows to 0
  }
}

TEST_CASE("donor_transfers_fpe_safety") {
  FpeTrap trap;

  // Degenerate donors, one per lane of a pack (and one per pack of one):
  // tr, q, two sinks' amounts and betas, and the expected total moved
  struct Case { Real tr, q, a0, a1, b0, b1, total; };
  const Case cases[] = {
    {0,      0,      0,    0,    1.1, 1,   0},       // empty donor, no sink
    {1e-3,   5e-3,   0,    0,    1.1, 0.9, 0},       // no sink
    {1e-3,   0,      1e-4, 0,    1.1, 1,   1e-3},    // sink from a dry donor: all leaves
    {1e-3,  -1e-3,   1e-4, 0,    1.1, 1,   1e-3},    // negative bulk water
    {1e-3,   1e-30,  1e-4, 0,    1.1, 1,   1e-3},    // below water_tracers_qmin
    {0,      5e-3,   1e-3, 0,    1.1, 1,   0},       // no tracer
    {-1e-20, 5e-3,   5e-3, 0,    1.1, 1,  -1e-20},   // negative tracer, all taken
    {1e-3,   5e-3,   5e-3, 0,    0.9, 1,   1e-3},    // A = q exactly
    {1e-3,   5e-3,   4e-3, 4e-3, 1.1, 2,   1e-3},    // A > q
    {1e-3,   5e-3,   1e-3, 1e-3, 0,   1e30, -1},     // betas out of bounds (-1: just finite, in [0, tr])
  };

  for (const auto& c : cases) {
    INFO("tr = " << c.tr << ", q = " << c.q << ", amounts = " << c.a0 << ", " << c.a1);
    const Real amount[2] = {c.a0, c.a1}, beta[2] = {c.b0, c.b1};
    Real moved[2];
    transfer(c.tr, c.q, amount, beta, moved);
    REQUIRE(std::isfinite(moved[0]));
    REQUIRE(std::isfinite(moved[1]));
    const Real total = moved[0] + moved[1];
    if (c.total == -1) {
      REQUIRE(total >= 0);
      REQUIRE(total <= c.tr);
    } else {
      REQUIRE(total == c.total);
    }
  }

  // The same cases in the lanes of packs: each lane as in a pack of one
  constexpr int ncases = sizeof(cases)/sizeof(Case);
  for (int start=0; start<ncases; start+=PackN::n) {
    PackN tr, q, a[2], b[2], m[2];
    for (int s=0; s<PackN::n; ++s) {
      const auto& c = cases[std::min(start + s, ncases - 1)];
      tr[s] = c.tr; q[s] = c.q;
      a[0][s] = c.a0; a[1][s] = c.a1;
      b[0][s] = c.b0; b[1][s] = c.b1;
    }
    wt::donor_transfers(tr, q, a, b, m);
    for (int s=0; s<PackN::n; ++s) {
      const auto& c = cases[std::min(start + s, ncases - 1)];
      const Real amount[2] = {c.a0, c.a1}, beta[2] = {c.b0, c.b1};
      Real moved[2];
      transfer(c.tr, c.q, amount, beta, moved);
      REQUIRE(m[0][s] == moved[0]);
      REQUIRE(m[1][s] == moved[1]);
    }
  }
}

TEST_CASE("shoc_condensation") {
  const Real tr_v = 6e-3, q_v = 1.2e-2, tr_c = 1e-4, q_c = 2e-4;

  SECTION("beta = 1: the linear transfer, bit for bit, with the sign of the amount") {
    for (const Real cond : {3e-5, -5e-5}) {
      const Real tr_d = cond > 0 ? tr_v : tr_c, q_d = cond > 0 ? q_v : q_c;
      const Pack1 moved = wt::shoc_condensation(Pack1(tr_d), Pack1(q_d), Pack1(cond), Pack1(1));
      REQUIRE(moved[0] == cond*wt::tracer_ratio(Pack1(tr_d), Pack1(q_d))[0]);
    }
  }

  SECTION("beta: alpha_l for condensation, 1/alpha_l for evaporation") {
    const Real alpha = 1.08;
    REQUIRE(wt::shoc_condensation_beta(Pack1(alpha), Pack1(1e-5))[0] == alpha);
    REQUIRE(wt::shoc_condensation_beta(Pack1(alpha), Pack1(-1e-5))[0] == 1/alpha);
    // all of the cloud evaporates: all of its tracer goes to vapor
    const Pack1 moved = wt::shoc_condensation(Pack1(tr_c), Pack1(q_c), Pack1(-q_c),
                                              wt::shoc_condensation_beta(Pack1(alpha), Pack1(-q_c)));
    REQUIRE(moved[0] == -tr_c);
  }

  SECTION("kinetic: condensation keeps alpha_l, evaporation takes 1/akel") {
    using WIK = water_isotopes::WaterIsotopeKinetic;
    const water_isotopes::WaterIsotopeParameters<Real> iso_params;
    const Real alpha = 1.08;
    for (const Real rh : {0.5, 0.9, 1.0}) {
      const Pack1 akel = WIK::alpha_liquid_vapor_kinetic(Pack1(alpha), Pack1(rh),
                                                         water_isotopes::WaterIsotopologues::HDO, iso_params);
      INFO("rh = " << rh << ", akel = " << akel[0]);
      REQUIRE(wt::shoc_condensation_beta(Pack1(alpha), akel, Pack1(1e-5))[0] == alpha);
      REQUIRE(wt::shoc_condensation_beta(Pack1(alpha), akel, Pack1(-1e-5))[0] == 1/akel[0]);
      // the evaporate is at least as depleted as at equilibrium
      REQUIRE(1/akel[0] <= 1/alpha);
      if (rh < 1) REQUIRE(1/akel[0] < 1/alpha);
    }
    // The two-factor form with alpha_evap = alpha_l is the equilibrium one, bit for bit
    for (const Real cond : {3e-5, -5e-5, 0.0}) {
      REQUIRE(wt::shoc_condensation_beta(Pack1(alpha), Pack1(alpha), Pack1(cond))[0] ==
              wt::shoc_condensation_beta(Pack1(alpha), Pack1(cond))[0]);
    }
  }

  SECTION("new cloud from vapor: R_qc/R_qv = alpha_l(T), to first order in the condensed fraction") {
    using WIF = water_isotopes::WaterIsotopeFractionation;
    const water_isotopes::WaterIsotopeParameters<Real> iso_params;
    for (const auto species : {water_isotopes::WaterIsotopologues::HDO, water_isotopes::WaterIsotopologues::H218O}) {
      for (const Real T : {300.0, 280.0, 260.0}) {
        const Real alpha = WIF::alpha_liquid_vapor(Pack1(T), species, water_isotopes::CondensedOverVapor,
                                                   iso_params)[0];
        const Real x = 1e-2;  // condensed fraction of the vapor
        const Real cond = x*q_v;
        const Pack1 moved = wt::shoc_condensation(Pack1(tr_v), Pack1(q_v), Pack1(cond),
                                                  wt::shoc_condensation_beta(Pack1(alpha), Pack1(cond)));
        const Real R_c = moved[0]/cond, R_v = (tr_v - moved[0])/(q_v - cond);
        INFO("T = " << T << ", alpha = " << alpha << ", R_c/R_v = " << R_c/R_v);
        REQUIRE(R_c/R_v > 1);  // the heavy isotope goes to the condensate
        REQUIRE(std::abs(R_c/R_v/alpha - 1) < x*(alpha - 1));
      }
    }
  }
}

TEST_CASE("part2_beta") {
  using WT = p3::WaterTransfer;
  const Pack1 al(1.09), ai(1.13);
  for (int cmp=WT::Part2VapdepV2I; cmp<=WT::Part2MeltI2R; ++cmp) {
    INFO("component " << cmp);
    const Real b = wt::part2_beta(cmp, al, ai)[0];
    if (cmp==WT::Part2VapdepV2I or cmp==WT::Part2NucleatV2I) {
      REQUIRE(b == ai[0]);
    } else if (cmp==WT::Part2EvapR2V) {
      REQUIRE(b == 1/al[0]);
    } else if (cmp==WT::Part2BergC2I) {
      REQUIRE(b == ai[0]/al[0]);
    } else {
      REQUIRE(b == 1);
    }
    // a non-fractionating tracer, or alpha_override (Bergeron: alpha/alpha)
    REQUIRE(wt::part2_beta(cmp, Pack1(1), Pack1(1))[0] == 1);
    if (cmp==WT::Part2BergC2I) {
      REQUIRE(wt::part2_beta(cmp, Pack1(1.37), Pack1(1.37))[0] == 1);
    }
  }

  SECTION("Bergeron goes through the vapor") {
    // A small Bergeron transfer gives the ice the liquid's ratio times
    // alpha_i/alpha_l: what evaporating the liquid (1/alpha_l) and depositing
    // that vapor (alpha_i) would give
    const Real tr_c = 1e-4, q_c = 2e-4, a = 1e-9;
    const Pack1 amount[1] = {Pack1(a)};
    const Pack1 beta[1] = {wt::part2_beta(WT::Part2BergC2I, al, ai)};
    Pack1 moved[1];
    wt::donor_transfers(Pack1(tr_c), Pack1(q_c), amount, beta, moved);
    const Real R_ice_new = moved[0][0]/a, R_c = tr_c/q_c;
    REQUIRE(rel_close(R_ice_new, R_c*(1/al[0])*ai[0], 1e-4));
  }

  SECTION("kinetic factors: Bergeron is akci/akel, less fractionating than at equilibrium") {
    // The replay passes akel as alpha_l and akci as alpha_i (cold, the liquid
    // evaporating at RH_l < 1, the vapor depositing at S_i > 1)
    using WIK = water_isotopes::WaterIsotopeKinetic;
    using WIF = water_isotopes::WaterIsotopeFractionation;
    const water_isotopes::WaterIsotopeParameters<Real> iso_params;
    const Real T = 245, rh_l = 0.9, S_i = 1.1;
    for (const auto sp : {water_isotopes::WaterIsotopologues::HDO, water_isotopes::WaterIsotopologues::H218O}) {
      const Pack1 alpha_l = WIF::alpha_liquid_vapor(Pack1(T), sp, water_isotopes::CondensedOverVapor, iso_params);
      const Pack1 alpha_i = WIF::alpha_ice_vapor(Pack1(T), sp, water_isotopes::CondensedOverVapor, iso_params);
      const Pack1 akel = WIK::alpha_liquid_vapor_kinetic(alpha_l, Pack1(rh_l), sp, iso_params);
      const Pack1 akci = WIK::alpha_ice_vapor_kinetic(alpha_i, Pack1(T), Pack1(S_i), sp, iso_params);
      const Real berg_eq  = wt::part2_beta(WT::Part2BergC2I, alpha_l, alpha_i)[0];
      const Real berg_kin = wt::part2_beta(WT::Part2BergC2I, akel, akci)[0];
      INFO("alpha_l = " << alpha_l[0] << ", alpha_i = " << alpha_i[0] << ", akel = " << akel[0]
           << ", akci = " << akci[0] << ", Bergeron beta: eq " << berg_eq << ", kinetic " << berg_kin);
      REQUIRE(berg_kin == akci[0]/akel[0]);
      REQUIRE(berg_eq > 1);
      REQUIRE(berg_kin < berg_eq);
      REQUIRE(wt::part2_beta(WT::Part2EvapR2V, akel, akci)[0] < wt::part2_beta(WT::Part2EvapR2V, alpha_l, alpha_i)[0]);
      REQUIRE(wt::part2_beta(WT::Part2VapdepV2I, akel, akci)[0] < wt::part2_beta(WT::Part2VapdepV2I, alpha_l, alpha_i)[0]);
      REQUIRE(wt::part2_beta(WT::Part2VapdepV2I, akel, akci)[0] > 1);
    }
  }
}

TEST_CASE("donor_transfers_on_device") {
  using KT = ekat::KokkosTypes<DefaultDevice>;
  KT::view_1d<Real> out("donor_transfers_device", 2);
  Kokkos::parallel_for("donor_transfers_device", 1, KOKKOS_LAMBDA(const int) {
    const Pack1 a[2] = {Pack1(1e-3), Pack1(2e-3)}, b[2] = {Pack1(1.1), Pack1(0.9)};
    Pack1 m[2];
    wt::donor_transfers(Pack1(2e-3), Pack1(1e-2), a, b, m);
    out(0) = m[0][0];
    out(1) = m[1][0];
  });
  const auto out_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);
  const Real amount[2] = {1e-3, 2e-3}, beta[2] = {1.1, 0.9};
  Real moved[2];
  transfer(2e-3, 1e-2, amount, beta, moved);
  REQUIRE(rel_close(out_h(0), moved[0]));
  REQUIRE(rel_close(out_h(1), moved[1]));
}

} // namespace scream
