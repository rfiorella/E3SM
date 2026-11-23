#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_kinetic.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_parameters.hpp"
#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"

#include "share/core/eamxx_types.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_fpe.hpp>
#include <ekat_kokkos_types.hpp>
#include <ekat_pack.hpp>
#include <ekat_pack_kokkos.hpp>
#include <ekat_pack_utils.hpp>
#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <type_traits>

namespace scream {
namespace {

using Pack1 = ekat::Pack<Real, 1>;
using Pack4 = ekat::Pack<Real, 4>;
using WIK   = water_isotopes::WaterIsotopeKinetic;
using water_isotopes::WaterIsotopologues;
using Params = water_isotopes::WaterIsotopeParameters<Real>;

constexpr int num_species = etoi(WaterIsotopologues::Count);
const char* const species_names[num_species] = {"H216O", "HDO", "H218O", "H217O", "HTO"};

// Enables FE_DIVBYZERO|FE_INVALID for its lifetime (not just in the fpe
// build), so the masked-lane tests prove they raise nothing.
struct FpeTrap {
  int saved;
  FpeTrap () : saved(ekat::get_enabled_fpes()) {
    ekat::disable_all_fpes();
    ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);
  }
  ~FpeTrap () {
    ekat::disable_all_fpes();
    ekat::enable_fpes(saved);
  }
};

// Bottom-level (rbot [kg/m3], zbot [m], ustar [m/s]). The first three are in
// the smooth regime (Reynolds number < 1), the others in the rough one.
struct KineticState {
  Real rbot, zbot, ustar;
};
const KineticState states[] = {
  { 1.20, 10.0, 0.05 },
  { 1.29,  5.0, 0.10 },
  { 1.15, 30.0, 0.20 },
  { 1.20, 10.0, 0.30 },
  { 1.17, 20.0, 0.50 },
  { 1.25, 10.0, 1.00 },
};

// Same Reynolds number as wiso_kmol, to confirm the states span both regimes.
Real reynolds(const KineticState& st)
{
  using PC = scream::physics::Constants<Real>;
  const Real z0  = st.ustar*st.ustar/(Real(81.1)*PC::gravit.value);
  const Real vmu = Real(1.7e-5)/st.rbot;
  return st.ustar*z0/vmu;
}

Real alpkn(const KineticState& st, const WaterIsotopologues sp, const Params& p)
{
  return WIK::alpha_kinetic(Pack1(st.rbot), Pack1(st.zbot), Pack1(st.ustar), sp, p)[0];
}

water_isotopes::WaterIsotopeRuntimeOptions kinetic_none()
{
  water_isotopes::WaterIsotopeRuntimeOptions opts;
  opts.kinetic_fractionation = water_isotopes::KineticFractionationFormulation::None;
  return opts;
}

// |computed - expected| <= atol + rtol*|expected|
void require_within(const Real computed, const Real expected, const Real rtol, const Real atol)
{
  const Real diff  = std::abs(computed - expected);
  const Real bound = atol + rtol*std::abs(expected);
  INFO(std::setprecision(17) << "computed = " << computed << ", expected = " << expected
       << ", |diff| = " << diff << ", bound = " << bound);
  REQUIRE(diff <= bound);
}

// rtol 1e-12 is the spec's bound in double precision. Single precision gets
// 100 float ulps, as in water_isotopes_ocean_flux.
Real rtol_for_real()
{
  return std::is_same<Real, double>::value ? Real(1e-12)
                                           : Real(100)*std::numeric_limits<Real>::epsilon();
}

// Ocean-flux inputs (notation of water_tracers::ocean_flux)
struct OceanState {
  Real f, rho, u, r, s, q, q_i, e_ocn;
};

} // namespace

// -----------------------------------------------------------------------------
TEST_CASE("smooth-and-rough-regimes") {
  const Params p;

  bool any_smooth = false, any_rough = false;
  for (const auto& st : states) {
    const Real re = reynolds(st);
    any_smooth = any_smooth || re < 1;
    any_rough  = any_rough  || re >= 1;

    const Real a_hdo  = alpkn(st, WaterIsotopologues::HDO, p);
    const Real a_h218o = alpkn(st, WaterIsotopologues::H218O, p);
    INFO(std::setprecision(17) << "rbot = " << st.rbot << ", zbot = " << st.zbot
         << ", ustar = " << st.ustar << ", Re = " << re
         << ", alpkn(HDO) = " << a_hdo << ", alpkn(H218O) = " << a_h218o);

    // Kinetic fractionation depletes, more so for H218O (D/D_i is larger).
    REQUIRE(a_hdo > 0);   REQUIRE(a_hdo < 1);
    REQUIRE(a_h218o > 0); REQUIRE(a_h218o < 1);
    // (1-alpkn) of H218O is ~1.125x that of HDO (0.9727 vs 0.9757).
    const Real ratio = (1 - a_h218o)/(1 - a_hdo);
    REQUIRE(ratio > Real(1.12));
    REQUIRE(ratio < Real(1.13));
  }
  REQUIRE(any_smooth);
  REQUIRE(any_rough);

  // Cappa et al. (2003) has a smaller D/D_i for HDO and larger for H218O.
  water_isotopes::WaterIsotopeRuntimeOptions cappa;
  cappa.diffusivity = water_isotopes::DiffusivityRatioFormulation::Cappa2003;
  const Params pc(cappa);
  REQUIRE(pc.difrm(water_isotopes::IsoElement::Hydrogen) == Real(0.9839));
  REQUIRE(pc.difrm(water_isotopes::IsoElement::Oxygen)   == Real(0.9691));
  REQUIRE(p.difrm(water_isotopes::IsoElement::Hydrogen)  == Real(0.9757));
  REQUIRE(p.difrm(water_isotopes::IsoElement::Oxygen)    == Real(0.9727));
  for (const auto& st : states) {
    REQUIRE(alpkn(st, WaterIsotopologues::HDO,   pc) > alpkn(st, WaterIsotopologues::HDO,   p));
    REQUIRE(alpkn(st, WaterIsotopologues::H218O, pc) < alpkn(st, WaterIsotopologues::H218O, p));
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("h217o-hto-mass-scaling-unit") {
  const Real rtol = rtol_for_real();
  const Real atol = Real(1e-15);

  for (const auto& p : {Params(), Params([] {
         water_isotopes::WaterIsotopeRuntimeOptions o;
         o.diffusivity = water_isotopes::DiffusivityRatioFormulation::Cappa2003;
         return o; }())}) {
    for (const auto& st : states) {
      INFO(std::setprecision(17) << "rbot = " << st.rbot << ", zbot = " << st.zbot
           << ", ustar = " << st.ustar);
      const Real k_hdo   = 1 - alpkn(st, WaterIsotopologues::HDO,   p);
      const Real k_h218o = 1 - alpkn(st, WaterIsotopologues::H218O, p);
      const Real k_h217o = 1 - alpkn(st, WaterIsotopologues::H217O, p);
      const Real k_hto   = 1 - alpkn(st, WaterIsotopologues::HTO,   p);

      require_within(k_h217o, 1 - std::pow(1 - k_h218o, Real(0.529)), rtol, atol);
      require_within(k_hto,   1 - std::pow(1 - k_hdo,   Real(2.0)),   rtol, atol);

      // H216O does not fractionate, exactly.
      REQUIRE(alpkn(st, WaterIsotopologues::H216O, p) == Real(1));
    }
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("kinetic-formulation-none-is-identity") {
  const Params p_none(kinetic_none());
  const Params p_on;

  // alpkn == 1 exactly for every species and state.
  for (const auto& st : states) {
    for (int isp = 0; isp < num_species; ++isp) {
      INFO("species " << species_names[isp]);
      REQUIRE(alpkn(st, static_cast<WaterIsotopologues>(isp), p_none) == Real(1));
    }
  }

  // ocean_flux(..., alpha, alpkn=1) equals ocean_flux(..., alpha) bit-for-bit,
  // for both evaporation and dew.
  const OceanState ocean[] = {
    //  f     rho   ustar  re      ssq     qv      q_i     e_ocn (sign only matters)
    {  1.0,  1.20,  0.30, 1.5e-3, 1.4e-2, 1.0e-2, 0.98e-2,  0 },
    {  0.5,  1.24,  0.25, 1.2e-3, 7.6e-3, 6.0e-3, 5.9e-3,   0 },
    {  1.0,  1.20,  0.30, 1.5e-3, 1.0e-2, 1.2e-2, 1.17e-2,  0 },  // dew
  };
  for (const auto& oc : ocean) {
    const Pack1 f(oc.f), rho(oc.rho), u(oc.u), r(oc.r), s(oc.s), q(oc.q), q_i(oc.q_i);
    const Pack1 e = water_tracers::bulk_ocean_flux(f, rho, u, r, s, q);
    const auto evap = water_tracers::is_evaporating(f, e);
    for (int isp = 0; isp < num_species; ++isp) {
      const auto sp = static_cast<WaterIsotopologues>(isp);
      INFO("species " << species_names[isp] << ", E_ocn = " << e[0]);
      const Pack1 alpha(Real(1.01));
      const Pack1 a = WIK::alpha_kinetic(rho, Pack1(10.0), u, sp, p_none, evap);
      const Real with = water_tracers::ocean_flux(water_tracers::OceanFlux::BulkFormula, Real(1),
                            f, rho, u, r, s, q, q_i, e, alpha, a)[0];
      const Real without = water_tracers::ocean_flux(water_tracers::OceanFlux::BulkFormula, Real(1),
                               f, rho, u, r, s, q, q_i, e, alpha)[0];
      REQUIRE(with == without);
    }
    // H216O is the identity with kinetic fractionation on, too.
    REQUIRE((WIK::alpha_kinetic(rho, Pack1(10.0), u, WaterIsotopologues::H216O, p_on, evap)
             == Pack1(1)).all());
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("dew-branch-unchanged") {
  // Dew (E_ocn <= 0) and invalid ocean (f <= 0): the flux does not read alpkn.
  const OceanState ocean[] = {
    //  f     rho   ustar  re      ssq     qv      q_i      e_ocn
    {  1.0,  1.20,  0.30, 1.5e-3, 1.0e-2, 1.2e-2, 1.17e-2, -1.0e-5 },  // dew
    {  0.3,  1.25,  0.10, 1.0e-3, 5.0e-3, 5.0e-3, 4.9e-3,   0.0    },  // E_ocn == 0
    {  0.0,  1.20,  0.30, 1.5e-3, 1.4e-2, 1.0e-2, 0.98e-2,  2.0e-5 },  // invalid ocean
    { -1.0,  1.20,  0.0,  0.0,    0.0,    1.0e-2, 0.98e-2,  0.0    },  // invalid ocean
  };
  const Real alpkns[] = {1.0, 0.999, 0.99, 0.95, 0.5, 1e-3};
  for (const auto kind : {water_tracers::OceanFlux::BulkFormula,
                          water_tracers::OceanFlux::BulkProportional}) {
    for (const auto& oc : ocean) {
      const Pack1 f(oc.f), rho(oc.rho), u(oc.u), r(oc.r), s(oc.s), q(oc.q), q_i(oc.q_i),
                  e(oc.e_ocn), alpha(Real(1.01));
      const Real without = water_tracers::ocean_flux(kind, Real(1),
                               f, rho, u, r, s, q, q_i, e, alpha)[0];
      for (const Real a : alpkns) {
        INFO("f = " << oc.f << ", E_ocn = " << oc.e_ocn << ", alpkn = " << a);
        const Real with = water_tracers::ocean_flux(kind, Real(1),
                              f, rho, u, r, s, q, q_i, e, alpha, Pack1(a))[0];
        REQUIRE(with == without);
      }
    }
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("z-mid-matches-exporter-pattern") {
  // compute_bulk_state's zbot = 0.5*calculate_dz(...) at the bottom level must
  // equal the exporter's Sa_z (calculate_dz -> calculate_z_int(z_surf=0) ->
  // calculate_z_mid, bottom value) bit-for-bit, on the same column state.
  using KT     = ekat::KokkosTypes<DefaultDevice>;
  using PF     = PhysicsFunctions<DefaultDevice>;
  using Pack   = ekat::Pack<Real, SCREAM_PACK_SIZE>;
  using TPF    = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using Member = KT::MemberType;

  constexpr int ncols = 4;
  constexpr int nlevs = 72;
  const int nmid = ekat::PackInfo<Pack::n>::num_packs(nlevs);
  const int nint = ekat::PackInfo<Pack::n>::num_packs(nlevs + 1);

  KT::view_2d<Pack> dp("dp", ncols, nmid), p_mid("p_mid", ncols, nmid),
                    T_mid("T_mid", ncols, nmid), qv("qv", ncols, nmid),
                    dz("dz", ncols, nmid), z_mid("z_mid", ncols, nmid),
                    z_int("z_int", ncols, nint);
  auto dp_h = Kokkos::create_mirror_view(dp);
  auto p_h  = Kokkos::create_mirror_view(p_mid);
  auto T_h  = Kokkos::create_mirror_view(T_mid);
  auto qv_h = Kokkos::create_mirror_view(qv);
  // A plausible column: pressure rising to the surface, temperature and
  // humidity decreasing upward, varied per column.
  for (int i = 0; i < ncols; ++i) {
    Real p_top = 200;
    for (int k = 0; k < nlevs; ++k) {
      const Real dpk = Real(800) + Real(900)*k/nlevs + Real(13)*i;
      const Real pk  = p_top + 0.5*dpk;
      p_top += dpk;
      dp_h(i, k/Pack::n)[k%Pack::n] = dpk;
      p_h (i, k/Pack::n)[k%Pack::n] = pk;
      T_h (i, k/Pack::n)[k%Pack::n] = Real(200) + Real(100)*k/nlevs + Real(1.7)*i;
      qv_h(i, k/Pack::n)[k%Pack::n] = Real(1e-6) + Real(1.5e-2)*std::pow(Real(k)/nlevs, 3) + Real(1e-4)*i;
    }
  }
  Kokkos::deep_copy(dp, dp_h);
  Kokkos::deep_copy(p_mid, p_h);
  Kokkos::deep_copy(T_mid, T_h);
  Kokkos::deep_copy(qv, qv_h);

  KT::view_1d<Real> zbot_local("zbot_local", ncols), zbot_exporter("zbot_exporter", ncols);

  // compute_bulk_state's expression (scalar calculate_dz at the bottom level)
  const int kbot = nlevs - 1;
  Kokkos::parallel_for("zbot_local", KT::RangePolicy(0, ncols), KOKKOS_LAMBDA(const int i) {
    const auto s_dp = ekat::scalarize(ekat::subview(dp, i));
    const auto s_p  = ekat::scalarize(ekat::subview(p_mid, i));
    const auto s_T  = ekat::scalarize(ekat::subview(T_mid, i));
    const auto s_qv = ekat::scalarize(ekat::subview(qv, i));
    const Real dzb  = PF::calculate_dz(s_dp(kbot), s_p(kbot), s_T(kbot), s_qv(kbot));
    zbot_local(i) = 0.5*dzb;
  });

  // The exporter's sequence (atmosphere_surface_coupling_exporter.cpp)
  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(ncols, nlevs);
  Kokkos::parallel_for("zbot_exporter", policy, KOKKOS_LAMBDA(const Member& team) {
    const int i = team.league_rank();
    const auto dz_i    = ekat::subview(dz, i);
    const auto z_int_i = ekat::subview(z_int, i);
    const auto z_mid_i = ekat::subview(z_mid, i);
    PF::calculate_dz(team, ekat::subview(dp, i), ekat::subview(p_mid, i),
                     ekat::subview(T_mid, i), ekat::subview(qv, i), dz_i);
    team.team_barrier();
    const Real z_surf = 0.0;
    PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    team.team_barrier();
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      zbot_exporter(i) = ekat::scalarize(z_mid_i)(nlevs - 1);
    });
  });
  Kokkos::fence();

  auto local_h    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), zbot_local);
  auto exporter_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), zbot_exporter);
  for (int i = 0; i < ncols; ++i) {
    INFO(std::setprecision(17) << "column " << i << ": 0.5*dz = " << local_h(i)
         << ", Sa_z = " << exporter_h(i));
    REQUIRE(local_h(i) > 0);
    REQUIRE(local_h(i) == exporter_h(i));
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("zero-wind-safety") {
  FpeTrap trap;
  const Params p;

  SECTION("ustar*zbot sweep over [1e-4, 30] m2/s") {
    using KT = ekat::KokkosTypes<DefaultDevice>;
    constexpr int nsweep = 61;
    const Real rbots[] = {1.0, 1.2, 1.3};
    const Real zbots[] = {5.0, 10.0, 30.0};
    constexpr int nr = std::size(rbots), nz = std::size(zbots);
    constexpr int npts = nr*nz*nsweep;

    // ustar*zbot spaced logarithmically from 1e-4 to 30 m2/s.
    KT::view_2d<Real> in("in", npts, 3);
    auto in_h = Kokkos::create_mirror_view(in);
    int n = 0;
    for (const Real rb : rbots) {
      for (const Real zb : zbots) {
        for (int j = 0; j < nsweep; ++j) {
          const Real lg = Real(-4) + j*(std::log10(Real(30)) + Real(4))/(nsweep - 1);
          in_h(n, 0) = rb;
          in_h(n, 1) = zb;
          in_h(n, 2) = std::pow(Real(10), lg)/zb;
          ++n;
        }
      }
    }
    Kokkos::deep_copy(in, in_h);

    KT::view_2d<Real> out("alpkn", npts, num_species);
    Kokkos::parallel_for("zero_wind_sweep", KT::RangePolicy(0, npts), KOKKOS_LAMBDA(const int i) {
      for (int isp = 0; isp < num_species; ++isp) {
        out(i, isp) = WIK::alpha_kinetic(Pack1(in(i, 0)), Pack1(in(i, 1)), Pack1(in(i, 2)),
                                         static_cast<WaterIsotopologues>(isp), p)[0];
      }
    });
    Kokkos::fence();
    auto out_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);

    Real amin = 1;
    for (int i = 0; i < npts; ++i) {
      for (int isp = 0; isp < num_species; ++isp) {
        const Real a = out_h(i, isp);
        INFO(std::setprecision(17) << "species " << species_names[isp] << ", rbot = " << in_h(i, 0)
             << ", zbot = " << in_h(i, 1) << ", ustar = " << in_h(i, 2) << ", alpkn = " << a);
        REQUIRE(std::isfinite(a));
        REQUIRE(a > 0);
        REQUIRE(a <= 1);
        amin = std::min(amin, a);
      }
    }
    INFO("min alpkn over the sweep = " << amin);
    REQUIRE(amin > Real(0.9));
  }

  SECTION("masked lanes with ustar = 0") {
    // compute_bulk_state sets ustar = 0 (and f = 0) where the ocean state is
    // invalid. Those lanes must not form log(zbot/z0) with z0 = 0.
    const Pack1 zero(0), rbot(1.2), zbot(10.0);
    for (int isp = 0; isp < num_species; ++isp) {
      INFO("species " << species_names[isp]);
      const auto sp = static_cast<WaterIsotopologues>(isp);
      const Pack1 e = water_tracers::bulk_ocean_flux(zero, rbot, zero, zero, zero, zero);
      const auto evap = water_tracers::is_evaporating(zero, e);
      REQUIRE(!evap.any());
      REQUIRE(WIK::alpha_kinetic(rbot, zbot, zero, sp, p, evap)[0] == Real(1));
      // Masked lanes do not see rbot = zbot = 0 either.
      REQUIRE(WIK::alpha_kinetic(zero, zero, zero, sp, p, evap)[0] == Real(1));
    }

    // Mixed pack: valid smooth and rough lanes next to ustar = 0 lanes.
    Pack4 ustar, rb, zb;
    ustar[0] = 0.05; ustar[1] = 0;   ustar[2] = 0.5; ustar[3] = 0;
    rb[0]    = 1.2;  rb[1]    = 0;   rb[2]    = 1.2; rb[3]    = 1.2;
    zb[0]    = 10;   zb[1]    = 10;  zb[2]    = 10;  zb[3]    = 0;
    const auto mask = ustar > 0;
    for (int isp = 0; isp < num_species; ++isp) {
      INFO("species " << species_names[isp]);
      const auto sp = static_cast<WaterIsotopologues>(isp);
      const Pack4 a = WIK::alpha_kinetic(rb, zb, ustar, sp, p, mask);
      for (int l = 0; l < Pack4::n; ++l) {
        INFO("lane " << l);
        REQUIRE(std::isfinite(a[l]));
        if (mask[l]) {
          const Real a1 = WIK::alpha_kinetic(Pack1(rb[l]), Pack1(zb[l]), Pack1(ustar[l]), sp, p)[0];
          REQUIRE(a[l] == a1);
        } else {
          REQUIRE(a[l] == Real(1));
        }
      }
    }
  }
}

} // namespace scream
