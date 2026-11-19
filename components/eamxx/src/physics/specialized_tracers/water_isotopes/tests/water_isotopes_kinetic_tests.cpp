#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"
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
#include <string>
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

  // The rows of water_isotopes.F90's difrm, per species. Cappa et al. (2003)
  // has a smaller D/D_i for HDO and larger for H218O, H217O differs in the
  // fourth digit and HTO not at all.
  water_isotopes::WaterIsotopeRuntimeOptions cappa;
  cappa.diffusivity = water_isotopes::DiffusivityRatioFormulation::Cappa2003;
  const Params pc(cappa);
  const Real merlivat_row[num_species] = {1.0, 0.9757, 0.9727, 0.9858, 0.9679};
  const Real cappa_row[num_species]    = {1.0, 0.9839, 0.9691, 0.9856, 0.9679};
  for (int isp = 0; isp < num_species; ++isp) {
    INFO("species " << species_names[isp]);
    REQUIRE(p.difrm(static_cast<WaterIsotopologues>(isp))  == merlivat_row[isp]);
    REQUIRE(pc.difrm(static_cast<WaterIsotopologues>(isp)) == cappa_row[isp]);
  }
  for (const auto& st : states) {
    REQUIRE(alpkn(st, WaterIsotopologues::HDO,   pc) > alpkn(st, WaterIsotopologues::HDO,   p));
    REQUIRE(alpkn(st, WaterIsotopologues::H218O, pc) < alpkn(st, WaterIsotopologues::H218O, p));
    REQUIRE(alpkn(st, WaterIsotopologues::H217O, pc) < alpkn(st, WaterIsotopologues::H217O, p));
    REQUIRE(alpkn(st, WaterIsotopologues::HTO,   pc) == alpkn(st, WaterIsotopologues::HTO,   p));
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

// =============================================================================
// Kinetic fractionation of liquid evaporation (akel) and deposition on ice
// (akci), and their helpers heff and ssatf
// =============================================================================
namespace {

// Reference values from the Fortran (wiso_heff, wiso_ssatf, wiso_akel,
// wiso_akci and wiso_kmol of the unmodified water_isotopes.F90, compiled with
// gfortran 13.4 against a stub shr_sys_mod; Cappa2003 rows from a copy with
// its commented-out Cappa et al. (2003) difrm row active instead). Every
// species, H217O and HTO included, with its own D_i/D. The alpha_eq inputs are
// given, not computed.
using DRF = water_isotopes::DiffusivityRatioFormulation;
struct AkelCase { DRF form; WaterIsotopologues sp; Real alpha_eq, rh, akel; };
struct AkciCase { DRF form; WaterIsotopologues sp; Real alpha_eq, T, S, akci; };
struct KmolCase { DRF form; WaterIsotopologues sp; Real rbot, zbot, ustar, alpkn; };
// {rh, heff}
const Real heff_ref[][2] = {
  {0.00000000000000000e+00, 7.50000000000000000e-01},
  {2.99999999999999989e-01, 8.24999999999999956e-01},
  {7.50000000000000000e-01, 9.37500000000000000e-01},
  {9.89999999999999991e-01, 9.97500000000000053e-01},
  {1.00000000000000000e+00, 1.00000000000000000e+00},
  {1.50000000000000000e+00, 1.00000000000000000e+00},
};
// {T [K], ssatf}
const Real ssatf_ref[][2] = {
  {1.50000000000000000e+02, 1.24632000000000009e+00},
  {2.00000000000000000e+02, 1.14632000000000001e+00},
  {2.33150000000000006e+02, 1.08001999999999998e+00},
  {2.53149000000000001e+02, 1.04002200000000000e+00},
  {2.53150000000000006e+02, 1.04001999999999994e+00},
  {2.60000000000000000e+02, 1.02632000000000012e+00},
  {2.73160000000000025e+02, 1.00000000000000000e+00},
  {3.00000000000000000e+02, 1.00000000000000000e+00},
};
// {diffusivity row, species, alpha_eq, rh, akel}
const AkelCase akel_ref[] = {
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 0.00000000000000000e+00, 1.14416768609835762e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 2.99999999999999989e-01, 1.12770223214731624e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 5.99999999999999978e-01, 1.11433873845749076e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 9.00000000000000022e-01, 1.10327608999536042e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 1.00000000000000000e+00, 1.10000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.10000000000000009e+00, 1.30000000000000004e+00, 1.10000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 0.00000000000000000e+00, 1.37260216510693978e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 2.99999999999999989e-01, 1.32533273589730660e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 5.99999999999999978e-01, 1.28835916606935386e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 9.00000000000000022e-01, 1.25864797427364983e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.00000000000000000e+00, 1.25000000000000000e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.30000000000000004e+00, 1.25000000000000000e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 0.00000000000000000e+00, 1.01894824982190002e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 2.99999999999999989e-01, 1.01567605429643937e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 5.99999999999999978e-01, 1.01296523596558319e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 9.00000000000000022e-01, 1.01068274333775432e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 1.00000000000000000e+00, 1.01000000000000001e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01000000000000001e+00, 1.30000000000000004e+00, 1.01000000000000001e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 0.00000000000000000e+00, 1.03940908671653953e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 2.99999999999999989e-01, 1.03412278432191629e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 5.99999999999999978e-01, 1.02975843398652334e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 9.00000000000000022e-01, 1.02609419272570523e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 1.00000000000000000e+00, 1.02499999999999991e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02499999999999991e+00, 1.30000000000000004e+00, 1.02499999999999991e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 0.00000000000000000e+00, 1.01623726788336932e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 2.99999999999999989e-01, 1.01396027711057579e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 5.99999999999999978e-01, 1.01207056476172474e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 9.00000000000000022e-01, 1.01047707031060829e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 1.00000000000000000e+00, 1.01000000000000001e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01000000000000001e+00, 1.30000000000000004e+00, 1.01000000000000001e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 0.00000000000000000e+00, 1.03658828742714837e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 2.99999999999999989e-01, 1.03234417177790516e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 5.99999999999999978e-01, 1.02883386625192141e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 9.00000000000000022e-01, 1.02588220010151532e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 1.00000000000000000e+00, 1.02499999999999991e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02499999999999991e+00, 1.30000000000000004e+00, 1.02499999999999991e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 0.00000000000000000e+00, 1.14623686790526924e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 2.99999999999999989e-01, 1.12898049696126157e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 5.99999999999999978e-01, 1.11499217345706136e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 9.00000000000000022e-01, 1.10342383670424593e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 1.00000000000000000e+00, 1.10000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.10000000000000009e+00, 1.30000000000000004e+00, 1.10000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 0.00000000000000000e+00, 1.37558113071560206e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 2.99999999999999989e-01, 1.32709864394393051e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 5.99999999999999978e-01, 1.28923270368404252e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 9.00000000000000022e-01, 1.25884026846062302e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.00000000000000000e+00, 1.25000000000000000e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.30000000000000004e+00, 1.25000000000000000e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 0.00000000000000000e+00, 1.14202827057157563e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 2.99999999999999989e-01, 1.12637876593005970e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 5.99999999999999978e-01, 1.11366144361557162e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 9.00000000000000022e-01, 1.10312280722954559e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 1.00000000000000000e+00, 1.10000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.10000000000000009e+00, 1.30000000000000004e+00, 1.10000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 0.00000000000000000e+00, 1.36952434631356001e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 2.99999999999999989e-01, 1.32350512172678414e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 5.99999999999999978e-01, 1.28745390051720698e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 9.00000000000000022e-01, 1.25844848250503816e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.00000000000000000e+00, 1.25000000000000000e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.30000000000000004e+00, 1.25000000000000000e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 0.00000000000000000e+00, 1.01970595739216940e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 2.99999999999999989e-01, 1.01615500969534023e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 5.99999999999999978e-01, 1.01321472376970911e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 9.00000000000000022e-01, 1.01074004741371692e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 1.00000000000000000e+00, 1.01000000000000001e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01000000000000001e+00, 1.30000000000000004e+00, 1.01000000000000001e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 0.00000000000000000e+00, 1.04019754164880229e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 2.99999999999999989e-01, 1.03461929955721010e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 5.99999999999999978e-01, 1.03001626355815601e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 9.00000000000000022e-01, 1.02615325778546373e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 1.00000000000000000e+00, 1.02499999999999991e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02499999999999991e+00, 1.30000000000000004e+00, 1.02499999999999991e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 0.00000000000000000e+00, 1.01627812123722516e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 2.99999999999999989e-01, 1.01398615794920599e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 5.99999999999999978e-01, 1.01208407074323814e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 9.00000000000000022e-01, 1.01048017723676353e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 1.00000000000000000e+00, 1.01000000000000001e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01000000000000001e+00, 1.30000000000000004e+00, 1.01000000000000001e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 0.00000000000000000e+00, 1.03663079344547260e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 2.99999999999999989e-01, 1.03237099961642254e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 5.99999999999999978e-01, 1.02884782335102076e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 9.00000000000000022e-01, 1.02588540248262006e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 1.00000000000000000e+00, 1.02499999999999991e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02499999999999991e+00, 1.30000000000000004e+00, 1.02499999999999991e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 0.00000000000000000e+00, 1.14623686790526924e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 2.99999999999999989e-01, 1.12898049696126157e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 5.99999999999999978e-01, 1.11499217345706136e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 9.00000000000000022e-01, 1.10342383670424593e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 1.00000000000000000e+00, 1.10000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.10000000000000009e+00, 1.30000000000000004e+00, 1.10000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 0.00000000000000000e+00, 1.37558113071560206e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 2.99999999999999989e-01, 1.32709864394393051e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 5.99999999999999978e-01, 1.28923270368404252e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 9.00000000000000022e-01, 1.25884026846062302e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.00000000000000000e+00, 1.25000000000000000e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.30000000000000004e+00, 1.25000000000000000e+00},
};
// {diffusivity row, species, alpha_eq, T [K], S (< 0: wiso_ssatf(T)), akci}
const AkciCase akci_ref[] = {
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 1.50000000000000000e+02, -1, 1.14784840271725286e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.00000000000000000e+02, -1, 1.16579168024570845e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.33150000000000006e+02, -1, 1.17990323702744382e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.53149000000000001e+02, -1, 1.18947734324327303e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.53150000000000006e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.60000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.73160000000000025e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 3.00000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 1.50000000000000000e+02, -1, 1.25485882656898773e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.00000000000000000e+02, -1, 1.28698072708552114e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.33150000000000006e+02, -1, 1.31268939862160838e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.53149000000000001e+02, -1, 1.33036055800857289e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.53150000000000006e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.60000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.73160000000000025e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 3.00000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 1.50000000000000000e+02, -1, 1.00635070088298573e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.00000000000000000e+02, -1, 1.00939699151313378e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.33150000000000006e+02, -1, 1.01174015322600774e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.53149000000000001e+02, -1, 1.01330426656146022e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.53150000000000006e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.60000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.73160000000000025e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 3.00000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 1.50000000000000000e+02, -1, 1.01814624406071186e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.00000000000000000e+02, -1, 1.02231299066698345e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.33150000000000006e+02, -1, 1.02552401841078633e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.53149000000000001e+02, -1, 1.02767037614289025e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.53150000000000006e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.60000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.73160000000000025e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 3.00000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 1.50000000000000000e+02, -1, 1.00909261633701086e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.00000000000000000e+02, -1, 1.01117687068880691e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.33150000000000006e+02, -1, 1.01277733117167501e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.53149000000000001e+02, -1, 1.01384436259682698e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.53150000000000006e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.60000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.73160000000000025e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 3.00000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 1.50000000000000000e+02, -1, 1.02095290240591918e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.00000000000000000e+02, -1, 1.02413875226394890e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.33150000000000006e+02, -1, 1.02658966460673007e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.53149000000000001e+02, -1, 1.02822589934675834e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.53150000000000006e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.60000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.73160000000000025e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 3.00000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 1.50000000000000000e+02, -1, 1.14570168961379815e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.00000000000000000e+02, -1, 1.16436063259302292e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.33150000000000006e+02, -1, 1.17905191498640871e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.53149000000000001e+02, -1, 1.18902781988931494e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.53150000000000006e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.60000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.73160000000000025e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 3.00000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 1.50000000000000000e+02, -1, 1.25229363995964427e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.00000000000000000e+02, -1, 1.28523691026160281e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.33150000000000006e+02, -1, 1.31163576460113007e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.53149000000000001e+02, -1, 1.32979826951893099e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.53150000000000006e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.60000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.73160000000000025e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 3.00000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 8.00000000000000044e-01, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.00000000000000000e+00, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.05000000000000004e+00, 1.18700588973641019e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.30000000000000004e+00, 1.13954630530807832e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 2.00000000000000000e+00, 1.07628809118904267e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 2.50000000000000000e+00, 1.07628809118904267e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.20000000000000018e+00, 260.0, 1.30000000000000004e+00, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 8.00000000000000044e-01, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.00000000000000000e+00, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.05000000000000004e+00, 1.32578101035276097e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.30000000000000004e+00, 1.24020677916998623e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 2.00000000000000000e+00, 1.13272993077353057e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 2.50000000000000000e+00, 1.13272993077353057e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.35000000000000009e+00, 260.0, 1.30000000000000004e+00, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 8.00000000000000044e-01, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.00000000000000000e+00, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.05000000000000004e+00, 1.01290246226000291e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.30000000000000004e+00, 1.00491516598809905e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 2.00000000000000000e+00, 9.93399909443074991e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 2.50000000000000000e+00, 9.93399909443074991e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.01499999999999990e+00, 260.0, 1.30000000000000004e+00, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 8.00000000000000044e-01, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.00000000000000000e+00, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.05000000000000004e+00, 1.02711877562975973e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.30000000000000004e+00, 1.01618576890067858e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 2.00000000000000000e+00, 1.00053028411644251e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 2.50000000000000000e+00, 1.00053028411644251e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.02999999999999980e+00, 260.0, 1.30000000000000004e+00, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 8.00000000000000044e-01, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.00000000000000000e+00, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.05000000000000004e+00, 1.01357035357676906e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.30000000000000004e+00, 1.00810904440827698e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 2.00000000000000000e+00, 1.00018692522990782e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 2.50000000000000000e+00, 1.00018692522990782e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.01499999999999990e+00, 260.0, 1.30000000000000004e+00, 1.01499999999999990e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 8.00000000000000044e-01, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.00000000000000000e+00, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.05000000000000004e+00, 1.02780555287766306e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.30000000000000004e+00, 1.01945180722891560e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 2.00000000000000000e+00, 1.00741541819624958e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 2.50000000000000000e+00, 1.00741541819624958e+00},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.02999999999999980e+00, 260.0, 1.30000000000000004e+00, 1.02999999999999980e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 8.00000000000000044e-01, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.00000000000000000e+00, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.05000000000000004e+00, 1.18645198949314157e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.30000000000000004e+00, 1.13707658709240178e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 2.00000000000000000e+00, 1.07152543936528466e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 2.50000000000000000e+00, 1.07152543936528466e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.20000000000000018e+00, 260.0, 1.30000000000000004e+00, 1.20000000000000018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 8.00000000000000044e-01, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.00000000000000000e+00, 1.35000000000000009e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.05000000000000004e+00, 1.32509006181186018e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.30000000000000004e+00, 1.23728203073785425e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 2.00000000000000000e+00, 1.12745588679408071e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 2.50000000000000000e+00, 1.12745588679408071e+00},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.35000000000000009e+00, 260.0, 1.30000000000000004e+00, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 1.50000000000000000e+02, -1, 1.15007698228348310e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.00000000000000000e+02, -1, 1.16727535489924850e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.33150000000000006e+02, -1, 1.18078495671917660e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.53149000000000001e+02, -1, 1.18994259179724460e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.53150000000000006e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.60000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 2.73160000000000025e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 3.00000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 1.50000000000000000e+02, -1, 1.25752278492166703e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.00000000000000000e+02, -1, 1.28878914300611491e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.33150000000000006e+02, -1, 1.31378083452465222e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.53149000000000001e+02, -1, 1.33094256945340339e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.53150000000000006e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.60000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 2.73160000000000025e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 3.00000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 1.50000000000000000e+02, -1, 1.00558687555921700e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.00000000000000000e+02, -1, 1.00890055561182623e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.33150000000000006e+02, -1, 1.01145059533744130e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.53149000000000001e+02, -1, 1.01315338842016089e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.53150000000000006e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.60000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 2.73160000000000025e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 3.00000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 1.50000000000000000e+02, -1, 1.01736441499968566e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.00000000000000000e+02, -1, 1.02180377214126383e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.33150000000000006e+02, -1, 1.02522651811100984e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.53149000000000001e+02, -1, 1.02751518985627799e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.53150000000000006e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.60000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 2.73160000000000025e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 3.00000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 1.50000000000000000e+02, -1, 1.00905119207053851e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.00000000000000000e+02, -1, 1.01115000604877037e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.33150000000000006e+02, -1, 1.01276168786814402e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.53149000000000001e+02, -1, 1.01383622050788924e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.53150000000000006e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.60000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 2.73160000000000025e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 3.00000000000000000e+02, -1, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 1.50000000000000000e+02, -1, 1.02091049868410577e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.00000000000000000e+02, -1, 1.02411119448431998e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.33150000000000006e+02, -1, 1.02657359170787532e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.53149000000000001e+02, -1, 1.02821752462689098e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.53150000000000006e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.60000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 2.73160000000000025e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 3.00000000000000000e+02, -1, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 1.50000000000000000e+02, -1, 1.14570168961379815e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.00000000000000000e+02, -1, 1.16436063259302292e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.33150000000000006e+02, -1, 1.17905191498640871e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.53149000000000001e+02, -1, 1.18902781988931494e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.53150000000000006e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.60000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 2.73160000000000025e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 3.00000000000000000e+02, -1, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 1.50000000000000000e+02, -1, 1.25229363995964427e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.00000000000000000e+02, -1, 1.28523691026160281e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.33150000000000006e+02, -1, 1.31163576460113007e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.53149000000000001e+02, -1, 1.32979826951893099e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.53150000000000006e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.60000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 2.73160000000000025e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 3.00000000000000000e+02, -1, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 8.00000000000000044e-01, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.00000000000000000e+00, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.05000000000000004e+00, 1.18757927004502362e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 1.30000000000000004e+00, 1.14211176426817484e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 2.00000000000000000e+00, 1.08125829937268203e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 240.0, 2.50000000000000000e+00, 1.08125829937268203e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.20000000000000018e+00, 260.0, 1.30000000000000004e+00, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 8.00000000000000044e-01, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.00000000000000000e+00, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.05000000000000004e+00, 1.32649633821571244e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 1.30000000000000004e+00, 1.24324609403124775e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 2.00000000000000000e+00, 1.13823642829598537e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 240.0, 2.50000000000000000e+00, 1.13823642829598537e+00},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.35000000000000009e+00, 260.0, 1.30000000000000004e+00, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 8.00000000000000044e-01, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.00000000000000000e+00, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.05000000000000004e+00, 1.01271591410501527e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 1.30000000000000004e+00, 1.00402595006281392e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 2.00000000000000000e+00, 9.91519076659442389e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 240.0, 2.50000000000000000e+00, 9.91519076659442389e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.01499999999999990e+00, 260.0, 1.30000000000000004e+00, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 8.00000000000000044e-01, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.00000000000000000e+00, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.05000000000000004e+00, 1.02692695473251017e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 1.30000000000000004e+00, 1.01527650418590087e+00},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 2.00000000000000000e+00, 9.98622380071031879e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 240.0, 2.50000000000000000e+00, 9.98622380071031879e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.02999999999999980e+00, 260.0, 1.30000000000000004e+00, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 8.00000000000000044e-01, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.00000000000000000e+00, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.05000000000000004e+00, 1.01356028368794315e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 1.30000000000000004e+00, 1.00806077048290832e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 2.00000000000000000e+00, 1.00008397480755784e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 240.0, 2.50000000000000000e+00, 1.00008397480755784e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.01499999999999990e+00, 260.0, 1.30000000000000004e+00, 1.01499999999999990e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 8.00000000000000044e-01, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.00000000000000000e+00, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.05000000000000004e+00, 1.02779519814868370e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 1.30000000000000004e+00, 1.01940244090838839e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 2.00000000000000000e+00, 1.00731097439968220e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 240.0, 2.50000000000000000e+00, 1.00731097439968220e+00},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.02999999999999980e+00, 260.0, 1.30000000000000004e+00, 1.02999999999999980e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 8.00000000000000044e-01, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.00000000000000000e+00, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.05000000000000004e+00, 1.18645198949314157e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 1.30000000000000004e+00, 1.13707658709240178e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 2.00000000000000000e+00, 1.07152543936528466e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 240.0, 2.50000000000000000e+00, 1.07152543936528466e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.20000000000000018e+00, 260.0, 1.30000000000000004e+00, 1.20000000000000018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 8.00000000000000044e-01, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.00000000000000000e+00, 1.35000000000000009e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.05000000000000004e+00, 1.32509006181186018e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 1.30000000000000004e+00, 1.23728203073785425e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 2.00000000000000000e+00, 1.12745588679408071e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 240.0, 2.50000000000000000e+00, 1.12745588679408071e+00},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.35000000000000009e+00, 260.0, 1.30000000000000004e+00, 1.35000000000000009e+00},
};
// {diffusivity row, species, rbot [kg/m3], zbot [m], ustar [m/s], alpkn} (wiso_kmol)
const KmolCase kmol_ref[] = {
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.94183761011705247e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.94400080965430977e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.95143999809103841e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.97153619682395354e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.96123966424890761e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.93897572033770493e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.93453751637355831e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.93697050268776083e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.94533849544473991e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.96795986269803258e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.95637528231441360e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H218O, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.93133649646904981e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.96621958059428881e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.96747892219608800e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.97180801326420396e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.98347483853695539e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.97748725515724133e-01},
  {DRF::Merlivat1978, WaterIsotopologues::H217O, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.96452269196754559e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.92280149873662309e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.92566740363336697e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.93552629457721048e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.96220856580812164e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.94855525203735103e-01},
  {DRF::Merlivat1978, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.91906457786925988e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.96165566043958006e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.96308451578339249e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.96799669652107090e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.98124082907810850e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.97444586603472683e-01},
  {DRF::Cappa2003, WaterIsotopologues::HDO, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.95973690583753180e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.92574197890904841e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.92849950847241813e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.93798511820377928e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.96364977167570331e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.95051455131067741e-01},
  {DRF::Cappa2003, WaterIsotopologues::H218O, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.92213825400145133e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.96573966186037663e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.96701683572829622e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.97140726166415581e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.98323993685514166e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.97716743589635913e-01},
  {DRF::Cappa2003, WaterIsotopologues::H217O, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.96401936044860048e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.19999999999999996e+00, 1.00000000000000000e+01, 5.00000000000000028e-02, 9.92280149873662309e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.29000000000000004e+00, 5.00000000000000000e+00, 1.00000000000000006e-01, 9.92566740363336697e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.14999999999999991e+00, 3.00000000000000000e+01, 2.00000000000000011e-01, 9.93552629457721048e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.19999999999999996e+00, 1.00000000000000000e+01, 2.99999999999999989e-01, 9.96220856580812164e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.16999999999999993e+00, 2.00000000000000000e+01, 5.00000000000000000e-01, 9.94855525203735103e-01},
  {DRF::Cappa2003, WaterIsotopologues::HTO, 1.25000000000000000e+00, 1.00000000000000000e+01, 1.00000000000000000e+00, 9.91906457786925988e-01},
};


Real akel1(const Real alpha_eq, const Real rh, const WaterIsotopologues sp, const Params& p)
{
  return WIK::alpha_liquid_vapor_kinetic(Pack1(alpha_eq), Pack1(rh), sp, p)[0];
}
Real akci1(const Real alpha_eq, const Real T, const Real S, const WaterIsotopologues sp, const Params& p)
{
  return WIK::alpha_ice_vapor_kinetic(Pack1(alpha_eq), Pack1(T), Pack1(S), sp, p)[0];
}

// The species' equilibrium alphas (the default formulations)
Real alpha_l(const Real T, const WaterIsotopologues sp, const Params& p)
{
  return water_isotopes::WaterIsotopeFractionation::alpha_liquid_vapor(
      Pack1(T), sp, water_isotopes::CondensedOverVapor, p)[0];
}
Real alpha_i(const Real T, const WaterIsotopologues sp, const Params& p)
{
  return water_isotopes::WaterIsotopeFractionation::alpha_ice_vapor(
      Pack1(T), sp, water_isotopes::CondensedOverVapor, p)[0];
}

// The fractionating species
const WaterIsotopologues fractionating_species[] = {
  WaterIsotopologues::HDO, WaterIsotopologues::H218O, WaterIsotopologues::H217O, WaterIsotopologues::HTO};

Params params_of(const DRF form)
{
  water_isotopes::WaterIsotopeRuntimeOptions o;
  o.diffusivity = form;
  return Params(o);
}

} // namespace

// -----------------------------------------------------------------------------
TEST_CASE("akel-akci-fortran-reference") {
  // One ulp of pow/division differences between compilers at most
  const Real rtol = std::is_same<Real, double>::value ? Real(1e-14) : Real(4)*std::numeric_limits<Real>::epsilon();

  SECTION("heff and ssatf") {
    for (const auto& c : heff_ref) {
      INFO("rh = " << c[0]);
      require_within(WIK::effective_humidity(Pack1(c[0]))[0], c[1], rtol, 0);
    }
    for (const auto& c : ssatf_ref) {
      INFO("T = " << c[0]);
      require_within(WIK::ice_supersaturation_empirical(Pack1(c[0]))[0], c[1], rtol, 0);
    }
  }
  SECTION("akel, every species, both diffusivity rows") {
    for (const auto& c : akel_ref) {
      INFO("species " << species_names[etoi(c.sp)] << ", row " << etoi(c.form)
           << ", alpha_eq = " << c.alpha_eq << ", rh = " << c.rh);
      require_within(akel1(c.alpha_eq, c.rh, c.sp, params_of(c.form)), c.akel, rtol, 0);
    }
  }
  SECTION("akci, with the empirical supersaturation and with a given RH over ice") {
    for (const auto& c : akci_ref) {
      const Real S = c.S < 0 ? WIK::ice_supersaturation_empirical(Pack1(c.T))[0] : c.S;
      INFO("species " << species_names[etoi(c.sp)] << ", row " << etoi(c.form)
           << ", alpha_eq = " << c.alpha_eq << ", T = " << c.T << ", S = " << S);
      require_within(akci1(c.alpha_eq, c.T, S, c.sp, params_of(c.form)), c.akci, rtol, 0);
    }
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("alpha-kinetic-fortran-reference") {
  // wiso_kmol's alpkn, every species with its own D_i/D, both diffusivity
  // rows, smooth and rough regimes (log and pow of different libms: a little
  // looser than akel/akci)
  const Real rtol = std::is_same<Real, double>::value ? Real(1e-13) : Real(4)*std::numeric_limits<Real>::epsilon();
  for (const auto& c : kmol_ref) {
    INFO(std::setprecision(17) << "species " << species_names[etoi(c.sp)] << ", row " << etoi(c.form)
         << ", rbot = " << c.rbot << ", zbot = " << c.zbot << ", ustar = " << c.ustar);
    const KineticState st = {c.rbot, c.zbot, c.ustar};
    require_within(alpkn(st, c.sp, params_of(c.form)), c.alpkn, rtol, 0);
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("akel-properties") {
  for (const auto& p : {Params(), Params([] {
         water_isotopes::WaterIsotopeRuntimeOptions o;
         o.diffusivity = water_isotopes::DiffusivityRatioFormulation::Cappa2003;
         return o; }())}) {
    // heff: in [0.75, 1]; rh < 0 counts as 0, rh > 1 as 1
    REQUIRE(WIK::effective_humidity(Pack1(-0.3))[0] == Real(0.75));
    REQUIRE(WIK::effective_humidity(Pack1(0.0))[0]  == Real(0.75));
    REQUIRE(WIK::effective_humidity(Pack1(1.0))[0]  == Real(1));
    REQUIRE(WIK::effective_humidity(Pack1(3.0))[0]  == Real(1));

    for (const auto sp : fractionating_species) {
      for (const Real T : {255.0, 270.0, 300.0}) {
        const Real a = alpha_l(T, sp, p);
        INFO("species " << species_names[etoi(sp)] << ", T = " << T << ", alpha_l = " << a);
        // = alpha_eq at and above saturation, exactly
        for (const Real rh : {1.0, 1.01, 1.5}) {
          REQUIRE(akel1(a, rh, sp, p) == a);
        }
        // > alpha_eq below it, growing as rh falls
        Real prev = a;
        for (const Real rh : {0.99, 0.9, 0.7, 0.5, 0.2, 0.0}) {
          const Real ak = akel1(a, rh, sp, p);
          INFO("rh = " << rh << ", akel = " << ak);
          REQUIRE(ak > prev);
          prev = ak;
        }
        // heff >= 0.75 bounds the kinetic part (HTO's, the largest, is ~15% here)
        REQUIRE(prev < a*Real(1.2));
      }
    }
    // H216O: alpha_eq (1), exactly
    REQUIRE(akel1(1.0, 0.3, WaterIsotopologues::H216O, p) == Real(1));
  }

  SECTION("diffusivity_formulation applies, per species") {
    const Params pm = params_of(DRF::Merlivat1978), pc = params_of(DRF::Cappa2003);
    // Cappa et al. (2003): D/D_i smaller for HDO, larger for H218O; H217O's
    // differs in the fourth digit (0.9856 vs 0.9858), HTO's not at all
    for (const auto sp : fractionating_species) {
      const Real a = alpha_l(280.0, sp, pm);
      const Real km = akel1(a, 0.5, sp, pm), kc = akel1(a, 0.5, sp, pc);
      INFO("species " << species_names[etoi(sp)] << ", akel Merlivat " << km << ", Cappa " << kc);
      if (sp == WaterIsotopologues::HDO) {
        REQUIRE(kc < km);
      } else if (sp == WaterIsotopologues::HTO) {
        REQUIRE(kc == km);
      } else {
        REQUIRE(kc > km);
      }
    }
  }
}

// -----------------------------------------------------------------------------
TEST_CASE("akci-properties") {
  const Params p;
  const Real tkini = WIK::tkini;
  REQUIRE(tkini == Real(253.15));

  SECTION("ssatf: 1 at and above the triple point, capped at ssatmx") {
    REQUIRE(WIK::ice_supersaturation_empirical(Pack1(273.16))[0] == Real(1));
    REQUIRE(WIK::ice_supersaturation_empirical(Pack1(300.0))[0]  == Real(1));
    REQUIRE(WIK::ice_supersaturation_empirical(Pack1(-300.0))[0] == Real(2));  // (not a temperature)
    Real prev = 1;
    for (const Real T : {270.0, 253.15, 230.0, 200.0, 180.0}) {
      const Real S = WIK::ice_supersaturation_empirical(Pack1(T))[0];
      REQUIRE(S > prev);
      REQUIRE(S < Real(1.2));
      prev = S;
    }
  }

  for (const auto sp : fractionating_species) {
    INFO("species " << species_names[etoi(sp)]);
    {  // hard step at tkini
      for (const Real T : {tkini, Real(253.2), Real(260), Real(273.15)}) {
        const Real a = alpha_i(T, sp, p);
        for (const Real S : {1.0, 1.04, 1.3, 2.5}) {
          INFO("T = " << T << ", S = " << S);
          REQUIRE(akci1(a, T, S, sp, p) == a);
        }
      }
      const Real T = tkini - Real(1e-6);
      const Real a = alpha_i(T, sp, p);
      const Real S = WIK::ice_supersaturation_empirical(Pack1(T))[0];
      REQUIRE(akci1(a, T, S, sp, p) < a);
    }
    {  // below tkini: less fractionating as S grows
      for (const Real T : {240.0, 245.0, 250.0, 253.0}) {
        const Real a = alpha_i(T, sp, p);
        INFO("T = " << T << ", alpha_i = " << a);
        // = alpha_eq at S <= 1, exactly
        REQUIRE(akci1(a, T, 1.0, sp, p) == a);
        REQUIRE(akci1(a, T, 0.8, sp, p) == a);
        // < alpha_eq above, and falling as S grows (S is capped at ssatmx)
        Real prev = a;
        for (const Real S : {1.01, 1.1, 1.3, 1.6, 2.0}) {
          const Real ak = akci1(a, T, S, sp, p);
          INFO("S = " << S << ", akci = " << ak);
          REQUIRE(ak < prev);
          prev = ak;
        }
        REQUIRE(akci1(a, T, 3.0, sp, p) == akci1(a, T, 2.0, sp, p));
        // With the empirical S (<= 1.07 here) the ice is still enriched:
        // 1 < akci < alpha_eq
        const Real ak = akci1(a, T, WIK::ice_supersaturation_empirical(Pack1(T))[0], sp, p);
        REQUIRE(ak > 1);
        REQUIRE(ak < a);
      }
    }
  }
  // H216O: alpha_eq (1), exactly
  REQUIRE(akci1(1.0, 230.0, 1.3, WaterIsotopologues::H216O, p) == Real(1));
}

// -----------------------------------------------------------------------------
TEST_CASE("akel-akci-fpe-safety") {
  FpeTrap trap;
  const Params p;
  const Real nan = std::numeric_limits<Real>::quiet_NaN();

  // Lanes 1 and 3 are padding (outside the mask), with inputs that would
  // trap or give garbage if the math touched them
  Pack4 alpha, rh, T, S;
  alpha[0] = 1.08; alpha[1] = 0;   alpha[2] = 1.2;  alpha[3] = nan;
  rh[0]    = 0.6;  rh[1]    = nan; rh[2]    = 0.0;  rh[3]    = -1e30;
  T[0]     = 240;  T[1]     = 0;   T[2]     = 230;  T[3]     = nan;
  S[0]     = 1.2;  S[1]     = nan; S[2]     = 1.05; S[3]     = 0;
  ekat::Mask<4> mask(false);
  mask.set(0, true);
  mask.set(2, true);
  for (const auto sp : fractionating_species) {
    INFO("species " << species_names[etoi(sp)]);
    const Pack4 akel = WIK::alpha_liquid_vapor_kinetic(alpha, rh, sp, p, mask);
    const Pack4 akci = WIK::alpha_ice_vapor_kinetic(alpha, T, S, sp, p, mask);
    for (int l = 0; l < Pack4::n; ++l) {
      INFO("lane " << l);
      if (mask[l]) {
        REQUIRE(akel[l] == akel1(alpha[l], rh[l], sp, p));
        REQUIRE(akci[l] == akci1(alpha[l], T[l], S[l], sp, p));
      } else if (std::isnan(alpha[l])) {
        REQUIRE(std::isnan(akel[l]));
        REQUIRE(std::isnan(akci[l]));
      } else {
        REQUIRE(akel[l] == alpha[l]);
        REQUIRE(akci[l] == alpha[l]);
      }
    }
  }

  SECTION("an extrapolated fit far beyond liquid alphas: akel stays finite and positive") {
    // alpha_eq*(D/D_i)^0.58*(1 - heff) > 1 makes akel's denominator negative;
    // the floor keeps it finite
    for (const auto sp : fractionating_species) {
      const Real a = 6;
      const Real ak = akel1(a, 0.0, sp, p);
      INFO("species " << species_names[etoi(sp)] << ", alpha_eq = " << a << ", akel = " << ak);
      REQUIRE(std::isfinite(ak));
      REQUIRE(ak > a);
    }
  }
}

} // namespace scream
