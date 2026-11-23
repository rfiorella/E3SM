#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_water_transfers.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_team_policy_utils.hpp>

#include <cmath>
#include <limits>

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * p3::water_transfer_phases (p3_water_transfers.hpp) against P3's own updates.
 * For each transfer, P3 must take the mass from the table's source phase, give
 * it to the table's destination, and heat by
 *   dtheta/inv_exner = (H(to) - H(from))*amount/cp
 * (H(Vapor) = 0, H(CloudLiquid) = H(Rain) = LatVap, H(Ice) = LatVap + LatIce).
 *  - Part2* rates: each P3 rate argument is set alone in update_prognostic_ice
 *    and update_prognostic_liquid (P3 calls them in that order), and mapped to
 *    the component P3 exports it as (the export_transfers block after them in
 *    p3_main_impl_part2.hpp), with and without use_hetfrz_classnuc.
 *  - HomFreeze*: homogeneous_freezing on one level holding only cloud, and one
 *    holding only rain, with the amount from P3's own export.
 * The clips are inline in p3_main_part1/2/3 and cannot be called alone; the
 * water_tracer_transfers tests check all components end to end (P3's entry
 * temperature, rebuilt with this table, against T_mid before P3).
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestWaterTransferPhases : public UnitWrap::UnitTest<D>::Base
{
  using WT = WaterTransfer;
  using WE = WaterEntry;

  // Latent enthalpy of a WaterEntry phase relative to vapor [J/kg]
  static Scalar enthalpy (const int phase)
  {
    constexpr Scalar latvap = C::LatVap.value;
    constexpr Scalar latice = C::LatIce.value;
    return phase==WE::Vapor ? 0 : (phase==WE::Ice ? latvap + latice : latvap);
  }

  // One transfer as P3 applied it: the change of each species (WaterEntry
  // order) [kg/kg], the heating dtheta/inv_exner [K], and the mass moved [kg/kg]
  struct Moved {
    Scalar dq[WE::Count];
    Scalar dT;
    Scalar amount;
  };

  // P3 must have moved m.amount from the table's source of `cmp` to its
  // destination, changed nothing else, and heated as the table says. The
  // tolerances are roundoff of the inputs below; a wrong phase is off by at
  // least LatIce/LatVap ~ 0.13 of the heating scale, or by the amount itself.
  static void check (const std::string& label, const int cmp, const Moved& m,
                     const WaterTransferPhases& phases)
  {
    constexpr Scalar inv_cp = C::INV_CP.value;
    constexpr Scalar latvap = C::LatVap.value;
    const Scalar rtol = 1000*std::numeric_limits<Scalar>::epsilon();
    const int from = phases.from[cmp];
    const int to   = phases.to[cmp];

    INFO(label << " (WaterTransfer component " << cmp << ", from " << from << " to " << to << ")");
    REQUIRE(m.amount > 0);
    for (int p=0; p<WE::Count; ++p) {
      const Scalar expected = (p==to ? m.amount : 0) - (p==from ? m.amount : 0);
      INFO("species " << p << ": dq = " << m.dq[p] << ", expected " << expected);
      REQUIRE(std::abs(m.dq[p] - expected) <= rtol*m.amount);
    }
    const Scalar dT = (enthalpy(to) - enthalpy(from))*m.amount*inv_cp;
    INFO("heating = " << m.dT << " K, table " << dT << " K");
    REQUIRE(std::abs(m.dT - dT) <= rtol*latvap*m.amount*inv_cp);
  }

  // Every component has a source and a destination, and they differ
  static void check_table (const WaterTransferPhases& phases)
  {
    for (int cmp=0; cmp<WT::Count; ++cmp) {
      INFO("WaterTransfer component " << cmp);
      REQUIRE(phases.from[cmp] >= 0);
      REQUIRE(phases.from[cmp] <  WE::Count);
      REQUIRE(phases.to[cmp]   >= 0);
      REQUIRE(phases.to[cmp]   <  WE::Count);
      REQUIRE(phases.from[cmp] != phases.to[cmp]);
    }
  }

  // The P3 rate arguments of update_prognostic_ice/liquid that move mass
  // between species. The last three only act with use_hetfrz_classnuc.
  enum Rate : int {
    VapDep, Nucleat, Sublim, Evap, Accret, Autoconv, Shed, Berg, CollectC,
    HeteroFreeze, CollectR, ImmersFreezeR, Melt, QinucCnt, QchetiCnt, QiCnt,
    NumRates
  };

  void run_update_prognostics (const WaterTransferPhases& phases)
  {
    using results_t = typename KokkosTypes<Device>::template view_2d<Scalar>;

    // The component P3 exports each rate as (p3_main_impl_part2.hpp:
    // Part2NucleatV2I is qv2qi_nucleat_tend + qinuc_cnt with classnuc,
    // Part2FreezeC2I is qcheti_cnt + qicnt with classnuc)
    const int component[NumRates] = {
      WT::Part2VapdepV2I, WT::Part2NucleatV2I, WT::Part2SublimI2V, WT::Part2EvapR2V,
      WT::Part2AccretC2R, WT::Part2AutoconvC2R, WT::Part2ShedC2R,  WT::Part2BergC2I,
      WT::Part2CollectC2I, WT::Part2FreezeC2I, WT::Part2CollectR2I, WT::Part2FreezeR2I,
      WT::Part2MeltI2R,   WT::Part2NucleatV2I, WT::Part2FreezeC2I, WT::Part2FreezeC2I };
    const char* names[NumRates] = {
      "qv2qi_vapdep_tend", "qv2qi_nucleat_tend", "qi2qv_sublim_tend", "qr2qv_evap_tend",
      "qc2qr_accret_tend", "qc2qr_autoconv_tend", "qc2qr_ice_shed_tend", "qc2qi_berg_tend",
      "qc2qi_collect_tend", "qc2qi_hetero_freeze_tend", "qr2qi_collect_tend",
      "qr2qi_immers_freeze_tend", "qi2qr_melt_tend", "qinuc_cnt (classnuc)",
      "qcheti_cnt (classnuc)", "qicnt (classnuc)" };

    // Per rate: dq (WaterEntry order), heating [K], amount [kg/kg]
    results_t results("results", NumRates, WE::Count + 2);

    // A state with every species well above qsmall, so P3 applies every rate
    // in full, and a rate moving 1e-3 kg/kg (heating ~ 1-3 K, far above the
    // roundoff of theta)
    constexpr Scalar dt   = 100;    // [s]
    constexpr Scalar rate = 1e-5;   // [kg/kg/s]
    constexpr Scalar inv_exner_val = 1.1;
    constexpr Scalar th0 = 290;     // [K]
    constexpr Scalar qv0 = 1e-2, qc0 = 2e-3, qr0 = 2e-3, qi0 = 2e-3;  // [kg/kg]
    constexpr Scalar nmltratio = C::nmltratio;

    Kokkos::parallel_for(RangePolicy(0, NumRates), KOKKOS_LAMBDA(const int irate) {
      Pack r[NumRates];
      r[irate] = rate;
      const bool classnuc = irate >= QinucCnt;
      const Scalar nmlt = nmltratio;

      const Pack zero(0), inv_exner(inv_exner_val), rho_qm_cloud(400), inv_rho(1);
      Pack ncheti_cnt(0), nicnt(0), ninuc_cnt(0);
      Pack qcheti_cnt(r[QchetiCnt]), qicnt(r[QiCnt]), qinuc_cnt(r[QinucCnt]);
      Pack th(th0), qv(qv0), qc(qc0), qr(qr0), qi(qi0), ni(1e5), qm(1e-3), bm(1e-6),
           nc(1e8), nr(1e5);
      const Mask context(true), wetgrowth(false);

      Functions::update_prognostic_ice(
        r[HeteroFreeze], r[CollectC], r[Shed], zero, zero, zero, r[CollectR], zero,
        r[ImmersFreezeR], zero, zero, r[Melt], zero, r[Sublim], r[VapDep], r[Nucleat],
        zero, zero, zero, r[Berg], inv_exner, true, wetgrowth, dt, nmlt,
        rho_qm_cloud, ncheti_cnt, nicnt, ninuc_cnt, qcheti_cnt, qicnt, qinuc_cnt,
        th, qv, qi, ni, qm, bm, qc, nc, qr, nr, classnuc, context);
      Functions::update_prognostic_liquid(
        r[Accret], zero, r[Autoconv], zero, zero, zero, r[Evap], zero, zero,
        true, false, inv_rho, inv_exner, dt, th, qv, qc, nc, qr, nr, context);

      results(irate, WE::Vapor)       = qv[0] - qv0;
      results(irate, WE::CloudLiquid) = qc[0] - qc0;
      results(irate, WE::Rain)        = qr[0] - qr0;
      results(irate, WE::Ice)         = qi[0] - qi0;
      results(irate, WE::Count)       = (th[0] - th0)/inv_exner_val;
      results(irate, WE::Count + 1)   = rate*dt;
    });

    const auto results_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
    for (int irate=0; irate<NumRates; ++irate) {
      Moved m;
      for (int p=0; p<WE::Count; ++p) {
        m.dq[p] = results_h(irate, p);
      }
      m.dT     = results_h(irate, WE::Count);
      m.amount = results_h(irate, WE::Count + 1);
      check(names[irate], component[irate], m, phases);
    }
  }

  void run_homogeneous_freezing (const WaterTransferPhases& phases)
  {
    using TPF = ekat::TeamPolicyFactory<ExeSpace>;
    using view_1d_pack = view_1d<Pack>;
    using view_2d_pack = view_2d<Pack>;

    // Two levels below T_homogfrz: level 0 holds cloud only, level 1 rain only
    constexpr int nk = 2;
    const int npack = ekat::npack<Pack>(nk);
    constexpr Scalar inv_exner_val = 1.1;
    constexpr Scalar th0 = 260;   // [K]
    constexpr Scalar q0  = 1e-3;  // [kg/kg]
    constexpr Scalar qi0 = 2e-3;  // [kg/kg]
    const Scalar T_cold = C::T_homogfrz.value - 5;

    view_1d_pack T("T", npack), inv_exner("inv_exner", npack), qc("qc", npack), nc("nc", npack),
                 qr("qr", npack), nr("nr", npack), qi("qi", npack), ni("ni", npack),
                 qm("qm", npack), bm("bm", npack), th("th", npack);
    view_2d_pack transfers("transfers", WT::Count, npack);

    const auto set = [&] (const view_1d_pack& v, const Scalar v0, const Scalar v1) {
      const auto h = Kokkos::create_mirror_view(v);
      const auto s = ekat::scalarize(h);
      Kokkos::deep_copy(h, Pack(0));
      s(0) = v0;
      s(1) = v1;
      Kokkos::deep_copy(v, h);
    };
    set(T, T_cold, T_cold);
    set(inv_exner, inv_exner_val, inv_exner_val);
    set(th, th0, th0);
    set(qc, q0, 0);
    set(nc, 1e8, 0);
    set(qr, 0, q0);
    set(nr, 0, 1e5);
    set(qi, qi0, qi0);
    set(ni, 1e5, 1e5);
    set(qm, 1e-3, 1e-3);
    set(bm, 1e-6, 1e-6);

    const auto policy = TPF::get_default_team_policy(1, npack);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      // Level 0 is the top (kdir = -1)
      const Int nlev = nk, ktop = 0, kbot = nk - 1, kdir = -1;
      Functions::homogeneous_freezing(T, inv_exner, team, nlev, ktop, kbot, kdir,
                                      qc, nc, qr, nr, qi, ni, qm, bm, th, transfers);
    });

    const auto get = [&] (const view_1d_pack& v, const int k) {
      const auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
      return ekat::scalarize(h)(k);
    };
    const auto transfers_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), transfers);
    const auto transfer = [&] (const int cmp, const int k) {
      return transfers_h(cmp, k/Pack::n)[k%Pack::n];
    };

    // Level 0 freezes its cloud, level 1 its rain; neither exports the other
    const int cmp_of_level[nk] = {WT::HomFreezeC2I, WT::HomFreezeR2I};
    const int other[nk]        = {WT::HomFreezeR2I, WT::HomFreezeC2I};
    const Scalar qc0[nk] = {q0, 0}, qr0[nk] = {0, q0};
    for (int k=0; k<nk; ++k) {
      REQUIRE(transfer(other[k], k) == 0);
      Moved m;
      m.amount = transfer(cmp_of_level[k], k);
      REQUIRE(m.amount == q0);
      m.dq[WE::Vapor]       = 0;  // homogeneous_freezing does not touch qv
      m.dq[WE::CloudLiquid] = get(qc, k) - qc0[k];
      m.dq[WE::Rain]        = get(qr, k) - qr0[k];
      m.dq[WE::Ice]         = get(qi, k) - qi0;
      m.dT = (get(th, k) - th0)/inv_exner_val;
      check(k==0 ? "homogeneous freezing of cloud" : "homogeneous freezing of rain",
            cmp_of_level[k], m, phases);
    }
  }

  void run ()
  {
    const auto phases = water_transfer_phases();
    check_table(phases);
    run_update_prognostics(phases);
    run_homogeneous_freezing(phases);
  }
};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_water_transfer_phases", "[p3_unit_tests]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestWaterTransferPhases;

  T t;
  t.run();
}

} // namespace
