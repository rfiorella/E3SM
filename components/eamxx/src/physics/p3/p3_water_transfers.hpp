#ifndef EAMXX_P3_WATER_TRANSFERS_HPP
#define EAMXX_P3_WATER_TRANSFERS_HPP

namespace scream {
namespace p3 {

/*
 * Components of the p3_water_transfers field, written by P3 in builds with
 * water tracers (EAMXX_ENABLE_WATER_TRACERS).
 *
 * Every component is the rate at which one P3 call moved mass between two
 * bulk water species (v=qv, c=qc, r=qr, i=qi), from the first to the second,
 * averaged over the call [kg/kg/s, wet basis]; times the P3 dt it is the mass
 * moved. (Internally P3 accumulates the dry-basis mass; the P3 interface
 * converts it to this rate.) Together with p3_water_entry (the state at P3
 * entry) and
 * p3_sed_flux (sedimentation), they account for every change P3 makes to qv,
 * qc, qr and qi. They are listed in the order P3 applies them; processes
 * within one stage are applied simultaneously, from the same state.
 *
 * Process transfers (Part2*) use the post-limiter rates P3 applies, so they
 * are >= 0. The other stages move a whole reservoir (qx < qsmall is clipped, or
 * all of qc/qr below T_homogfrz freezes), so their amount can be negative
 * where P3 clips a small negative qx.
 *
 * Kept free of other includes so the water tracer code can use it. The enums
 * are unscoped (scoped by a struct instead) so they can index views directly,
 * on host or device.
 */
struct WaterTransfer {
enum : int {
  // p3_main_part1: clipping of small masses
  Part1ClipC2V,
  Part1ClipR2V,
  Part1ClipI2V,   // qi < qsmall, or qi < 1e-8 and subsaturated w.r.t. ice
  Part1ClipI2R,   // tiny qi above freezing

  // p3_main_part2: process rates
  Part2VapdepV2I,    // vapor deposition
  Part2NucleatV2I,   // deposition/condensation-freezing nucleation (incl. classnuc)
  Part2SublimI2V,    // sublimation
  Part2EvapR2V,      // rain evaporation
  Part2AccretC2R,    // accretion of cloud by rain
  Part2AutoconvC2R,  // autoconversion
  Part2ShedC2R,      // shedding of collected cloud water
  Part2BergC2I,      // Bergeron process
  Part2CollectC2I,   // riming of cloud water
  Part2FreezeC2I,    // immersion (or classnuc contact+immersion) freezing of cloud
  Part2CollectR2I,   // collection of rain by ice
  Part2FreezeR2I,    // immersion freezing of rain
  Part2MeltI2R,      // melting

  // p3_main_part2: clipping of small masses after the update
  Part2ClipC2V,
  Part2ClipR2V,
  Part2ClipI2V,

  // homogeneous_freezing (after sedimentation)
  HomFreezeC2I,
  HomFreezeR2I,

  // p3_main_part3: clipping of small masses
  Part3ClipC2V,
  Part3ClipR2V,
  Part3ClipI2V,

  Count
};
};

// Components of the p3_sed_flux field: downward mass flux of each species out
// of the bottom of a level, averaged over the P3 call [kg/m2/s]. The flux out
// of the lowest level is that species' surface precipitation.
struct SedFlux { enum : int { CloudLiquid, Rain, Ice, Count }; };

// Components of the p3_water_entry field: bulk water at P3 entry [kg/kg, wet].
struct WaterEntry { enum : int { Vapor, CloudLiquid, Rain, Ice, Count }; };

/*
 * The phases of each WaterTransfer component: the WaterEntry component it
 * moves mass from, and the one it moves it to (the names encode them: C2V is
 * cloud liquid to vapor).
 *
 * This is also P3's heating table. P3 heats the air by
 *   dT = (H(to) - H(from))*amount/cp   [K], amount [kg/kg],
 * for every transfer, with H(Vapor) = 0, H(CloudLiquid) = H(Rain) = LatVap and
 * H(Ice) = LatVap + LatIce [J/kg] (P3 adds dT*inv_exner to theta): see
 * update_prognostic_ice/liquid (Part2*), the clips in p3_main_part1/2/3 and
 * homogeneous_freezing. tests/p3_water_transfers_tests.cpp checks the table
 * against update_prognostic_ice/liquid and homogeneous_freezing; the clips are
 * inline in the part functions and are only checked end to end, by the
 * water_tracer_transfers tests.
 */
struct WaterTransferPhases {
  int from[WaterTransfer::Count];
  int to[WaterTransfer::Count];
};

// Build the table. Host only: device code needs a copy (e.g. captured by
// value in a lambda), since an array at namespace scope cannot be indexed in
// a CUDA kernel.
inline WaterTransferPhases water_transfer_phases ()
{
  using WT = WaterTransfer;
  constexpr int v = WaterEntry::Vapor;
  constexpr int c = WaterEntry::CloudLiquid;
  constexpr int r = WaterEntry::Rain;
  constexpr int i = WaterEntry::Ice;
  static_assert(WT::Count == 25, "Update water_transfer_phases to match WaterTransfer.");

  WaterTransferPhases p {};
  const auto set = [&p] (const int cmp, const int from, const int to) {
    p.from[cmp] = from;
    p.to[cmp]   = to;
  };
  set(WT::Part1ClipC2V,     c, v);
  set(WT::Part1ClipR2V,     r, v);
  set(WT::Part1ClipI2V,     i, v);
  set(WT::Part1ClipI2R,     i, r);
  set(WT::Part2VapdepV2I,   v, i);
  set(WT::Part2NucleatV2I,  v, i);
  set(WT::Part2SublimI2V,   i, v);
  set(WT::Part2EvapR2V,     r, v);
  set(WT::Part2AccretC2R,   c, r);
  set(WT::Part2AutoconvC2R, c, r);
  set(WT::Part2ShedC2R,     c, r);
  set(WT::Part2BergC2I,     c, i);
  set(WT::Part2CollectC2I,  c, i);
  set(WT::Part2FreezeC2I,   c, i);
  set(WT::Part2CollectR2I,  r, i);
  set(WT::Part2FreezeR2I,   r, i);
  set(WT::Part2MeltI2R,     i, r);
  set(WT::Part2ClipC2V,     c, v);
  set(WT::Part2ClipR2V,     r, v);
  set(WT::Part2ClipI2V,     i, v);
  set(WT::HomFreezeC2I,     c, i);
  set(WT::HomFreezeR2I,     r, i);
  set(WT::Part3ClipC2V,     c, v);
  set(WT::Part3ClipR2V,     r, v);
  set(WT::Part3ClipI2V,     i, v);
  return p;
}

} // namespace p3
} // namespace scream

#endif // EAMXX_P3_WATER_TRANSFERS_HPP
