#ifndef EAMXX_P3_WATER_TRANSFERS_HPP
#define EAMXX_P3_WATER_TRANSFERS_HPP

namespace scream {
namespace p3 {

/*
 * Components of the p3_water_transfers field, written by P3 when the runtime
 * option export_water_transfers is true.
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

} // namespace p3
} // namespace scream

#endif // EAMXX_P3_WATER_TRANSFERS_HPP
