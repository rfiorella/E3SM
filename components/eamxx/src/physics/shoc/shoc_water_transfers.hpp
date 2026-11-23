#ifndef EAMXX_SHOC_WATER_TRANSFERS_HPP
#define EAMXX_SHOC_WATER_TRANSFERS_HPP

namespace scream {
namespace shoc {

/*
 * Components of the shoc_water_stages field, written by SHOC when the runtime
 * option export_water_transfers is true.
 *
 * SHOC changes bulk water in two stages. First, the implicit solve diffuses
 * total water qw = qv + qc and every member of turbulence_advected_tracers
 * (bulk qc and the water tracers among them) with the same eddy diffusivity.
 * Then the assumed PDF replaces qc with its SGS cloud liquid, shoc_ql, and
 * qv = qw - qc follows. Each component is the bulk state at the end of one of
 * these stages, at the end of the SHOC call [kg/kg, wet basis]:
 *
 *   VaporAfterTransport     = qw - qc, after the solve (qc = diffused member)
 *   LiquidAfterTransport    = qc, after the solve
 *   LiquidAfterCondensation = shoc_ql, the qc SHOC returns
 *
 * Hence, at every level, SHOC's net condensation (< 0 for evaporation) is
 *   LiquidAfterCondensation - LiquidAfterTransport,
 * and the qv SHOC returns (before its qv postcondition repair) is
 *   VaporAfterTransport + LiquidAfterTransport - LiquidAfterCondensation.
 * The water tracers follow the first stage themselves (they are in the same
 * solve), and replay the second one from these components.
 *
 * Only valid with one SHOC iteration per call (nadv=1), which SHOC enforces
 * when exporting.
 *
 * Kept free of other includes so the water tracer code can use it. The enum
 * is unscoped (scoped by a struct instead) so it can index views directly,
 * on host or device.
 */
struct WaterStage { enum : int { VaporAfterTransport, LiquidAfterTransport, LiquidAfterCondensation, Count }; };

} // namespace shoc
} // namespace scream

#endif // EAMXX_SHOC_WATER_TRANSFERS_HPP
