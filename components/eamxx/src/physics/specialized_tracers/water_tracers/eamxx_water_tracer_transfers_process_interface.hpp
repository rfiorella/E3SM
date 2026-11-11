#ifndef EAMXX_WATER_TRACER_TRANSFERS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_HPP

#include "eamxx_water_tracers_registry.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include <ekat_parameter_list.hpp>
#include <ekat_pack.hpp>

#include <string>
#include <vector>

namespace scream
{
/*
 * Water tracer transfers: moves the water tracers between phases the way
 * SHOC moved bulk water, by replaying the transfers SHOC exports with
 * export_water_transfers=true. Bulk water is only read.
 *
 * Each transfer moves the tracer at the ratio R = tracer/bulk of its donor
 * phase (no fractionation yet). The tracers are those of the water_tracers
 * process, read from the same water_tracers_file.
 *
 * SHOC replay (replay_shoc, default true). The process must run right after
 * SHOC, in the same subcycle: shoc_water_stages holds only SHOC's last call.
 * The tracers are diffused in SHOC's own solve, so only SHOC's condensation is
 * replayed (see shoc_water_transfers.hpp). Afterwards <tracer>_residual_shoc
 * is the bulk water SHOC returned minus the tracer, per phase.
 *
 * P3 replay (replay_p3, default true): not implemented yet, so replay_p3 must
 * be set to false.
 */

class WaterTracerTransfers : public AtmosphereProcess
{
public:
  using KT   = ekat::KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  WaterTracerTransfers (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  std::string name () const override { return "water_tracer_transfers"; }

  void create_requests () override;

protected:

#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif

  void initialize_impl (const RunType /* run_type */) override {}
  void run_impl        (const double dt) override;
  void finalize_impl   () override {}

  // Steps of run_impl
  void replay_shoc ();

protected:

  int m_num_cols;
  int m_num_levs;

  std::shared_ptr<const AbstractGrid> m_grid;

  // Which exported transfers to replay
  bool m_replay_shoc;
  bool m_replay_p3;

  // One entry per tracer in water_tracers_file, sorted by name. Parsed in
  // create_requests, so the parameter can be set after construction.
  std::vector<water_tracers::TracerSpec> m_specs;

}; // class WaterTracerTransfers

} // namespace scream

#endif // EAMXX_WATER_TRACER_TRANSFERS_HPP
