/*
 * The AD test driver, plus a test-only atm process, wtt_dummy_modifier, that
 * claims to update bulk water (qc) or a tracer group. The driver link of
 * water_tracer_transfers must refuse it between SHOC/P3 and the replay.
 * No real process that updates bulk water is linked into these tests, and the
 * dummy never runs: the link has to stop during initialization.
 */

#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_units.hpp>

namespace scream {

class WttDummyModifier : public AtmosphereProcess
{
public:
  WttDummyModifier (const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params) {}

  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }
  std::string name () const override { return "wtt_dummy_modifier"; }

  // What it claims to update: parameter 'updates' is a bulk field name (qc),
  // or the name of a group (tracers)
  void create_requests () override
  {
    using namespace ekat::units;
    using namespace ShortFieldTagsNames;

    const auto grid = m_grids_manager->get_grid("physics");
    const auto what = m_params.get<std::string>("updates");
    if (what=="tracers") {
      add_group<Updated>(what, grid->name(), 1, MonolithicAlloc::NotRequired);
    } else {
      add_tracer<Updated>(what, grid, kg/kg, 1);
    }
  }

protected:
  void initialize_impl (const RunType) override {}
  void run_impl (const double) override {}
  void finalize_impl () override {}
};

// Register it before main, so it is in the factory when the driver creates the processes
namespace {
const bool wtt_dummy_registered = [] {
  AtmosphereProcessFactory::instance().register_product("wtt_dummy_modifier",
                                                        &create_atmosphere_process<WttDummyModifier>);
  return true;
}();
} // anonymous namespace

} // namespace scream

#include "eamxx_ad_test.cpp"
