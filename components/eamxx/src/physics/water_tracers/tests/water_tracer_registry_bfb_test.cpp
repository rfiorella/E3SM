// BFB test: bulk-only (N=1) with SCREAM_TRACE_WATER=ON vs OFF
// When configured with only bulk H2O, the ON build should produce
// bit-for-bit identical results to the OFF build

#include "catch2/catch.hpp"
#include "physics/water_tracers/water_tracer_registry.hpp"

namespace {

using namespace scream::WaterTracers;

TEST_CASE("water_tracer_registry_bulk_vs_default_bfb", "[water_tracers]") {

  SECTION("single_bulk_tracer") {
    // When built with bulk_only.cmake, WTRC_MAX_CNST should be 1
    REQUIRE(WTRC_MAX_CNST == 1);

    // Tracer 0 should be bulk H2O
    REQUIRE(tracer_isotopologue(0) == 0);  // catalog index 0 = H2O
    REQUIRE(tracer_name(0) == "bulk");
    REQUIRE_FALSE(tracer_is_tag(0));
  }

  SECTION("registry_lookup") {
    auto bulk_idx = find_tracer_by_name("bulk");
    REQUIRE(bulk_idx.has_value());
    REQUIRE(bulk_idx.value() == 0);
  }

  // Note: This test verifies the registry structure for bulk-only config.
  // The actual BFB comparison of physics results is done via a separate
  // test harness that runs the same physics problem with both builds and
  // compares outputs. This test ensures the registry metadata is correct.
}

} // anonymous namespace
