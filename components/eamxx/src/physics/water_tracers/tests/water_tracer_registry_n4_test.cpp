// Unit test for water tracer registry with 4-tracer configuration
// Requires build with registry_n4.cmake tracer config

#include "catch2/catch.hpp"
#include "physics/water_tracers/water_tracer_registry.hpp"
#include "physics/water_tracers/water_isotopes.hpp"

namespace {

using namespace scream::WaterTracers;
using namespace scream::WaterIsotopes;

TEST_CASE("water_tracer_registry_n4", "[water_tracers]") {

  SECTION("registry_size") {
    // registry_n4.cmake defines 4 tracers
    REQUIRE(WTRC_MAX_CNST == 4);
  }

  SECTION("tracer_0_bulk_H2O") {
    // Tracer 0: bulk H2O (catalog index 0, not a tag)
    REQUIRE(tracer_isotopologue(0) == 0);
    REQUIRE(tracer_name(0) == "bulk");
    REQUIRE_FALSE(tracer_is_tag(0));
  }

  SECTION("tracer_1_passive_H2O_tag") {
    // Tracer 1: passive H2O tag (catalog index 0, is a tag)
    REQUIRE(tracer_isotopologue(1) == 0);
    REQUIRE(tracer_name(1) == "passive");
    REQUIRE(tracer_is_tag(1));
  }

  SECTION("tracer_2_HDO") {
    // Tracer 2: HDO (catalog index 2, not a tag)
    REQUIRE(tracer_isotopologue(2) == 2);
    REQUIRE(tracer_name(2) == "hdo");
    REQUIRE_FALSE(tracer_is_tag(2));
  }

  SECTION("tracer_3_H218O") {
    // Tracer 3: H218O (catalog index 3, not a tag)
    REQUIRE(tracer_isotopologue(3) == 3);
    REQUIRE(tracer_name(3) == "h218o");
    REQUIRE_FALSE(tracer_is_tag(3));
  }

  SECTION("per_tracer_names_unique") {
    // All tracer names must be unique
    std::set<std::string> names;
    for (int i = 0; i < WTRC_MAX_CNST; ++i) {
      auto name = std::string(tracer_name(i));
      REQUIRE(names.find(name) == names.end());
      names.insert(name);
    }
    REQUIRE(names.size() == 4);
  }

  SECTION("find_tracer_by_name_roundtrip") {
    // All configured names should be findable
    auto bulk_idx = find_tracer_by_name("bulk");
    REQUIRE(bulk_idx.has_value());
    REQUIRE(bulk_idx.value() == 0);

    auto passive_idx = find_tracer_by_name("passive");
    REQUIRE(passive_idx.has_value());
    REQUIRE(passive_idx.value() == 1);

    auto hdo_idx = find_tracer_by_name("hdo");
    REQUIRE(hdo_idx.has_value());
    REQUIRE(hdo_idx.value() == 2);

    auto h218o_idx = find_tracer_by_name("h218o");
    REQUIRE(h218o_idx.has_value());
    REQUIRE(h218o_idx.value() == 3);

    // Non-existent name should return nullopt
    auto missing = find_tracer_by_name("not_a_tracer");
    REQUIRE_FALSE(missing.has_value());
  }

  SECTION("catalog_index_to_isotopologue_name") {
    // Verify catalog indices map to expected isotopologue names
    REQUIRE(WaterIsotopologueNames[tracer_isotopologue(0)] == "H2O");
    REQUIRE(WaterIsotopologueNames[tracer_isotopologue(1)] == "H2O");
    REQUIRE(WaterIsotopologueNames[tracer_isotopologue(2)] == "HDO");
    REQUIRE(WaterIsotopologueNames[tracer_isotopologue(3)] == "H218O");
  }

  SECTION("device_callable") {
    // Verify registry queries are constexpr and device-callable
    // This section compiles = proof that functions are constexpr/device-safe
    constexpr int bulk_cat = tracer_isotopologue(0);
    constexpr auto bulk_name = tracer_name(0);
    constexpr bool bulk_is_tag = tracer_is_tag(0);

    REQUIRE(bulk_cat == 0);
    REQUIRE(bulk_name == "bulk");
    REQUIRE_FALSE(bulk_is_tag);
  }
}

} // anonymous namespace
