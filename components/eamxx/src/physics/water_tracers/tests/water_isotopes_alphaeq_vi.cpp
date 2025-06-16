#include "catch2/catch.hpp"
#include "../water_isotopes.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include <string>

namespace {

/* test fractionation factor routines
*/

TEST_CASE("AlphaEqIceVapor computes expected values for all isotopologues", "[AlphaEqIceVapor]") {
    double T = 253.15; // -20°C as a test temperature

    for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
        UNSCOPED_INFO("Testing species: " << species);
        double result = scream::WaterIsotopes::AlphaEqIceVapor(species, T);
        // Example condition - this needs real expected values or ranges.
        REQUIRE(result >= 1.0);
        REQUIRE(result < 1.5);
    }
}

TEST_CASE("AlphaEqLiquidVapor computes expected values for all isotopologues", "[AlphaEqLiquidVapor]") {
    double T = 293.15; // 20°C as a test temperature

    for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
        UNSCOPED_INFO("Testing species: " << species);
        double result = scream::WaterIsotopes::AlphaEqLiquidVapor(species, T);
        // Example condition - this needs real expected values or ranges.
        REQUIRE(result >= 1.0);
        REQUIRE(result < 1.5);
    }
}

}
