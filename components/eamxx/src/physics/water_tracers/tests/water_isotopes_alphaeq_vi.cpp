#include "catch2/catch.hpp"
#include "water_isotopes.hpp"

namespace {

TEST_CASE("AlphaEqIceVapor computes expected values for all isotopologues", "[AlphaEqIceVapor]") {
    double T = 253.15; // 0°C as a test temperature

    for (const auto& species : WaterIsotopologues::isoname) {
        INFO("Testing species: " << species);
        double result = AlphaEqIceVapor(species, T);
        
        // Example condition - this needs real expected values or ranges.
        REQUIRE(result > 1.0);
        REQUIRE(result < 1.5);
    }
}

}

