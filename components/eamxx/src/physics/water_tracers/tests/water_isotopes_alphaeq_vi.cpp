#include "catch2/catch.hpp"
#include "../water_isotopes.hpp"
#include "../water_isotope_fractionation.cpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include <string>

namespace {

/* test fractionation factor routines
*/

TEST_CASE("AlphaEqIceVapor computes expected values for all isotopologues", "[AlphaEqIceVapor]") {
    double T = 253.15; // -20°C as a test temperature

    //define expected values
    std::map<std::string, double> expected_values = {
        {"H2O", 1.00000},
        {"H216O", 1.00000},
        {"HDO", 1.17313},
        {"H218O", 1.01871},
        {"H217O", 1.00983},
        {"HTO", 1.37624}
    };

    for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
        UNSCOPED_INFO("Testing species: " << species);
        double result = scream::WaterIsotopes::AlphaEqIceVapor(species, T);
        REQUIRE_THAT(result, Catch::Matchers::WithinRel(expected_values[species],1e-5));
    }
}

TEST_CASE("AlphaEqLiquidVapor computes expected values for all isotopologues", "[AlphaEqLiquidVapor]") {
    double T = 293.15; // 20°C as a test temperature

    //define expected values
    std::map<std::string, double> expected_values = {
        {"H2O", 1.00000},
        {"H216O", 1.00000},
        {"HDO", 1.08435},
        {"H218O", 1.00977},
        {"H217O", 1.00515},
        {"HTO", 1.17582}
    };

    for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
        UNSCOPED_INFO("Testing species: " << species);
        double result = scream::WaterIsotopes::AlphaEqLiquidVapor(species, T);
        REQUIRE_THAT(result, Catch::Matchers::WithinRel(expected_values[species],1e-5));
    }
}

}
