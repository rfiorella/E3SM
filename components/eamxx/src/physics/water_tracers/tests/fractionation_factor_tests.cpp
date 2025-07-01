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

/* Test kmol function across range of u* values */

TEST_CASE("AlphaKMol computes kinetic fractionations across wind speeds", "[AlphaKMol]") {
    std::vector<double> ufric = {0.01, 0.1, 0.2, 0.3, 0.5, 1.0};

    // Map of expected values per temperature per species
    std::map<double, std::map<std::string, double>> expected_values = {
        {0.01, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.98411},
            {"H218O", 0.98214},
            {"H217O", 0.99073},
            {"HTO", 0.97898}
        }},
        {0.1, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.98427},
            {"H218O", 0.98232},
            {"H217O", 0.99083},
            {"HTO", 0.97919}
        }},
        {0.2, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.98432},
            {"H218O", 0.98237},
            {"H217O", 0.99085},
            {"HTO", 0.97925}
        }},
        {0.3, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.96729},
            {"H218O", 0.96332},
            {"H217O", 0.98077},
            {"HTO", 0.95700}
        }},
        {0.5, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.97851},
            {"H218O", 0.97586},
            {"H217O", 0.98742},
            {"HTO", 0.97164}
        }},
        {1.0, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.98348},
            {"H218O", 0.98144},
            {"H217O", 0.99036},
            {"HTO", 0.97817}
        }}
    };

    for (double U : ufric) {
        DYNAMIC_SECTION("Friction velocity = " << U) {
            for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
                UNSCOPED_INFO("Testing species: " << species << " at u* = " << U);
                double result = scream::WaterIsotopes::AlphaKMol(species, U);
                REQUIRE_THAT(result, Catch::Matchers::WithinRel(expected_values[U][species], 1e-5));
            }
        }
    }
}

}
