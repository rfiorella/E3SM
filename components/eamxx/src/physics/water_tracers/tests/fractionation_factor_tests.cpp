#include "catch2/catch.hpp"
#include "../water_isotopes.hpp"
#include "../water_isotope_fractionation.hpp"

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
            {"HDO", 0.99308},
            {"H218O", 0.99222},
            {"H217O", 0.99598},
            {"HTO", 0.99082}
        }},
        {0.1, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.99444},
            {"H218O", 0.99374},
            {"H217O", 0.99677},
            {"HTO", 0.99262}
        }},
        {0.2, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.99475},
            {"H218O", 0.99409},
            {"H217O", 0.99695},
            {"HTO", 0.99303}
        }},
        {0.3, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.99612},
            {"H218O", 0.99563},
            {"H217O", 0.99775},
            {"HTO", 0.99485}
        }},
        {0.5, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.99450},
            {"H218O", 0.99381},
            {"H217O", 0.99680},
            {"HTO", 0.99270}
        }},
        {1.0, {
            {"H2O", 1.00000},
            {"H216O", 1.00000},
            {"HDO", 0.99202},
            {"H218O", 0.99102},
            {"H217O", 0.99535},
            {"HTO", 0.98942}
        }}
    };

    for (double U : ufric) {
        DYNAMIC_SECTION("Friction velocity = " << U) {
            for (const auto& species : scream::WaterIsotopes::WaterIsotopologueNames) {
                UNSCOPED_INFO("Testing species: " << species << " at u* = " << U);
                double result = scream::WaterIsotopes::AlphaKMol(species, 1.225, 10., U);
                REQUIRE_THAT(result, Catch::Matchers::WithinRel(expected_values[U][species], 1e-5));
            }
        }
    }
}

}
