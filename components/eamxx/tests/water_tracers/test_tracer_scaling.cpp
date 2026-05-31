#include "catch2/catch.hpp"
#include "water_tracer_ratio.hpp"
#include "share/scream_types.hpp"
#include <Kokkos_Core.hpp>

namespace scream {
namespace water_tracers {

TEST_CASE("compute_tracer_ratio", "[water_tracers]") {
  using scream::Real;

  SECTION("basic_ratio") {
    Real bulk = 10.0;
    Real tracer = 5.0;
    Real ratio = compute_tracer_ratio(tracer, bulk);
    REQUIRE(ratio == Approx(0.5).epsilon(1e-15));
  }

  SECTION("zero_bulk") {
    Real bulk = 0.0;
    Real tracer = 5.0;
    Real ratio = compute_tracer_ratio(tracer, bulk);
    REQUIRE(ratio == 0.0);
  }

  SECTION("below_threshold") {
    Real bulk = 1e-21;  // Below default threshold of 1e-20
    Real tracer = 5e-21;
    Real ratio = compute_tracer_ratio(tracer, bulk);
    REQUIRE(ratio == 0.0);
  }

  SECTION("custom_threshold") {
    Real bulk = 1e-10;
    Real tracer = 5e-11;
    Real min_bulk = 1e-9;  // Bulk below custom threshold
    Real ratio = compute_tracer_ratio(tracer, bulk, min_bulk);
    REQUIRE(ratio == 0.0);
  }

  SECTION("machine_precision") {
    Real bulk = 1.0;
    Real tracer = 0.5;
    Real ratio = compute_tracer_ratio(tracer, bulk);
    REQUIRE(std::abs(ratio - 0.5) < 1e-15);
  }
}

TEST_CASE("compute_tracer_ratio_field", "[water_tracers]") {
  using scream::Real;

  const int ntracers = 2;
  const int ncols = 4;
  const int nlevs = 72;

  // Create test fields
  Kokkos::View<Real***> tracer_field("tracer_field", ntracers, ncols, nlevs);
  Kokkos::View<Real**> ratio_field("ratio_field", ncols, nlevs);

  // Initialize: tracer 0 = 10.0, tracer 1 = 5.0 (ratio = 0.5)
  Kokkos::parallel_for(
    "init_fields",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(const int icol, const int ilev) {
      tracer_field(0, icol, ilev) = 10.0;
      tracer_field(1, icol, ilev) = 5.0;
    }
  );
  Kokkos::fence();

  // Compute ratios
  compute_tracer_ratio_field(tracer_field, 1, ratio_field);
  Kokkos::fence();

  // Verify results on host
  auto ratio_host = Kokkos::create_mirror_view(ratio_field);
  Kokkos::deep_copy(ratio_host, ratio_field);

  for (int icol = 0; icol < ncols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      REQUIRE(ratio_host(icol, ilev) == Approx(0.5).epsilon(1e-15));
    }
  }
}

TEST_CASE("compute_mean_tracer_ratio", "[water_tracers]") {
  using scream::Real;

  const int ntracers = 2;
  const int ncols = 4;
  const int nlevs = 72;

  Kokkos::View<Real***> tracer_field("tracer_field", ntracers, ncols, nlevs);

  SECTION("uniform_ratio") {
    // Initialize: tracer 0 = 10.0, tracer 1 = 5.0 everywhere
    Kokkos::parallel_for(
      "init_uniform",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
      KOKKOS_LAMBDA(const int icol, const int ilev) {
        tracer_field(0, icol, ilev) = 10.0;
        tracer_field(1, icol, ilev) = 5.0;
      }
    );
    Kokkos::fence();

    Real mean_ratio = compute_mean_tracer_ratio(tracer_field, 1);
    REQUIRE(mean_ratio == Approx(0.5).epsilon(1e-15));
  }

  SECTION("mixed_ratios") {
    // Initialize: half points with ratio 0.5, half with ratio 0.25
    Kokkos::parallel_for(
      "init_mixed",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
      KOKKOS_LAMBDA(const int icol, const int ilev) {
        if (icol < ncols / 2) {
          tracer_field(0, icol, ilev) = 10.0;
          tracer_field(1, icol, ilev) = 5.0;  // ratio = 0.5
        } else {
          tracer_field(0, icol, ilev) = 10.0;
          tracer_field(1, icol, ilev) = 2.5;  // ratio = 0.25
        }
      }
    );
    Kokkos::fence();

    Real mean_ratio = compute_mean_tracer_ratio(tracer_field, 1);
    Real expected = 0.375;  // Average of 0.5 and 0.25
    REQUIRE(mean_ratio == Approx(expected).epsilon(1e-14));
  }

  SECTION("with_zeros") {
    // Initialize: some points below threshold
    Kokkos::parallel_for(
      "init_with_zeros",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
      KOKKOS_LAMBDA(const int icol, const int ilev) {
        if (ilev < nlevs / 2) {
          tracer_field(0, icol, ilev) = 10.0;
          tracer_field(1, icol, ilev) = 5.0;  // ratio = 0.5
        } else {
          tracer_field(0, icol, ilev) = 0.0;  // Below threshold
          tracer_field(1, icol, ilev) = 0.0;
        }
      }
    );
    Kokkos::fence();

    Real mean_ratio = compute_mean_tracer_ratio(tracer_field, 1);
    // Should average only over points with bulk > threshold
    REQUIRE(mean_ratio == Approx(0.5).epsilon(1e-15));
  }
}

} // namespace water_tracers
} // namespace scream
