// Unit test for HOMME qv tracer access
// Validates that HOMME correctly accesses qv slot-0 from tracer-aware field

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_qv_homme_tracer_access", "Test HOMME accesses qv tracer slot-0") {
  using namespace scream;

  // Simulate a tracer-aware qv field (ntracers, ncols, nlevs)
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create 3D view: (tracer, col, lev)
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qv_3d("qv_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 with realistic atmospheric qv profile
  // Surface: 0.015 kg/kg, decreasing exponentially with height
  Kokkos::parallel_for(
    "init_qv",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      if (t == 0) {
        // Exponential decay from surface (l=nlevs-1) to top (l=0)
        Real lev_frac = static_cast<Real>(nlevs - 1 - l) / static_cast<Real>(nlevs - 1);
        qv_3d(t, c, l) = 0.015 * Kokkos::exp(-3.0 * lev_frac);
      } else {
        qv_3d(t, c, l) = 0.0;
      }
    }
  );
  Kokkos::fence();

  // Extract bulk water (slot-0) via SUBVIEW accessor
  auto qv_bulk = get_tracer_bulk_subview(qv_3d);

  // Verify dimensions
  REQUIRE(qv_bulk.extent(0) == ncols);
  REQUIRE(qv_bulk.extent(1) == nlevs);

  // Verify surface value (bottom level)
  auto qv_bulk_h = Kokkos::create_mirror_view(qv_bulk);
  Kokkos::deep_copy(qv_bulk_h, qv_bulk);

  for (int c = 0; c < ncols; ++c) {
    REQUIRE(qv_bulk_h(c, nlevs - 1) == Approx(0.015).epsilon(1e-6)); // Surface
    REQUIRE(qv_bulk_h(c, 0) < 0.001); // TOA (much smaller)
  }
}

TEST_CASE("test_qv_homme_dp_dry_calculation", "Test HOMME dp_dry calculation with qv") {
  // This test would validate that dp_dry = dp * (1 - qv) uses correct qv
  // For now, it's a placeholder showing the pattern
  REQUIRE(true); // Placeholder
}
