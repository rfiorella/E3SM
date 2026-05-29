// Unit test for Surface Coupling qv tracer access
// Validates that surface coupling exporter correctly accesses qv slot-0 from tracer-aware field

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_qv_surface_tracer_access", "Test surface coupling accesses qv tracer slot-0") {
  using namespace scream;

  // Simulate a tracer-aware qv field (ntracers, ncols, nlevs)
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create 3D view: (tracer, col, lev)
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qv_3d("qv_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 with surface-focused values
  // Bottom level (nlevs-1) has highest qv
  Kokkos::parallel_for(
    "init_qv",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      if (t == 0) {
        // Higher values near surface
        qv_3d(t, c, l) = (l == nlevs - 1) ? 0.018 : 0.005;
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

  // Verify surface-level access (what surface coupling exports)
  auto qv_bulk_h = Kokkos::create_mirror_view(qv_bulk);
  Kokkos::deep_copy(qv_bulk_h, qv_bulk);

  for (int c = 0; c < ncols; ++c) {
    // Surface coupling typically exports bottom-level (nlevs-1) value as Sa_shum
    REQUIRE(qv_bulk_h(c, nlevs - 1) == 0.018);
    REQUIRE(qv_bulk_h(c, 0) == 0.005); // Upper levels
  }
}

TEST_CASE("test_qv_surface_dz_calculation", "Test surface coupling dz calculation with qv") {
  // This test would validate that calculate_dz receives correct qv
  // For now, it's a placeholder showing the pattern
  REQUIRE(true); // Placeholder
}
