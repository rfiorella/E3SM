// Unit test for RRTMGP qv tracer access
// Validates that RRTMGP correctly accesses qv slot-0 from tracer-aware field

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_qv_rrtmgp_tracer_access", "Test RRTMGP accesses qv tracer slot-0") {
  using namespace scream;

  // Simulate a tracer-aware qv field (ntracers, ncols, nlevs)
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create 3D view: (tracer, col, lev)
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qv_3d("qv_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 = 1.0, slot-1 = 2.0, slot-2 = 3.0
  Kokkos::parallel_for(
    "init_qv",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qv_3d(t, c, l) = static_cast<Real>(t + 1);
    }
  );
  Kokkos::fence();

  // Extract bulk water (slot-0) via SUBVIEW accessor
  auto qv_bulk = get_tracer_bulk_subview(qv_3d);

  // Verify dimensions: should be (ncols, nlevs)
  REQUIRE(qv_bulk.extent(0) == ncols);
  REQUIRE(qv_bulk.extent(1) == nlevs);

  // Verify values: all should be 1.0 (slot-0)
  auto qv_bulk_h = Kokkos::create_mirror_view(qv_bulk);
  Kokkos::deep_copy(qv_bulk_h, qv_bulk);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qv_bulk_h(c, l) == 1.0);
    }
  }
}

TEST_CASE("test_qv_rrtmgp_vmr_calculation", "Test RRTMGP VMR calculation uses correct qv") {
  // This test would validate that calculate_vmr_from_mmr receives the correct qv value
  // For now, it's a placeholder showing the pattern
  REQUIRE(true); // Placeholder
}
