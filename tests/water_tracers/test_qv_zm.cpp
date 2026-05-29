// Unit test for ZM qv tracer access
// Validates that ZM correctly accesses and modifies qv slot-0 from tracer-aware field

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_qv_zm_tracer_access", "Test ZM accesses qv tracer slot-0") {
  using namespace scream;

  // Simulate a tracer-aware qv field (ntracers, ncols, nlevs)
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create 3D view: (tracer, col, lev)
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qv_3d("qv_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 = 0.01 kg/kg (10 g/kg)
  Kokkos::parallel_for(
    "init_qv",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qv_3d(t, c, l) = (t == 0) ? 0.01 : 0.0;
    }
  );
  Kokkos::fence();

  // Extract bulk water (slot-0) via SUBVIEW accessor
  auto qv_bulk = get_tracer_bulk_subview(qv_3d);

  // Verify dimensions
  REQUIRE(qv_bulk.extent(0) == ncols);
  REQUIRE(qv_bulk.extent(1) == nlevs);

  // Simulate ZM tendency application: add 0.001 kg/kg
  const Real qv_tend = 0.001;
  Kokkos::parallel_for(
    "apply_zm_tend",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(int c, int l) {
      qv_bulk(c, l) += qv_tend;
    }
  );
  Kokkos::fence();

  // Verify: qv_3d slot-0 should be updated, other slots unchanged
  auto qv_3d_h = Kokkos::create_mirror_view(qv_3d);
  Kokkos::deep_copy(qv_3d_h, qv_3d);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qv_3d_h(0, c, l) == 0.011); // slot-0 updated
      REQUIRE(qv_3d_h(1, c, l) == 0.0);   // slot-1 unchanged
      REQUIRE(qv_3d_h(2, c, l) == 0.0);   // slot-2 unchanged
    }
  }
}
