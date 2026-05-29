// Unit test for cloud tracer (qc, qi, qm) access
// Validates that processes correctly access cloud tracer slot-0 from tracer-aware fields

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_cloud_tracer_access_qc", "Test qc tracer slot-0 access") {
  using namespace scream;

  // Simulate a tracer-aware qc field (ntracers, ncols, nlevs)
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create 3D view: (tracer, col, lev)
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qc_3d("qc_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 = 1.0e-5, slot-1 = 2.0e-5, slot-2 = 3.0e-5
  Kokkos::parallel_for(
    "init_qc",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qc_3d(t, c, l) = static_cast<Real>((t + 1) * 1.0e-5);
    }
  );
  Kokkos::fence();

  // Extract bulk water (slot-0) via SUBVIEW accessor
  auto qc_bulk = get_tracer_bulk_subview(qc_3d);

  // Verify dimensions: should be (ncols, nlevs)
  REQUIRE(qc_bulk.extent(0) == ncols);
  REQUIRE(qc_bulk.extent(1) == nlevs);

  // Verify values: all should be 1.0e-5 (slot-0)
  auto qc_bulk_h = Kokkos::create_mirror_view(qc_bulk);
  Kokkos::deep_copy(qc_bulk_h, qc_bulk);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qc_bulk_h(c, l) == 1.0e-5);
    }
  }
}

TEST_CASE("test_cloud_tracer_access_qi", "Test qi tracer slot-0 access") {
  using namespace scream;

  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qi_3d("qi_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 = 2.0e-5, slot-1 = 4.0e-5, slot-2 = 6.0e-5
  Kokkos::parallel_for(
    "init_qi",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qi_3d(t, c, l) = static_cast<Real>((t + 1) * 2.0e-5);
    }
  );
  Kokkos::fence();

  auto qi_bulk = get_tracer_bulk_subview(qi_3d);

  REQUIRE(qi_bulk.extent(0) == ncols);
  REQUIRE(qi_bulk.extent(1) == nlevs);

  auto qi_bulk_h = Kokkos::create_mirror_view(qi_bulk);
  Kokkos::deep_copy(qi_bulk_h, qi_bulk);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qi_bulk_h(c, l) == 2.0e-5);
    }
  }
}

TEST_CASE("test_cloud_tracer_access_qm", "Test qm tracer slot-0 access") {
  using namespace scream;

  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qm_3d("qm_3d", ntracers, ncols, nlevs);

  // Initialize: slot-0 = 3.0e-5, slot-1 = 6.0e-5, slot-2 = 9.0e-5
  Kokkos::parallel_for(
    "init_qm",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qm_3d(t, c, l) = static_cast<Real>((t + 1) * 3.0e-5);
    }
  );
  Kokkos::fence();

  auto qm_bulk = get_tracer_bulk_subview(qm_3d);

  REQUIRE(qm_bulk.extent(0) == ncols);
  REQUIRE(qm_bulk.extent(1) == nlevs);

  auto qm_bulk_h = Kokkos::create_mirror_view(qm_bulk);
  Kokkos::deep_copy(qm_bulk_h, qm_bulk);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qm_bulk_h(c, l) == 3.0e-5);
    }
  }
}

TEST_CASE("test_cloud_tracer_all_species", "Test all cloud species access patterns") {
  using namespace scream;

  constexpr int ntracers = 2;
  constexpr int ncols = 5;
  constexpr int nlevs = 32;

  // Create all three cloud species
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qc_3d("qc_3d", ntracers, ncols, nlevs);
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qi_3d("qi_3d", ntracers, ncols, nlevs);
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qm_3d("qm_3d", ntracers, ncols, nlevs);

  // Initialize with unique values for each species
  Kokkos::parallel_for(
    "init_all",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qc_3d(t, c, l) = 1.0e-4 + t * 0.1e-4;
      qi_3d(t, c, l) = 2.0e-4 + t * 0.2e-4;
      qm_3d(t, c, l) = 3.0e-4 + t * 0.3e-4;
    }
  );
  Kokkos::fence();

  // Extract bulk water for all species
  auto qc_bulk = get_tracer_bulk_subview(qc_3d);
  auto qi_bulk = get_tracer_bulk_subview(qi_3d);
  auto qm_bulk = get_tracer_bulk_subview(qm_3d);

  // Verify dimensions
  REQUIRE(qc_bulk.extent(0) == ncols);
  REQUIRE(qi_bulk.extent(0) == ncols);
  REQUIRE(qm_bulk.extent(0) == ncols);

  // Verify slot-0 values
  auto qc_bulk_h = Kokkos::create_mirror_view(qc_bulk);
  auto qi_bulk_h = Kokkos::create_mirror_view(qi_bulk);
  auto qm_bulk_h = Kokkos::create_mirror_view(qm_bulk);
  Kokkos::deep_copy(qc_bulk_h, qc_bulk);
  Kokkos::deep_copy(qi_bulk_h, qi_bulk);
  Kokkos::deep_copy(qm_bulk_h, qm_bulk);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      REQUIRE(qc_bulk_h(c, l) == 1.0e-4);
      REQUIRE(qi_bulk_h(c, l) == 2.0e-4);
      REQUIRE(qm_bulk_h(c, l) == 3.0e-4);
    }
  }
}
