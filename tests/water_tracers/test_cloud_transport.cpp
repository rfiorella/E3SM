// Component test for cloud tracer transport
// Validates that cloud tracer mass is conserved during transport operations

#include "catch2/catch.hpp"
#include "share/field/field_tracer_access.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("test_cloud_transport_conservation", "Test cloud tracer mass conservation") {
  using namespace scream;

  // Simulate a simple transport: move cloud water from one level to another
  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create qc field
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qc_3d("qc_3d", ntracers, ncols, nlevs);

  // Initialize: tracer 0 (bulk) = 1.0e-4, tracer 1 = 0.5e-4 (half of bulk), tracer 2 = 0.25e-4
  Kokkos::parallel_for(
    "init_qc",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      if (t == 0) {
        qc_3d(t, c, l) = 1.0e-4;
      } else if (t == 1) {
        qc_3d(t, c, l) = 0.5e-4;
      } else {
        qc_3d(t, c, l) = 0.25e-4;
      }
    }
  );
  Kokkos::fence();

  // Compute initial total mass for each tracer
  Kokkos::View<Real*, Kokkos::DefaultExecutionSpace> initial_mass("initial_mass", ntracers);
  Kokkos::parallel_for(
    "compute_initial_mass",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Kokkos::atomic_add(&initial_mass(t), qc_3d(t, c, l));
    }
  );
  Kokkos::fence();

  // Simulate a transport operation: proportional transport
  // Move 10% of cloud water from level l to level l+1 (where l+1 exists)
  Kokkos::parallel_for(
    "transport",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs - 1}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Real transport_amount = 0.1 * qc_3d(t, c, l);
      qc_3d(t, c, l) -= transport_amount;
      qc_3d(t, c, l + 1) += transport_amount;
    }
  );
  Kokkos::fence();

  // Compute final total mass for each tracer
  Kokkos::View<Real*, Kokkos::DefaultExecutionSpace> final_mass("final_mass", ntracers);
  Kokkos::parallel_for(
    "compute_final_mass",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Kokkos::atomic_add(&final_mass(t), qc_3d(t, c, l));
    }
  );
  Kokkos::fence();

  // Verify mass conservation: initial_mass == final_mass for each tracer
  auto initial_mass_h = Kokkos::create_mirror_view(initial_mass);
  auto final_mass_h = Kokkos::create_mirror_view(final_mass);
  Kokkos::deep_copy(initial_mass_h, initial_mass);
  Kokkos::deep_copy(final_mass_h, final_mass);

  for (int t = 0; t < ntracers; ++t) {
    // Mass should be conserved to machine precision
    REQUIRE(std::abs(initial_mass_h(t) - final_mass_h(t)) < 1.0e-14 * initial_mass_h(t));
  }
}

TEST_CASE("test_cloud_tracer_ratio_preservation", "Test tracer ratios preserved during transport") {
  using namespace scream;

  constexpr int ntracers = 3;
  constexpr int ncols = 10;
  constexpr int nlevs = 72;

  // Create qc field with known tracer ratios
  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qc_3d("qc_3d", ntracers, ncols, nlevs);

  // Initialize with specific ratios: tracer_1 / tracer_0 = 0.5, tracer_2 / tracer_0 = 0.25
  Kokkos::parallel_for(
    "init_qc",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Real base_value = 1.0e-4 * (1 + c * 0.1) * (1 + l * 0.01);  // Vary by column and level
      if (t == 0) {
        qc_3d(t, c, l) = base_value;
      } else if (t == 1) {
        qc_3d(t, c, l) = 0.5 * base_value;
      } else {
        qc_3d(t, c, l) = 0.25 * base_value;
      }
    }
  );
  Kokkos::fence();

  // Simulate proportional transport (preserves ratios)
  // Move 20% from level l to level l+1
  Kokkos::parallel_for(
    "transport",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs - 1}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Real transport_amount = 0.2 * qc_3d(t, c, l);
      qc_3d(t, c, l) -= transport_amount;
      qc_3d(t, c, l + 1) += transport_amount;
    }
  );
  Kokkos::fence();

  // Verify tracer ratios are preserved
  auto qc_3d_h = Kokkos::create_mirror_view(qc_3d);
  Kokkos::deep_copy(qc_3d_h, qc_3d);

  for (int c = 0; c < ncols; ++c) {
    for (int l = 0; l < nlevs; ++l) {
      Real bulk = qc_3d_h(0, c, l);
      if (bulk > 1.0e-20) {  // Avoid division by zero
        Real ratio_1 = qc_3d_h(1, c, l) / bulk;
        Real ratio_2 = qc_3d_h(2, c, l) / bulk;

        // Ratios should be preserved: 0.5 and 0.25
        REQUIRE(std::abs(ratio_1 - 0.5) < 1.0e-12);
        REQUIRE(std::abs(ratio_2 - 0.25) < 1.0e-12);
      }
    }
  }
}

TEST_CASE("test_cloud_qi_transport", "Test qi tracer transport") {
  using namespace scream;

  constexpr int ntracers = 2;
  constexpr int ncols = 5;
  constexpr int nlevs = 32;

  Kokkos::View<Real***, Kokkos::DefaultExecutionSpace> qi_3d("qi_3d", ntracers, ncols, nlevs);

  // Initialize with uniform distribution
  Kokkos::parallel_for(
    "init_qi",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      qi_3d(t, c, l) = (t == 0) ? 2.0e-4 : 1.0e-4;
    }
  );
  Kokkos::fence();

  // Compute initial total for each tracer
  Kokkos::View<Real*, Kokkos::DefaultExecutionSpace> initial_total("initial_total", ntracers);
  Kokkos::parallel_for(
    "sum_initial",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Kokkos::atomic_add(&initial_total(t), qi_3d(t, c, l));
    }
  );
  Kokkos::fence();

  // Simulate sedimentation: move 5% downward
  Kokkos::parallel_for(
    "sediment",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs - 1}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Real sed_amount = 0.05 * qi_3d(t, c, l);
      qi_3d(t, c, l) -= sed_amount;
      qi_3d(t, c, l + 1) += sed_amount;
    }
  );
  Kokkos::fence();

  // Compute final total
  Kokkos::View<Real*, Kokkos::DefaultExecutionSpace> final_total("final_total", ntracers);
  Kokkos::parallel_for(
    "sum_final",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ntracers, ncols, nlevs}),
    KOKKOS_LAMBDA(int t, int c, int l) {
      Kokkos::atomic_add(&final_total(t), qi_3d(t, c, l));
    }
  );
  Kokkos::fence();

  // Verify conservation
  auto initial_h = Kokkos::create_mirror_view(initial_total);
  auto final_h = Kokkos::create_mirror_view(final_total);
  Kokkos::deep_copy(initial_h, initial_total);
  Kokkos::deep_copy(final_h, final_total);

  for (int t = 0; t < ntracers; ++t) {
    REQUIRE(std::abs(initial_h(t) - final_h(t)) < 1.0e-14 * initial_h(t));
  }
}
