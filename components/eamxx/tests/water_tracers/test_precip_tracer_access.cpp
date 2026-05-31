/*
 * Test precipitation tracer field access patterns
 *
 * Validates extension of P3 precipitation fields to tracer dimension:
 * - qr (rain mixing ratio): (col, lev) → (tracer, col, lev)
 * - precip_liq_surf_mass: (col) → (tracer, col)
 * - precip_ice_surf_mass: (col) → (tracer, col)
 *
 * Tests subview access pattern for extracting slot 0 (bulk water).
 */

#include "catch2/catch.hpp"

#include "share/field/field_tracer_access.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

using Real = double;

// Kernel for 3D field (qr) using subview
template<typename ViewSubT, typename ViewOutT>
struct QrSubviewKernel {
  ViewSubT qr_bulk;     // 2D subview of slot 0 (col, lev)
  ViewOutT output;      // (col, lev)

  KOKKOS_INLINE_FUNCTION
  void operator()(int icol, int ilev) const {
    // Simple sedimentation-like calculation
    output(icol, ilev) = qr_bulk(icol, ilev) * 0.95;  // 5% falls out
  }
};

// Kernel for 2D field (surface precip) using subview
template<typename ViewSubT, typename ViewOutT>
struct SurfPrecipSubviewKernel {
  ViewSubT precip_bulk; // 1D subview of slot 0 (col)
  ViewOutT output;      // (col)

  KOKKOS_INLINE_FUNCTION
  void operator()(int icol) const {
    // Accumulate precipitation
    output(icol) = precip_bulk(icol) + 0.1;
  }
};

TEST_CASE("precip_3d_tracer_access") {
  const int ncols = 10;
  const int nlevs = 72;
  const int ntracers = 1;  // Only bulk water

  // Create 3D tracer field for qr (tracer, col, lev)
  using view_3d = Kokkos::View<Real***>;
  view_3d qr_tracer("qr_tracer", ntracers, ncols, nlevs);

  // Output view
  using view_2d = Kokkos::View<Real**>;
  view_2d output("output", ncols, nlevs);

  // Initialize qr with typical rain mixing ratio values
  auto init_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs});

  Kokkos::parallel_for("init_qr", init_policy,
    KOKKOS_LAMBDA(int icol, int ilev) {
      Real height_frac = static_cast<Real>(ilev) / nlevs;
      // Rain typically in lower atmosphere
      qr_tracer(0, icol, ilev) = 0.001 * Kokkos::exp(-(1.0 - height_frac) * 3.0);
    });

  Kokkos::fence();

  SECTION("qr_subview_access") {
    // Extract slot 0 using subview (as done in P3 interface)
    auto qr_bulk = Kokkos::subview(qr_tracer, 0, Kokkos::ALL(), Kokkos::ALL());

    // Run kernel with subview
    Kokkos::parallel_for("qr_kernel", init_policy,
      QrSubviewKernel<decltype(qr_bulk), view_2d>{qr_bulk, output});

    Kokkos::fence();

    // Verify output
    auto h_input = Kokkos::create_mirror_view(qr_tracer);
    auto h_output = Kokkos::create_mirror_view(output);
    Kokkos::deep_copy(h_input, qr_tracer);
    Kokkos::deep_copy(h_output, output);

    bool correct = true;
    Real max_error = 0.0;

    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real expected = h_input(0, icol, ilev) * 0.95;
        Real actual = h_output(icol, ilev);
        Real error = std::abs(expected - actual);

        if (error > 1e-15) {
          correct = false;
          max_error = std::max(max_error, error);
        }
      }
    }

    INFO("qr subview access validation");
    INFO("  Max error: " << max_error);
    REQUIRE(correct);
  }

  SECTION("qr_values_reasonable") {
    // Check that initialized values are in reasonable range for rain
    auto h_qr = Kokkos::create_mirror_view(qr_tracer);
    Kokkos::deep_copy(h_qr, qr_tracer);

    bool values_reasonable = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real val = h_qr(0, icol, ilev);
        // Rain mixing ratio should be non-negative and < 0.1 kg/kg
        if (val < 0.0 || val > 0.1 || std::isnan(val) || std::isinf(val)) {
          values_reasonable = false;
        }
      }
    }
    REQUIRE(values_reasonable);
  }
}

TEST_CASE("precip_2d_surface_tracer_access") {
  const int ncols = 10;
  const int ntracers = 1;  // Only bulk water

  // Create 2D tracer fields for surface precip (tracer, col)
  using view_2d = Kokkos::View<Real**>;
  view_2d precip_liq_tracer("precip_liq_tracer", ntracers, ncols);
  view_2d precip_ice_tracer("precip_ice_tracer", ntracers, ncols);

  // Output views
  using view_1d = Kokkos::View<Real*>;
  view_1d output_liq("output_liq", ncols);
  view_1d output_ice("output_ice", ncols);

  // Initialize with typical accumulated precipitation values
  Kokkos::parallel_for("init_precip", ncols,
    KOKKOS_LAMBDA(int icol) {
      precip_liq_tracer(0, icol) = static_cast<Real>(icol) * 0.5;  // kg/m2
      precip_ice_tracer(0, icol) = static_cast<Real>(icol) * 0.2;  // kg/m2
    });

  Kokkos::fence();

  SECTION("precip_liq_subview_access") {
    // Extract slot 0 using subview (as done in P3 interface)
    auto precip_liq_bulk = Kokkos::subview(precip_liq_tracer, 0, Kokkos::ALL());

    // Run kernel with subview
    Kokkos::parallel_for("precip_liq_kernel", ncols,
      SurfPrecipSubviewKernel<decltype(precip_liq_bulk), view_1d>{precip_liq_bulk, output_liq});

    Kokkos::fence();

    // Verify output
    auto h_input = Kokkos::create_mirror_view(precip_liq_tracer);
    auto h_output = Kokkos::create_mirror_view(output_liq);
    Kokkos::deep_copy(h_input, precip_liq_tracer);
    Kokkos::deep_copy(h_output, output_liq);

    bool correct = true;
    Real max_error = 0.0;

    for (int icol = 0; icol < ncols; ++icol) {
      Real expected = h_input(0, icol) + 0.1;
      Real actual = h_output(icol);
      Real error = std::abs(expected - actual);

      if (error > 1e-15) {
        correct = false;
        max_error = std::max(max_error, error);
      }
    }

    INFO("precip_liq_surf_mass subview access validation");
    INFO("  Max error: " << max_error);
    REQUIRE(correct);
  }

  SECTION("precip_ice_subview_access") {
    // Extract slot 0 using subview (as done in P3 interface)
    auto precip_ice_bulk = Kokkos::subview(precip_ice_tracer, 0, Kokkos::ALL());

    // Run kernel with subview
    Kokkos::parallel_for("precip_ice_kernel", ncols,
      SurfPrecipSubviewKernel<decltype(precip_ice_bulk), view_1d>{precip_ice_bulk, output_ice});

    Kokkos::fence();

    // Verify output
    auto h_input = Kokkos::create_mirror_view(precip_ice_tracer);
    auto h_output = Kokkos::create_mirror_view(output_ice);
    Kokkos::deep_copy(h_input, precip_ice_tracer);
    Kokkos::deep_copy(h_output, output_ice);

    bool correct = true;
    Real max_error = 0.0;

    for (int icol = 0; icol < ncols; ++icol) {
      Real expected = h_input(0, icol) + 0.1;
      Real actual = h_output(icol);
      Real error = std::abs(expected - actual);

      if (error > 1e-15) {
        correct = false;
        max_error = std::max(max_error, error);
      }
    }

    INFO("precip_ice_surf_mass subview access validation");
    INFO("  Max error: " << max_error);
    REQUIRE(correct);
  }

  SECTION("surface_precip_values_reasonable") {
    // Check that initialized values are in reasonable range
    auto h_liq = Kokkos::create_mirror_view(precip_liq_tracer);
    auto h_ice = Kokkos::create_mirror_view(precip_ice_tracer);
    Kokkos::deep_copy(h_liq, precip_liq_tracer);
    Kokkos::deep_copy(h_ice, precip_ice_tracer);

    bool values_reasonable = true;
    for (int icol = 0; icol < ncols; ++icol) {
      Real liq_val = h_liq(0, icol);
      Real ice_val = h_ice(0, icol);

      // Accumulated precip should be non-negative
      if (liq_val < 0.0 || std::isnan(liq_val) || std::isinf(liq_val) ||
          ice_val < 0.0 || std::isnan(ice_val) || std::isinf(ice_val)) {
        values_reasonable = false;
      }
    }
    REQUIRE(values_reasonable);
  }
}

TEST_CASE("multi_tracer_independence") {
  const int ncols = 5;
  const int nlevs = 10;
  const int ntracers = 3;  // Bulk + 2 isotopes

  // Create 3D tracer field
  using view_3d = Kokkos::View<Real***>;
  view_3d qr_multi("qr_multi", ntracers, ncols, nlevs);

  // Initialize with different values per tracer
  auto init_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs});

  Kokkos::parallel_for("init_multi", init_policy,
    KOKKOS_LAMBDA(int icol, int ilev) {
      qr_multi(0, icol, ilev) = 1.0;  // Bulk
      qr_multi(1, icol, ilev) = 2.0;  // Isotope 1
      qr_multi(2, icol, ilev) = 3.0;  // Isotope 2
    });

  Kokkos::fence();

  SECTION("slot_independence") {
    // Verify each slot maintains independent values
    auto h_qr = Kokkos::create_mirror_view(qr_multi);
    Kokkos::deep_copy(h_qr, qr_multi);

    bool slots_independent = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        if (h_qr(0, icol, ilev) != 1.0 ||
            h_qr(1, icol, ilev) != 2.0 ||
            h_qr(2, icol, ilev) != 3.0) {
          slots_independent = false;
        }
      }
    }
    REQUIRE(slots_independent);
  }

  SECTION("bulk_slot_zero_extraction") {
    // Extract bulk water (slot 0) via subview
    auto qr_bulk = Kokkos::subview(qr_multi, 0, Kokkos::ALL(), Kokkos::ALL());

    auto h_bulk = Kokkos::create_mirror_view(qr_bulk);
    Kokkos::deep_copy(h_bulk, qr_bulk);

    bool bulk_correct = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        if (h_bulk(icol, ilev) != 1.0) {
          bulk_correct = false;
        }
      }
    }
    REQUIRE(bulk_correct);
  }
}
