/*
 * Test tracer access patterns for BFB validation
 *
 * Compares tracer-dimension field access patterns against scalar baseline:
 * - EXPLICIT indexing: qv(0, icol, ilev)
 * - SUBVIEW accessor: auto qv_bulk = subview(qv, 0, ALL, ALL); qv_bulk(icol, ilev)
 *
 * Both patterns must produce BFB-identical results vs scalar (col, lev) layout.
 */

#include "catch2/catch.hpp"

#include "share/field/field_tracer_access.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "ekat_comm.hpp"
#include "ekat_pack.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

using namespace scream;

// Simple kernel that performs standard atmospheric calculation
// This mimics a typical physics operation on qv
template<typename ViewT>
struct AtmosphericKernel {
  ViewT qv;
  ViewT temperature;
  ViewT output;

  KOKKOS_INLINE_FUNCTION
  void operator()(int icol, int ilev) const {
    // Simple calculation: output = qv * exp(-temperature / 300)
    // This is representative of typical atmospheric physics operations
    output(icol, ilev) = qv(icol, ilev) * Kokkos::exp(-temperature(icol, ilev) / 300.0);
  }
};

// Kernel using explicit indexing into tracer dimension
template<typename ViewT>
struct TracerExplicitKernel {
  ViewT qv_tracer;      // (tracer, col, lev)
  ViewT temperature;    // (col, lev)
  ViewT output;         // (col, lev)

  KOKKOS_INLINE_FUNCTION
  void operator()(int icol, int ilev) const {
    // Access slot 0 explicitly
    output(icol, ilev) = qv_tracer(0, icol, ilev) *
                         Kokkos::exp(-temperature(icol, ilev) / 300.0);
  }
};

// Kernel using subview accessor
template<typename ViewT, typename SubviewT>
struct TracerSubviewKernel {
  SubviewT qv_bulk;     // 2D subview of slot 0
  ViewT temperature;    // (col, lev)
  ViewT output;         // (col, lev)

  KOKKOS_INLINE_FUNCTION
  void operator()(int icol, int ilev) const {
    output(icol, ilev) = qv_bulk(icol, ilev) *
                         Kokkos::exp(-temperature(icol, ilev) / 300.0);
  }
};

TEST_CASE("tracer_access_bfb_validation") {
  auto comm = ekat::Comm::self();

  const int ncols = 10;
  const int nlevs = 72;
  const int ntracers = 1;  // Only bulk water for BFB test

  auto grid = std::make_shared<PointGrid>("test_grid", ncols, nlevs, comm);

  // Create scalar baseline views (col, lev)
  using view_2d = Kokkos::View<Real**, DefaultDevice>;
  view_2d qv_scalar("qv_scalar", ncols, nlevs);
  view_2d temp_scalar("temp_scalar", ncols, nlevs);
  view_2d output_scalar("output_scalar", ncols, nlevs);

  // Create tracer-dimension views (tracer, col, lev)
  using view_3d = Kokkos::View<Real***, DefaultDevice>;
  view_3d qv_tracer("qv_tracer", ntracers, ncols, nlevs);
  view_2d output_explicit("output_explicit", ncols, nlevs);
  view_2d output_subview("output_subview", ncols, nlevs);

  // Initialize with representative atmospheric values
  auto init_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs});

  Kokkos::parallel_for("init_scalar", init_policy,
    KOKKOS_LAMBDA(int icol, int ilev) {
      Real height_frac = static_cast<Real>(ilev) / nlevs;
      Real col_offset = static_cast<Real>(icol) * 0.01;

      // Typical atmospheric profile
      qv_scalar(icol, ilev) = 0.01 * Kokkos::exp(-height_frac * 5.0) + col_offset;
      temp_scalar(icol, ilev) = 288.15 - height_frac * 60.0 + col_offset * 10.0;
    });

  // Copy scalar qv to tracer array slot 0
  Kokkos::parallel_for("init_tracer", init_policy,
    KOKKOS_LAMBDA(int icol, int ilev) {
      qv_tracer(0, icol, ilev) = qv_scalar(icol, ilev);
    });

  Kokkos::fence();

  SECTION("scalar_baseline") {
    // Run baseline calculation with scalar layout
    Kokkos::parallel_for("scalar_kernel", init_policy,
      AtmosphericKernel<view_2d>{qv_scalar, temp_scalar, output_scalar});
    Kokkos::fence();

    // Verify output is reasonable
    auto h_output = Kokkos::create_mirror_view(output_scalar);
    Kokkos::deep_copy(h_output, output_scalar);

    bool values_reasonable = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real val = h_output(icol, ilev);
        if (val < 0.0 || val > 0.1 || std::isnan(val) || std::isinf(val)) {
          values_reasonable = false;
        }
      }
    }
    REQUIRE(values_reasonable);
  }

  SECTION("explicit_indexing_bfb") {
    // Run baseline
    Kokkos::parallel_for("scalar_kernel", init_policy,
      AtmosphericKernel<view_2d>{qv_scalar, temp_scalar, output_scalar});

    // Run with explicit indexing
    Kokkos::parallel_for("explicit_kernel", init_policy,
      TracerExplicitKernel<view_2d>{qv_tracer, temp_scalar, output_explicit});

    Kokkos::fence();

    // Compare BFB
    auto h_scalar = Kokkos::create_mirror_view(output_scalar);
    auto h_explicit = Kokkos::create_mirror_view(output_explicit);
    Kokkos::deep_copy(h_scalar, output_scalar);
    Kokkos::deep_copy(h_explicit, output_explicit);

    bool is_bfb = true;
    int diff_count = 0;
    Real max_diff = 0.0;

    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real scalar_val = h_scalar(icol, ilev);
        Real explicit_val = h_explicit(icol, ilev);

        if (scalar_val != explicit_val) {
          is_bfb = false;
          diff_count++;
          Real diff = std::abs(scalar_val - explicit_val);
          max_diff = std::max(max_diff, diff);
        }
      }
    }

    INFO("Explicit indexing BFB check:");
    INFO("  Differences found: " << diff_count);
    INFO("  Max difference: " << max_diff);

    REQUIRE(is_bfb);
  }

  SECTION("subview_bfb") {
    // Run baseline
    Kokkos::parallel_for("scalar_kernel", init_policy,
      AtmosphericKernel<view_2d>{qv_scalar, temp_scalar, output_scalar});

    // Run with subview accessor
    auto qv_bulk_subview = Kokkos::subview(qv_tracer, 0, Kokkos::ALL(), Kokkos::ALL());
    Kokkos::parallel_for("subview_kernel", init_policy,
      TracerSubviewKernel<view_2d, decltype(qv_bulk_subview)>{
        qv_bulk_subview, temp_scalar, output_subview});

    Kokkos::fence();

    // Compare BFB
    auto h_scalar = Kokkos::create_mirror_view(output_scalar);
    auto h_subview = Kokkos::create_mirror_view(output_subview);
    Kokkos::deep_copy(h_scalar, output_scalar);
    Kokkos::deep_copy(h_subview, output_subview);

    bool is_bfb = true;
    int diff_count = 0;
    Real max_diff = 0.0;

    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real scalar_val = h_scalar(icol, ilev);
        Real subview_val = h_subview(icol, ilev);

        if (scalar_val != subview_val) {
          is_bfb = false;
          diff_count++;
          Real diff = std::abs(scalar_val - subview_val);
          max_diff = std::max(max_diff, diff);
        }
      }
    }

    INFO("Subview BFB check:");
    INFO("  Differences found: " << diff_count);
    INFO("  Max difference: " << max_diff);

    REQUIRE(is_bfb);
  }
}

TEST_CASE("tracer_accessor_helpers") {
  const int ntracers = 1;
  const int ncols = 5;
  const int nlevs = 10;

  using view_3d = Kokkos::View<Real***, DefaultDevice>;
  view_3d qv_tracer("qv_tracer", ntracers, ncols, nlevs);

  // Initialize
  Kokkos::parallel_for("init", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(int icol, int ilev) {
      qv_tracer(0, icol, ilev) = static_cast<Real>(icol * 100 + ilev);
    });
  Kokkos::fence();

  SECTION("explicit_accessor_helper") {
    using namespace scream;

    view_3d output("output", ntracers, ncols, nlevs);

    Kokkos::parallel_for("test_explicit",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
      KOKKOS_LAMBDA(int icol, int ilev) {
        auto accessor = get_tracer_bulk_explicit(qv_tracer);
        output(0, icol, ilev) = accessor(icol, ilev) * 2.0;
      });
    Kokkos::fence();

    auto h_input = Kokkos::create_mirror_view(qv_tracer);
    auto h_output = Kokkos::create_mirror_view(output);
    Kokkos::deep_copy(h_input, qv_tracer);
    Kokkos::deep_copy(h_output, output);

    bool correct = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real expected = h_input(0, icol, ilev) * 2.0;
        Real actual = h_output(0, icol, ilev);
        if (expected != actual) {
          correct = false;
        }
      }
    }
    REQUIRE(correct);
  }

  SECTION("subview_accessor_helper") {
    using namespace scream;

    using view_2d = Kokkos::View<Real**, DefaultDevice>;
    view_2d output("output", ncols, nlevs);

    Kokkos::parallel_for("test_subview",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncols, nlevs}),
      KOKKOS_LAMBDA(int icol, int ilev) {
        auto bulk_view = get_tracer_bulk_subview(qv_tracer);
        output(icol, ilev) = bulk_view(icol, ilev) * 2.0;
      });
    Kokkos::fence();

    auto h_input = Kokkos::create_mirror_view(qv_tracer);
    auto h_output = Kokkos::create_mirror_view(output);
    Kokkos::deep_copy(h_input, qv_tracer);
    Kokkos::deep_copy(h_output, output);

    bool correct = true;
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        Real expected = h_input(0, icol, ilev) * 2.0;
        Real actual = h_output(icol, ilev);
        if (expected != actual) {
          correct = false;
        }
      }
    }
    REQUIRE(correct);
  }
}
