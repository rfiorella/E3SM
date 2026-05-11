#include "catch2/catch.hpp"

#include "physics/water_tracers/water_tracers.hpp"
#include "share/core/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {

/*
 * Water Tracer Subview Accessor Tests
 *
 * Purpose: Unit tests for the WaterTracers::get_bulk_water_subview() accessor
 * that extracts CMP=0 from rank-3 water tracer fields.
 *
 * These tests validate the fundamental infrastructure without requiring
 * full physics integration.
 */

namespace {

using namespace scream;
using KT = KokkosTypes<DefaultDevice>;
using Pack = ekat::Pack<Real, SCREAM_PACK_SIZE>;

constexpr int NCOLS = 4;
constexpr int NLEVS = 72;

TEST_CASE("water_tracer_subview_accessor_basic", "[water_tracers]") {
  const int num_cols = NCOLS;
  const int num_levs = NLEVS;

#ifdef SCREAM_TRACE_WATER
  const int num_cmp = SCREAM_NUM_WATER_TRACERS;
#else
  const int num_cmp = 1; // Always 1 when OFF
#endif

  INFO("Testing with num_cmp = " << num_cmp);

  // Create a rank-3 field (ncols, CMP, nlevs)
  using view_3d = typename KT::template view_3d<Pack>;
  view_3d qv_rank3("qv", num_cols, num_cmp, num_levs);

  // Initialize with unique values: qv(i, c, k) = i*1000 + c*100 + k
  auto qv_h = Kokkos::create_mirror_view(qv_rank3);
  for (int i = 0; i < num_cols; ++i) {
    for (int c = 0; c < num_cmp; ++c) {
      for (int k = 0; k < num_levs; ++k) {
        qv_h(i, c, k) = i * 1000.0 + c * 100.0 + k;
      }
    }
  }
  Kokkos::deep_copy(qv_rank3, qv_h);

  // Extract bulk water subview (CMP=0)
  auto qv_bulk = WaterTracers::get_bulk_water_subview(qv_rank3);

  // Verify dimensions: should be (ncols, nlevs)
  REQUIRE(qv_bulk.extent(0) == num_cols);
  REQUIRE(qv_bulk.extent(1) == num_levs);

  // Verify values: bulk water should match CMP=0
  auto qv_bulk_h = Kokkos::create_mirror_view(qv_bulk);
  Kokkos::deep_copy(qv_bulk_h, qv_bulk);

  for (int i = 0; i < num_cols; ++i) {
    for (int k = 0; k < num_levs; ++k) {
      Real expected = i * 1000.0 + 0 * 100.0 + k; // CMP=0
      Real actual = qv_bulk_h(i, k);
      REQUIRE(actual == expected);
    }
  }

  INFO("Subview accessor correctly extracts CMP=0 slice");
}

#ifdef SCREAM_TRACE_WATER
#if SCREAM_NUM_WATER_TRACERS >= 2

TEST_CASE("water_tracer_subview_accessor_n2_modification", "[water_tracers]") {
  const int num_cols = NCOLS;
  const int num_levs = NLEVS;
  const int num_cmp = SCREAM_NUM_WATER_TRACERS;

  REQUIRE(num_cmp >= 2);

  INFO("Testing modification through subview with num_cmp = " << num_cmp);

  // Create a rank-3 field
  using view_3d = typename KT::template view_3d<Pack>;
  view_3d qv_rank3("qv", num_cols, num_cmp, num_levs);

  // Initialize all CMP slices to zero
  Kokkos::deep_copy(qv_rank3, 0.0);

  // Extract bulk water subview (CMP=0)
  auto qv_bulk = WaterTracers::get_bulk_water_subview(qv_rank3);

  // Modify through bulk water view
  Kokkos::parallel_for(
    "modify_bulk_water",
    Kokkos::RangePolicy<>(0, num_cols * num_levs),
    KOKKOS_LAMBDA(const int idx) {
      const int i = idx / num_levs;
      const int k = idx % num_levs;
      qv_bulk(i, k) = i * 1000.0 + k;
    }
  );
  Kokkos::fence();

  // Verify: CMP=0 should be modified, other CMP slices should remain zero
  auto qv_h = Kokkos::create_mirror_view(qv_rank3);
  Kokkos::deep_copy(qv_h, qv_rank3);

  for (int i = 0; i < num_cols; ++i) {
    for (int k = 0; k < num_levs; ++k) {
      // CMP=0 should have the modified values
      Real expected_cmp0 = i * 1000.0 + k;
      REQUIRE(qv_h(i, 0, k) == expected_cmp0);

      // Other CMP slices should still be zero
      for (int c = 1; c < num_cmp; ++c) {
        REQUIRE(qv_h(i, c, k) == 0.0);
      }
    }
  }

  INFO("Modifications through bulk water subview correctly update only CMP=0");
}

TEST_CASE("water_tracer_all_species_subview", "[water_tracers]") {
  const int num_cols = NCOLS;
  const int num_levs = NLEVS;
  const int num_cmp = SCREAM_NUM_WATER_TRACERS;

  REQUIRE(num_cmp >= 2);

  INFO("Testing subview accessor for all 5 water mass species");

  // Create all 5 mass tracer fields
  using view_3d = typename KT::template view_3d<Pack>;
  view_3d qv_rank3("qv", num_cols, num_cmp, num_levs);
  view_3d qc_rank3("qc", num_cols, num_cmp, num_levs);
  view_3d qi_rank3("qi", num_cols, num_cmp, num_levs);
  view_3d qr_rank3("qr", num_cols, num_cmp, num_levs);
  view_3d qm_rank3("qm", num_cols, num_cmp, num_levs);

  // Initialize each species with unique pattern
  auto init_species = [&](view_3d& field, Real offset) {
    auto field_h = Kokkos::create_mirror_view(field);
    for (int i = 0; i < num_cols; ++i) {
      for (int c = 0; c < num_cmp; ++c) {
        for (int k = 0; k < num_levs; ++k) {
          field_h(i, c, k) = offset + i * 1000.0 + c * 100.0 + k;
        }
      }
    }
    Kokkos::deep_copy(field, field_h);
  };

  init_species(qv_rank3, 10000.0);
  init_species(qc_rank3, 20000.0);
  init_species(qi_rank3, 30000.0);
  init_species(qr_rank3, 40000.0);
  init_species(qm_rank3, 50000.0);

  // Extract bulk water for all species
  auto qv_bulk = WaterTracers::get_bulk_water_subview(qv_rank3);
  auto qc_bulk = WaterTracers::get_bulk_water_subview(qc_rank3);
  auto qi_bulk = WaterTracers::get_bulk_water_subview(qi_rank3);
  auto qr_bulk = WaterTracers::get_bulk_water_subview(qr_rank3);
  auto qm_bulk = WaterTracers::get_bulk_water_subview(qm_rank3);

  // Verify all dimensions
  REQUIRE(qv_bulk.extent(0) == num_cols);
  REQUIRE(qv_bulk.extent(1) == num_levs);
  REQUIRE(qc_bulk.extent(0) == num_cols);
  REQUIRE(qc_bulk.extent(1) == num_levs);
  REQUIRE(qi_bulk.extent(0) == num_cols);
  REQUIRE(qi_bulk.extent(1) == num_levs);
  REQUIRE(qr_bulk.extent(0) == num_cols);
  REQUIRE(qr_bulk.extent(1) == num_levs);
  REQUIRE(qm_bulk.extent(0) == num_cols);
  REQUIRE(qm_bulk.extent(1) == num_levs);

  // Verify values for each species
  auto verify_species = [&](auto bulk_view, Real offset, const char* name) {
    auto bulk_h = Kokkos::create_mirror_view(bulk_view);
    Kokkos::deep_copy(bulk_h, bulk_view);

    for (int i = 0; i < num_cols; ++i) {
      for (int k = 0; k < num_levs; ++k) {
        Real expected = offset + i * 1000.0 + 0 * 100.0 + k; // CMP=0
        Real actual = bulk_h(i, k);
        INFO("Checking " << name << " at (" << i << ", " << k << ")");
        REQUIRE(actual == expected);
      }
    }
  };

  verify_species(qv_bulk, 10000.0, "qv");
  verify_species(qc_bulk, 20000.0, "qc");
  verify_species(qi_bulk, 30000.0, "qi");
  verify_species(qr_bulk, 40000.0, "qr");
  verify_species(qm_bulk, 50000.0, "qm");

  INFO("All 5 water mass species correctly extracted via subview accessor");
}

#endif // SCREAM_NUM_WATER_TRACERS >= 2
#endif // SCREAM_TRACE_WATER

} // anonymous namespace

} // namespace scream
