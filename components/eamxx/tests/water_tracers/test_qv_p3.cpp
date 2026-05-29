#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_tracer_access.hpp"
#include "share/grid/point_grid.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {

TEST_CASE("p3_qv_tracer_access", "[water_tracers][p3]") {
  using namespace ShortFieldTagsNames;
  using Pack = ekat::Pack<Real, 1>;

  // Create a simple grid
  constexpr int ncols = 4;
  constexpr int nlevs = 72;
  #ifndef SCREAM_NUM_TRACERS
  constexpr int ntracers = 1;
  #else
  constexpr int ntracers = SCREAM_NUM_TRACERS;
  #endif

  auto grid = create_point_grid("test_grid", ncols, nlevs, ncols);

  // Create tracer-aware qv field
  FieldIdentifier qv_fid("qv", grid->get_3d_tracer_layout(ntracers), ekat::units::kg/ekat::units::kg, grid->name());
  Field qv_field(qv_fid);
  qv_field.allocate_view();

  // Get 3D tracer view and 2D bulk subview
  auto qv_3d = qv_field.get_view<Pack***>();
  auto qv_bulk = get_tracer_bulk_subview(qv_3d);

  // Initialize slot-0 with test data
  const Real test_value = 0.012; // kg/kg, typical atmospheric qv
  Kokkos::deep_copy(qv_3d, 0.0);
  Kokkos::parallel_for("init_qv_bulk", ncols * nlevs, KOKKOS_LAMBDA(int idx) {
    int icol = idx / nlevs;
    int ilev = idx % nlevs;
    qv_bulk(icol, ilev)[0] = test_value;
  });
  Kokkos::fence();

  // Verify P3 can correctly read qv slot-0
  bool all_correct = true;
  Kokkos::parallel_reduce("verify_qv_access", ncols * nlevs,
    KOKKOS_LAMBDA(int idx, bool& lcorrect) {
      int icol = idx / nlevs;
      int ilev = idx % nlevs;
      const Real val = qv_bulk(icol, ilev)[0];
      if (std::abs(val - test_value) > 1e-15) {
        lcorrect = false;
      }
    }, Kokkos::LAnd<bool>(all_correct));
  Kokkos::fence();

  REQUIRE(all_correct);

  // Test that slot-0 subview is independent from other tracer slots
  if (ntracers > 1) {
    const Real tracer1_value = 0.008;
    Kokkos::parallel_for("init_qv_tracer1", ncols * nlevs, KOKKOS_LAMBDA(int idx) {
      int icol = idx / nlevs;
      int ilev = idx % nlevs;
      qv_3d(1, icol, ilev)[0] = tracer1_value;
    });
    Kokkos::fence();

    // Verify bulk (slot-0) unchanged
    bool bulk_unchanged = true;
    Kokkos::parallel_reduce("verify_bulk_unchanged", ncols * nlevs,
      KOKKOS_LAMBDA(int idx, bool& lcorrect) {
        int icol = idx / nlevs;
        int ilev = idx % nlevs;
        const Real val = qv_bulk(icol, ilev)[0];
        if (std::abs(val - test_value) > 1e-15) {
          lcorrect = false;
        }
      }, Kokkos::LAnd<bool>(bulk_unchanged));
    Kokkos::fence();

    REQUIRE(bulk_unchanged);
  }
}

} // namespace scream
