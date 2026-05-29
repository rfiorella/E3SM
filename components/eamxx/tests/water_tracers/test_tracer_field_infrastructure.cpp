/*
 * Test tracer field infrastructure
 *
 * Validates:
 * - TRACER FieldTag is recognized
 * - get_3d_tracer_layout() produces correct layout
 * - Layout dimensions match expected (tracer, col, lev)
 */

#include "catch2/catch.hpp"

#include "share/field/field_tag.hpp"
#include "share/field/field_layout.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"

#include "ekat_comm.hpp"

using namespace scream;
using namespace ShortFieldTagsNames;

TEST_CASE("tracer_fieldtag_exists") {
  // Verify TRACER is a valid FieldTag
  auto tracer_tag = FieldTag::Tracer;
  REQUIRE(e2str(tracer_tag) == "tracer");
}

TEST_CASE("point_grid_3d_tracer_layout") {
  auto comm = ekat::Comm::self();

  const int ncols = 10;
  const int nlevs = 72;
  const int ntracers = 3;

  auto grid = std::make_shared<PointGrid>("test_grid", ncols, nlevs, comm);

  SECTION("tracer_layout_bulk_only") {
    // Single tracer (bulk water only)
    auto layout = grid->get_3d_tracer_layout(1);

    REQUIRE(layout.rank() == 3);
    REQUIRE(layout.has_tag(TRACER));
    REQUIRE(layout.has_tag(COL));
    REQUIRE(layout.has_tag(LEV));

    // Check dimension order: (tracer, col, lev)
    REQUIRE(layout.tag(0) == TRACER);
    REQUIRE(layout.tag(1) == COL);
    REQUIRE(layout.tag(2) == LEV);

    // Check dimensions
    REQUIRE(layout.dim(TRACER) == 1);
    REQUIRE(layout.dim(COL) == ncols);
    REQUIRE(layout.dim(LEV) == nlevs);
  }

  SECTION("tracer_layout_multi_tracer") {
    // Multiple tracers (bulk + isotopes)
    auto layout = grid->get_3d_tracer_layout(ntracers);

    REQUIRE(layout.rank() == 3);
    REQUIRE(layout.dim(TRACER) == ntracers);
    REQUIRE(layout.dim(COL) == ncols);
    REQUIRE(layout.dim(LEV) == nlevs);
  }

  SECTION("tracer_layout_custom_name") {
    // Test custom dimension name
    auto layout = grid->get_3d_tracer_layout(ntracers, "water_isotope");

    REQUIRE(layout.name(0) == "water_isotope");
  }
}

TEST_CASE("se_grid_3d_tracer_layout") {
  auto comm = ekat::Comm::self();

  const int nelem = 5;
  const int ngp = 4;  // 4x4 GLL points per element
  const int nlevs = 72;
  const int ntracers = 3;

  auto grid = std::make_shared<SEGrid>("test_se_grid", nelem, ngp, nlevs, comm);

  SECTION("tracer_layout_se") {
    // SE grid tracer layout: (tracer, elem, gp, gp, lev)
    auto layout = grid->get_3d_tracer_layout(ntracers);

    REQUIRE(layout.rank() == 5);
    REQUIRE(layout.has_tag(TRACER));
    REQUIRE(layout.has_tag(EL));
    REQUIRE(layout.has_tag(GP));
    REQUIRE(layout.has_tag(LEV));

    // Check dimension order
    REQUIRE(layout.tag(0) == TRACER);
    REQUIRE(layout.tag(1) == EL);
    REQUIRE(layout.tag(2) == GP);
    REQUIRE(layout.tag(3) == GP);
    REQUIRE(layout.tag(4) == LEV);

    // Check dimensions
    REQUIRE(layout.dim(TRACER) == ntracers);
    REQUIRE(layout.dims()[1] == nelem);  // Element dim
    REQUIRE(layout.dims()[2] == ngp);    // GP dim 1
    REQUIRE(layout.dims()[3] == ngp);    // GP dim 2
    REQUIRE(layout.dims()[4] == nlevs);  // Level dim
  }
}

TEST_CASE("tracer_layout_consistency") {
  // Ensure tracer layouts are consistent across grid types
  auto comm = ekat::Comm::self();

  const int ncols = 10;
  const int nlevs = 72;
  const int ntracers = 2;

  auto pg = std::make_shared<PointGrid>("point", ncols, nlevs, comm);
  auto pg_layout = pg->get_3d_tracer_layout(ntracers);

  // All tracer layouts should have TRACER as first dimension
  REQUIRE(pg_layout.tag(0) == TRACER);
  REQUIRE(pg_layout.dim(0) == ntracers);
}
