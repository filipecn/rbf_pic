#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(FastSweep2, propagate) {
  {
    SimRegularDomain2 rd(10);
    auto mask = rd.addCellCenteredByteField(0);
    rd.byteAtCell(mask, 0) = 1;
    auto field = rd.addCellCenteredScalarField(INFINITY);
    rd.scalarAtCell(field, 0) = 0.;
    FastSweep2::propagate(&rd, field, mask);
    FastSweep2::propagate(&rd, field, mask);
    FastSweep2::propagate(&rd, field, mask);
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++)
        ASSERT_NEAR(rd.cellCenterPosition(0).distance(
                        rd.cellCenterPosition(i * 10 + j)),
                    rd.scalarAtCell(field, i * 10 + j) / 10., 1e-1);
  }
}

TEST(FastMarchingMethod2, propagate) {
  {
    SimRegularDomain2 rd(10);
    auto frozen = rd.addCellCenteredByteField(0);
    rd.byteAtCell(frozen, 0) = 1;
    auto field = rd.addCellCenteredScalarField(INFINITY);
    rd.scalarAtCell(field, 0) = 0.;
    FastMarchingMethod2::propagate(&rd, field, frozen, {0});
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++)
        ASSERT_NEAR(rd.cellCenterPosition(0).distance(
                        rd.cellCenterPosition(i * 10 + j)),
                    rd.scalarAtCell(field, i * 10 + j) / 10., 1e-1);
  }
}
