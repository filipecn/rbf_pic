#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

// tests for a simple tree
//       |10          |11
//  ------------------------
//  :           :           :
// 3_    c0    _4_    c1   _5_
//  :           :           :
//  -----|8-----------|9----
//  :           :           :
// 0_    c2    _1_    c3   _2_
//  :           :           :
//  ------------------------
//       |6           |7

TEST(DomainBoundaryHandler2, SimQuadTree) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  qt.setBoundaryHandler(&bh);
  for (size_t i = 0; i < qt.faceCount(); i++)
    EXPECT_EQ(bh.faceBoundaryType(i), DomainBoundaryHandler2::BoundaryType::NONE);
  EXPECT_EQ(qt.faceCount(), 12u);
  bh.markFaces(&qt);
  EXPECT_EQ(bh.faceBoundaryType(8), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(9), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(10), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(11), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(4), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(5), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(6), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(7), DomainBoundaryHandler2::BoundaryType::NEUMANN);
}

TEST(DomainBoundaryHandler2, SimRegularDomain) {
  SimRegularDomain2 rd(10);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  rd.setBoundaryHandler(&bh);
  for (size_t i = 0; i < rd.faceCount(); i++)
    EXPECT_EQ(bh.faceBoundaryType(i), DomainBoundaryHandler2::BoundaryType::NONE);
  return;
  EXPECT_EQ(rd.faceCount(), 12u);
  bh.markFaces(&rd);
  EXPECT_EQ(bh.faceBoundaryType(8), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(9), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(10), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(11), DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  EXPECT_EQ(bh.faceBoundaryType(4), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(5), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(6), DomainBoundaryHandler2::BoundaryType::NEUMANN);
  EXPECT_EQ(bh.faceBoundaryType(7), DomainBoundaryHandler2::BoundaryType::NEUMANN);
}

TEST(DomainBoundaryHandler2, iterateBoundaryFaces) {
  SimRegularDomain2 rd(2);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  rd.setBoundaryHandler(&bh);
  std::vector<size_t> faces;
  bh.markFaces(&rd);
  bh.iterateBoundaryFaces(
      [&](size_t i, DomainBoundaryHandler2::BoundaryType b) {
        EXPECT_EQ(b, DomainBoundaryHandler2::BoundaryType::DIRICHLET);
        faces.emplace_back(i);
      });
  std::sort(faces.begin(), faces.end());
  std::vector<size_t> expected = {0, 2, 3, 5, 6, 7, 10, 11};
  EXPECT_EQ(expected.size(), faces.size());
  for (size_t i = 0; i < faces.size(); i++)
    EXPECT_EQ(faces[i], expected[i]);
}

TEST(DomainBoundaryHandler2, iterateBoundaryCells) {
  // 5 AAAAAA  30 31 32 33 34 35
  // 4 ASSSAA  24 25 26 27 28 29
  // 3 AAFAAA  18 19 20 21 22 23
  // 2 AFFFAA  12 13 14 15 16 17
  // 1 AAFAAS  06 07 08 09 10 11
  // 0 AAAAAA  00 01 02 03 04 05
  //   012345
  SimRegularDomain2 rd(6);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  rd.setBoundaryHandler(&bh);
  for (size_t i = 0; i < rd.cellCount(); i++)
    bh.setCellType(i, DomainBoundaryHandler2::MaterialType::AIR);
  bh.setCellType(1 * 6 + 5, DomainBoundaryHandler2::MaterialType::SOLID);
  bh.setCellType(4 * 6 + 1, DomainBoundaryHandler2::MaterialType::SOLID);
  bh.setCellType(4 * 6 + 2, DomainBoundaryHandler2::MaterialType::SOLID);
  bh.setCellType(4 * 6 + 3, DomainBoundaryHandler2::MaterialType::SOLID);
  bh.setCellType(1 * 6 + 2, DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(2 * 6 + 1, DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(2 * 6 + 2, DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(2 * 6 + 3, DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(3 * 6 + 2, DomainBoundaryHandler2::MaterialType::FLUID);
  {
    std::vector<size_t> cells;
    bh.iterateBoundaryCells(&rd, DomainBoundaryHandler2::MaterialType::SOLID,
                            [&](size_t i) {
                              auto ns = rd.cellNeighborhood(i);
                              cells.emplace_back(i);
                            });
    std::sort(cells.begin(), cells.end());
    std::vector<size_t> expected = {11, 25, 26, 27};
    EXPECT_EQ(cells.size(), expected.size());
    for (size_t i = 0; i < cells.size(); i++)
      EXPECT_EQ(cells[i], expected[i]);
  }
  {
    std::vector<size_t> cells;
    bh.iterateBoundaryCells(&rd, DomainBoundaryHandler2::MaterialType::FLUID,
                            [&](size_t i) {
                              auto ns = rd.cellNeighborhood(i);
                              cells.emplace_back(i);
                            });
    std::sort(cells.begin(), cells.end());
    std::vector<size_t> expected = {8, 13, 15, 20};
    EXPECT_EQ(cells.size(), expected.size());
    for (size_t i = 0; i < cells.size(); i++)
      EXPECT_EQ(cells[i], expected[i]);
  }
  {
    std::vector<size_t> cells;
    bh.iterateBoundaryCells(&rd, DomainBoundaryHandler2::MaterialType::AIR,
                            [&](size_t i) {
                              auto ns = rd.cellNeighborhood(i);
                              cells.emplace_back(i);
                            });
    std::sort(cells.begin(), cells.end());
    std::vector<size_t> expected = {0,  1,  2,  3,  4,  5,  6,  7,  9,
                                    10, 12, 16, 17, 18, 19, 21, 23, 24,
                                    28, 29, 30, 31, 32, 33, 34, 35};
    EXPECT_EQ(cells.size(), expected.size());
    for (size_t i = 0; i < cells.size(); i++)
      EXPECT_EQ(cells[i], expected[i]);
  }
}
