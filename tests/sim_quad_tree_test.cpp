#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(SimQuadTree, Contructors) {
  { SimQuadTree qt; }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
    EXPECT_EQ(qt.cellCount(), 4u * 4u * 4u * 4u);
  }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
    EXPECT_EQ(qt.faceCount(), 12u);
  }
}

TEST(SimQuadTree, RandomAccessors) {
  SimQuadTree gd(BBox2D::make_unit_bbox(), 3);
  {
    auto field = gd.addCellCenteredScalarField(10.);
    for (unsigned i = 0; i < gd.cellCount(); i++) {
      ASSERT_NEAR(gd.scalarAtCell(field, i), 10., 1e-8);
      gd.scalarAtCell(field, i) = i * 1.;
    }
    for (unsigned i = 0; i < gd.cellCount(); i++)
      ASSERT_NEAR(gd.scalarAtCell(field, i), i * 1., 1e-8);
  }
  {
    auto field = gd.addFaceCenteredScalarField(10.);
    for (unsigned i = 0; i < gd.faceCount(); i++) {
      ASSERT_NEAR(gd.scalarAtFace(field, i), 10., 1e-8);
      gd.scalarAtFace(field, i) = i * 1.;
    }
    for (unsigned i = 0; i < gd.faceCount(); i++)
      ASSERT_NEAR(gd.scalarAtFace(field, i), i * 1., 1e-8);
  }
  {
    auto field = gd.addCellCenteredByteField(255u);
    for (unsigned i = 0; i < gd.cellCount(); i++) {
      EXPECT_EQ(gd.byteAtCell(field, i), 255u);
      gd.byteAtCell(field, i) = i;
    }
    for (unsigned i = 0; i < gd.cellCount(); i++)
      EXPECT_EQ(gd.byteAtCell(field, i), i);
  }
}

TEST(SimQuadTree, CellRegion) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  auto r = qt.cellRegion(2);
  EXPECT_EQ(Point2d(0., 0.), r.lower());
  EXPECT_EQ(Point2d(0.5, 0.5), r.upper());
  r = qt.cellRegion(1);
  EXPECT_EQ(Point2d(0.5, 0.5), r.lower());
  EXPECT_EQ(Point2d(1., 1.), r.upper());
}

TEST(SimQuadTree, IterateCell) {
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
    size_t id = qt.addCellCenteredScalarField(1.0);
    size_t k = 0;
    qt.iterateCellScalar(id, [&k](double &v, Point2d p) {
      EXPECT_EQ(IS_EQUAL(v, 1.0), true);
      UNUSED_VARIABLE(p);
      k++;
    });
    EXPECT_EQ(k, qt.cellCount());
  }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
    size_t id = qt.addCellCenteredScalarField(1.0);
    size_t k = 0;
    qt.iterateCellScalar(id, [&k](double &v, Point2d p) {
      static int x = 1;
      static int y = 2;
      EXPECT_EQ(IS_EQUAL(v, 1.0), true);
      EXPECT_EQ(p == Point2d(x * 0.5 - 0.25, y * 0.5 - 0.25), true);
      x++;
      if (x > 2) {
        x = 1;
        y--;
      }
      v = 1.0 * k++;
    });
    EXPECT_EQ(k, qt.cellCount());
    EXPECT_EQ(k, 4u);
    k = 0;
    qt.iterateCellScalar(id, [&k](double &v, Point2d p) {
      EXPECT_EQ(IS_EQUAL(v, 1.0 * k++), true);
      UNUSED_VARIABLE(p);
    });
  }
}

// tests for a simple tree
//       |6           |7
//  ------------------------
//  :           :           :
// 10_   c0    _0_    c1   _8_
//  :           :           :
//  -----|2------------|3----
//  :           :           :
// 11_   c2    _1_    c3   _9_
//  :           :           :
//  ------------------------
//       |4           |5
TEST(SimQuadTree, CellNeighbors) {
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
    {
      auto neighbors = qt.cellNeighborhood(0);
      EXPECT_EQ(neighbors.size(), 4u);
      for (size_t i = 0; i < neighbors.size(); i++)
        switch (neighbors[i].id) {
        case -1:
          EXPECT_EQ(neighbors[i].orientation ==
                            SimulationDomain2::Orientation::LEFT ||
                        neighbors[i].orientation ==
                            SimulationDomain2::Orientation::TOP,
                    true);
          break;
        case 1:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::RIGHT,
                    true);
          break;
        case 2:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::BOTTOM,
                    true);
          break;
        default:
          EXPECT_TRUE(false);
        }
    }
    {
      auto neighbors = qt.cellNeighborhood(1);
      EXPECT_EQ(neighbors.size(), 4u);
      for (size_t i = 0; i < neighbors.size(); i++)
        switch (neighbors[i].id) {
        case -1:
          EXPECT_EQ(neighbors[i].orientation ==
                            SimulationDomain2::Orientation::RIGHT ||
                        neighbors[i].orientation ==
                            SimulationDomain2::Orientation::TOP,
                    true);
          break;
        case 0:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::LEFT,
                    true);
          break;
        case 3:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::BOTTOM,
                    true);
          break;
        default:
          EXPECT_TRUE(false);
        }
    }
    {
      auto neighbors = qt.cellNeighborhood(2);
      EXPECT_EQ(neighbors.size(), 4u);
      for (size_t i = 0; i < neighbors.size(); i++)
        switch (neighbors[i].id) {
        case -1:
          EXPECT_EQ(neighbors[i].orientation ==
                            SimulationDomain2::Orientation::LEFT ||
                        neighbors[i].orientation ==
                            SimulationDomain2::Orientation::BOTTOM,
                    true);
          break;
        case 3:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::RIGHT,
                    true);
          break;
        case 0:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::TOP,
                    true);
          break;
        default:
          EXPECT_TRUE(false);
        }
    }
    {
      auto neighbors = qt.cellNeighborhood(3);
      EXPECT_EQ(neighbors.size(), 4u);
      for (size_t i = 0; i < neighbors.size(); i++)
        switch (neighbors[i].id) {
        case -1:
          EXPECT_EQ(neighbors[i].orientation ==
                            SimulationDomain2::Orientation::RIGHT ||
                        neighbors[i].orientation ==
                            SimulationDomain2::Orientation::BOTTOM,
                    true);
          break;
        case 2:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::LEFT,
                    true);
          break;
        case 1:
          EXPECT_EQ(neighbors[i].orientation ==
                        SimulationDomain2::Orientation::TOP,
                    true);
          break;
        default:
          EXPECT_TRUE(false);
        }
    }
  }
}

TEST(SimQuadTree, FaceNeighbors) {
  int expected[12][2] = {{0, 1},  {2, 3},  {0, 2},  {1, 3},  {2, -1}, {3, -1},
                         {0, -1}, {1, -1}, {1, -1}, {3, -1}, {0, -1}, {2, -1}};
  SimulationDomain2::Orientation expectedOrientation[12][2] = {
      {SimulationDomain2::Orientation::LEFT, // 0
       SimulationDomain2::Orientation::RIGHT},
      {SimulationDomain2::Orientation::LEFT, // 1
       SimulationDomain2::Orientation::RIGHT},
      {SimulationDomain2::Orientation::TOP, // 2
       SimulationDomain2::Orientation::BOTTOM},
      {SimulationDomain2::Orientation::TOP, // 3
       SimulationDomain2::Orientation::BOTTOM},
      {SimulationDomain2::Orientation::TOP, // 4
       SimulationDomain2::Orientation::BOTTOM},
      {SimulationDomain2::Orientation::TOP, // 5
       SimulationDomain2::Orientation::BOTTOM},
      {SimulationDomain2::Orientation::BOTTOM, // 6
       SimulationDomain2::Orientation::TOP},
      {SimulationDomain2::Orientation::BOTTOM, // 7
       SimulationDomain2::Orientation::TOP},
      {SimulationDomain2::Orientation::LEFT, // 8
       SimulationDomain2::Orientation::RIGHT},
      {SimulationDomain2::Orientation::LEFT, // 9
       SimulationDomain2::Orientation::RIGHT},
      {SimulationDomain2::Orientation::RIGHT, // 10
       SimulationDomain2::Orientation::LEFT},
      {SimulationDomain2::Orientation::RIGHT, // 11
       SimulationDomain2::Orientation::LEFT}};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  for (size_t i = 0; i < qt.faceCount(); i++) {
    auto neighbors = qt.faceNeighborhood(i);
    for (size_t j = 0; j < 2; j++) {
      EXPECT_EQ(expected[i][j], neighbors[j].id);
      EXPECT_EQ(expectedOrientation[i][j], neighbors[j].orientation);
    }
  }
}

TEST(SimQuadTree, CellCenterStencil) {
  SimulationDomain2::Orientation expectedOrientation[4][5] = {
      {SimulationDomain2::Orientation::OTHER,
       SimulationDomain2::Orientation::RIGHT,
       SimulationDomain2::Orientation::BOTTOM,
       SimulationDomain2::Orientation::TOP,
       SimulationDomain2::Orientation::LEFT},
      {SimulationDomain2::Orientation::OTHER,
       SimulationDomain2::Orientation::LEFT,
       SimulationDomain2::Orientation::BOTTOM,
       SimulationDomain2::Orientation::TOP,
       SimulationDomain2::Orientation::RIGHT},
      {SimulationDomain2::Orientation::OTHER,
       SimulationDomain2::Orientation::RIGHT,
       SimulationDomain2::Orientation::TOP,
       SimulationDomain2::Orientation::BOTTOM,
       SimulationDomain2::Orientation::LEFT},
      {SimulationDomain2::Orientation::OTHER,
       SimulationDomain2::Orientation::LEFT,
       SimulationDomain2::Orientation::TOP,
       SimulationDomain2::Orientation::BOTTOM,
       SimulationDomain2::Orientation::RIGHT}};
  double expectedWeights[5] = {-37.771, 6.29517, 6.29517, 12.5903, 12.5903};
  int expectedIds[4][5] = {{0, 1, 2, -1, -1},
                           {1, 0, 3, -1, -1},
                           {2, 3, 0, -1, -1},
                           {3, 2, 1, -1, -1}};
  Point2d expectedCenters[4][5] = {
      {Point2d(0.25, 0.75), Point2d(0.75, 0.75), Point2d(0.25, 0.25),
       Point2d(0.25, 1), Point2d(0, 0.75)},
      {Point2d(0.75, 0.75), Point2d(0.25, 0.75), Point2d(0.75, 0.25),
       Point2d(0.75, 1), Point2d(1, 0.75)},
      {Point2d(0.25, 0.25), Point2d(0.75, 0.25), Point2d(0.25, 0.75),
       Point2d(0.25, 0), Point2d(0, 0.25)},
      {Point2d(0.75, 0.25), Point2d(0.25, 0.25), Point2d(0.75, 0.75),
       Point2d(0.75, 0), Point2d(1, 0.25)}};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  DomainBoundaryHandler2 bh =
      DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  for (size_t i = 0; i < qt.cellCount(); i++) {
    double wsum = 0.;
    int k = 0;
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          EXPECT_EQ(expectedIds[i][k], id);
          EXPECT_EQ(expectedOrientation[i][k], o);
          EXPECT_EQ(IS_EQUAL_ERROR(w, expectedWeights[k], 1e-4), true);
          EXPECT_EQ(p == expectedCenters[i][k], true);
          wsum += w;
          k++;
        });
    ASSERT_NEAR(wsum, 0., 1e-8);
  }
}

TEST(SimQuadTree, cellCenterPosition) {
  Point2d expectedPosition[4] = {Point2d(0.25, 0.75), Point2d(0.75, 0.75),
                                 Point2d(0.25, 0.25), Point2d(0.75, 0.25)};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  for (size_t i = 0; i < qt.cellCount(); i++)
    EXPECT_EQ(qt.cellCenterPosition(i) == expectedPosition[i], true);
}

TEST(SimQuadTree, faceCenterPosition) {
  Point2d expectedPosition[12] = {
      Point2d(0.5, 0.75), Point2d(0.5, 0.25), Point2d(0.25, 0.5),
      Point2d(0.75, 0.5), Point2d(0.25, 0.),  Point2d(0.75, 0.),
      Point2d(0.25, 1.),  Point2d(0.75, 1.),  Point2d(1., 0.75),
      Point2d(1., 0.25),  Point2d(0., 0.75),  Point2d(0., 0.25)};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  for (size_t i = 0; i < qt.faceCount(); i++)
    EXPECT_EQ(qt.faceCenterPosition(i) == expectedPosition[i], true);
}

TEST(SimQuadTree, IterateDyScalar) {
  Point2d expectedPosition[6] = {Point2d(0.5, 0.75), Point2d(0.5, 0.25),
                                 Point2d(1., 0.75),  Point2d(1., 0.25),
                                 Point2d(0., 0.75),  Point2d(0., 0.25)};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  size_t id = qt.addFaceCenteredScalarField(2.);
  size_t k = 0;
  qt.iterateDyScalar(id, [&](double &v, Point2d p) {
    EXPECT_EQ(IS_EQUAL(v, 2.), true);
    EXPECT_EQ(p == expectedPosition[k++], true);
  });
  EXPECT_EQ(k, 6u);
}

TEST(SimQuadTree, IterateDxScalar) {
  Point2d expectedPosition[6] = {Point2d(0.25, 0.5), Point2d(0.75, 0.5),
                                 Point2d(0.25, 0.),  Point2d(0.75, 0.),
                                 Point2d(0.25, 1.),  Point2d(0.75, 1.)};
  SimQuadTree qt(BBox2D::make_unit_bbox(), 1);
  size_t id = qt.addFaceCenteredScalarField(2.);
  size_t k = 0;
  qt.iterateDxScalar(id, [&](double &v, Point2d p) {
    EXPECT_EQ(IS_EQUAL(v, 2.), true);
    EXPECT_EQ(p == expectedPosition[k++], true);
  });
  EXPECT_EQ(k, 6u);
}

TEST(SimQuadTree, IterateCells) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
  std::cout << qt.cellCount() << std::endl;
  qt.setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                           new PointZGrid2(16));
  BoxSampler sampler;
  for (size_t i = 0; i < 10000; i++) {
    Point2d center = sampler.sample(BBox2D::make_unit_bbox());
    BBox2D b(center - Vector2d(0.3), center + Vector2d(0.3));
    std::vector<size_t> r;
    qt.iterateCells(center, 0.3, [&](size_t id) { r.emplace_back(id); });
    // brute force
    std::vector<size_t> ids;
    for (size_t j = 0; j < qt.cellCount(); j++)
      if (b.inside(qt.cellCenterPosition(j)))
        ids.emplace_back(j);
    EXPECT_EQ(r.size(), ids.size());
    std::sort(r.begin(), r.end());
    std::sort(ids.begin(), ids.end());
    for (size_t j = 0; j < ids.size(); j++)
      EXPECT_EQ(ids[j], r[j]);
  }
}

TEST(SimQuadTree, sampleCellCenteredScalar) {
  SimQuadTree grid(BBox2D::make_unit_bbox(), 6);
  grid.setPointSetStructures(new PointZGrid2(8, 8, BBox2D::make_unit_bbox()),
                             new PointZGrid2(8, 8, BBox2D::make_unit_bbox()),
                             new PointZGrid2(8, 8, BBox2D::make_unit_bbox()));
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  {
    grid.setInterpolationMethod(SimulationDomain2::InterpolationMethod::RBF);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      ASSERT_NEAR(grid.sampleCellCenteredScalar(id, target), _f(target), 1e-2);
    }
  }
}

TEST(SimQuadTree, sampleFaceCenteredScalar) {
  SimQuadTree grid(BBox2D::make_unit_bbox(), 6);
  grid.setPointSetStructures(new PointZGrid2(8, 8, BBox2D::make_unit_bbox()),
                             new PointZGrid2(8, 8, BBox2D::make_unit_bbox()),
                             new PointZGrid2(8, 8, BBox2D::make_unit_bbox()));
  size_t id = grid.addFaceCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  grid.iterateDxScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  grid.iterateDyScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  grid.setBoundaryHandler(new DomainBoundaryHandler2(
      DomainBoundaryHandler2::BoundaryType::NEUMANN));
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::BILINEAR);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      Vector2d s = grid.sampleFaceCenteredScalar(id, target);
      ASSERT_NEAR(s[0], _f(target), 1e-2);
      ASSERT_NEAR(s[1], _f(target), 1e-2);
    }
  }
}

TEST(SimQuadTree, regularInterpolation) {
  SimQuadTree domain(BBox2D::make_unit_bbox(), 1);
  domain.setInterpolationMethod(
      SimulationDomain2::InterpolationMethod::BILINEAR);
  size_t id = domain.addFaceCenteredScalarField(2.);
  domain.iterateDxScalar(id, [](double &v, Point2d p) { v = p.y(); });
  domain.iterateDyScalar(id, [](double &v, Point2d p) { v = p.x(); });
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d()),
            Vector2d(0.25, 0.25));
  return;
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0, .25)),
            Vector2d(0.25, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0, .5)),
            Vector2d(0.5, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0.25, 0)),
            Vector2d(0.25, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(.5, 0)),
            Vector2d(0.25, 0.5));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(.5, .5)),
            Vector2d(0.5, 0.5));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0)),
            Vector2d(0.25, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0.25)),
            Vector2d(0.25, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0.5)),
            Vector2d(0.5, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 1)),
            Vector2d(0.75, 0.75));
}

TEST(SimQuadTree, RBFInterpolation) {
  SimQuadTree domain(BBox2D::make_unit_bbox(), 1);
  domain.setInterpolationMethod(SimulationDomain2::InterpolationMethod::RBF);
  size_t id = domain.addFaceCenteredScalarField(2.);
  domain.iterateDxScalar(id, [](double &v, Point2d p) { v = p.y(); });
  domain.iterateDyScalar(id, [](double &v, Point2d p) { v = p.x(); });

  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1., 1.)),
            Vector2d(0.75, 0.75));
  // return;
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d()),
            Vector2d(0.25, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0, .25)),
            Vector2d(0.25, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0, .5)),
            Vector2d(0.5, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(0.25, 0)),
            Vector2d(0.25, 0.25));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(.5, 0)),
            Vector2d(0.25, 0.5));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0)),
            Vector2d(0.25, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0.25)),
            Vector2d(0.25, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 0.5)),
            Vector2d(0.5, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(1, 1)),
            Vector2d(0.75, 0.75));
  EXPECT_EQ(domain.sampleFaceCenteredScalar(id, Point2d(.5, .5)),
            Vector2d(0.5, 0.5));
}
