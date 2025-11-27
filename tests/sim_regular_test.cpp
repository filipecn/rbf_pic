#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(SimRegularDomain2, Constructors) {
  {
    SimRegularDomain2 grid;
    EXPECT_EQ(grid.cellCount(), 64u);
  }
  {
    SimRegularDomain2 grid(4u);
    EXPECT_EQ(grid.cellCount(), 16u);
  }
}

TEST(SimRegularDomain2, RandomAccessors) {
  SimRegularDomain2 gd(10);
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

TEST(SimRegularDomain2, idFunctions) {
  SimRegularDomain2 g(4);
  EXPECT_EQ(g.faceIdFromCellId(5, SimulationDomain2::Orientation::LEFT), 6);
  EXPECT_EQ(g.faceIdFromCellId(5, SimulationDomain2::Orientation::RIGHT), 7);
  EXPECT_EQ(g.faceIdFromCellId(5, SimulationDomain2::Orientation::BOTTOM), 25);
  EXPECT_EQ(g.faceIdFromCellId(5, SimulationDomain2::Orientation::TOP), 29);
  EXPECT_EQ(g.cellIdFromFaceId(5, SimulationDomain2::Orientation::LEFT), -1);
  EXPECT_EQ(g.cellIdFromFaceId(5, SimulationDomain2::Orientation::RIGHT), 4);
  EXPECT_EQ(g.cellIdFromFaceId(25, SimulationDomain2::Orientation::BOTTOM), 1);
  EXPECT_EQ(g.cellIdFromFaceId(25, SimulationDomain2::Orientation::TOP), 5);
  for (unsigned i = 0; i < 4; i++)
    for (unsigned j = 0; j < 4; j++)
      EXPECT_EQ(g.ijFromCellId(j * 4 + i), Point2i(i, j));
  for (unsigned i = 0; i < 5; i++)
    for (unsigned j = 0; j < 4; j++)
      EXPECT_EQ(g.ijFromFaceId(j * 5 + i), Point2i(i, j));
  for (unsigned i = 0; i < 4; i++)
    for (unsigned j = 0; j < 5; j++)
      EXPECT_EQ(g.ijFromFaceId(5 * 4 + j * 4 + i), Point2i(i, j));
  EXPECT_EQ(g.uFaceIdFromIJ(Point2i(0, 0)), 0u);
  EXPECT_EQ(g.uFaceIdFromIJ(Point2i(3, 3)), 18u);
  EXPECT_EQ(g.vFaceIdFromIJ(Point2i(0, 0)), 20u);
  EXPECT_EQ(g.vFaceIdFromIJ(Point2i(3, 3)), 35u);
}

TEST(SimRegularDomain2, CellRegion) {
  SimRegularDomain2 rd(10);
  auto r = rd.cellRegion(5);
  EXPECT_EQ(Point2d(5 * 0.1 - 0.05 + 0.05, 0.), r.lower());
  EXPECT_EQ(Point2d(5 * 0.1 + 0.05 + 0.05, 0.1), r.upper());
}

TEST(SimRegularDomain2, Transforms) {
  SimRegularDomain2 grid(10);
  for (size_t c = 0; c < grid.cellCount(); c++) {
    size_t i = c % 10;
    size_t j = c / 10;
    EXPECT_EQ(grid.cellCenterPosition(c),
              Point2d(i * 0.1 + 0.05, j * 0.1 + 0.05));
  }
  for (size_t f = 0; f < grid.faceCount(); f++) {
    if (f < 110) {
      size_t i = f % 11;
      size_t j = f / 11;
      EXPECT_EQ(grid.faceCenterPosition(f), Point2d(i * 0.1, j * 0.1 + 0.05));
    } else {
      size_t i = (f - 110) % 10;
      size_t j = (f - 110) / 10;
      EXPECT_EQ(grid.faceCenterPosition(f), Point2d(i * 0.1 + 0.05, j * 0.1));
    }
  }
}

TEST(SimRegularDomain2, IterateCell) {
  SimRegularDomain2 grid;
  unsigned id = grid.addCellCenteredScalarField(1.0);
  unsigned k = 0;
  double h = 1. / 8.;
  grid.iterateCellScalar(id, [&](double &v, Point2d p) {
    EXPECT_NEAR(v, 1.0, 1e-8);
    auto ij = grid.ijFromCellId(k);
    EXPECT_EQ(p, grid.cellCenterPosition(ij[1] * 8 + ij[0]));
    ASSERT_NEAR(ij[0] * h + 1. / 16., p.x(), 1e-8);
    ASSERT_NEAR(ij[1] * h + 1. / 16., p.y(), 1e-8);
    k++;
  });
  EXPECT_EQ(k, grid.cellCount());
}

TEST(SimRegularDomain2, IterateFaces) {
  SimRegularDomain2 grid(4u);
  size_t cellId = grid.addCellCenteredScalarField(1.0);
  size_t faceId = grid.addFaceCenteredScalarField(1.0);
  size_t k = 0;
  grid.iterateCellScalar(cellId, [&](double &v, Point2d p) {
    auto ij = grid.ijFromCellId(k++);
    EXPECT_EQ(v, 1.0);
    EXPECT_EQ(p, Point2d(ij[0] * 0.25 + 0.125, ij[1] * 0.25 + 0.125));
  });
  k = 0;
  grid.iterateDxScalar(faceId, [&](double &v, Point2d p) {
    auto ij = grid.ijFromFaceId(k++);
    EXPECT_EQ(v, 1.0);
    EXPECT_EQ(p, Point2d(ij[0] * 0.25, ij[1] * 0.25 + 0.125));
  });
  grid.iterateDyScalar(faceId, [&](double &v, Point2d p) {
    auto ij = grid.ijFromFaceId(k++);
    EXPECT_EQ(v, 1.0);
    EXPECT_EQ(p, Point2d(ij[0] * 0.25 + 0.125, ij[1] * 0.25));
  });
}

TEST(SimRegularDomain2, FaceNeibors) {
  SimRegularDomain2 grid(4);
  DomainBoundaryHandler2 bh =
      DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  grid.setBoundaryHandler(&bh);
  bh.markFaces(&grid);

  for (size_t i = 0; i < 4; i++) {
    auto ng = grid.faceNeighborhood(i);
    // std::cout << "Face " << i << ": " << '\n';
    // std::cout << ng[0];
    // std::cout << ng[1] << '\n';
  }
}

TEST(SimRegularDomain2, CellStencils) {
  SimRegularDomain2 grid(4);
  DomainBoundaryHandler2 bh =
      DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  grid.setBoundaryHandler(&bh);
  bh.markFaces(&grid);
  // bh.printFaces();
  for (size_t i = 0; i < grid.cellCount(); i++) {
    double wsum = 0.;
    int k = 0;
    grid.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          UNUSED_VARIABLE(o);
          UNUSED_VARIABLE(p);
          UNUSED_VARIABLE(id);
          //   std::cout << id << ' ';
          //   EXPECT_EQ(expectedIds[i][k], id);
          //   EXPECT_EQ(expectedOrientation[i][k], o);
          //   EXPECT_EQ(IS_EQUAL_ERROR(w, expectedWeights[k], 1e-4), true);
          //   EXPECT_EQ(p == expectedCenters[i][k], true);
          wsum += w;
          //   std::cout << id << ": " << w << '\n';
          k++;
        });
    // std::cout << '\n';
    // ASSERT_NEAR(wsum, 0., 1e-8);
  }
}

TEST(SimRegularDomain2, CellBoundaryStencil) {
  SimRegularDomain2 qt(5);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  EXPECT_EQ(bh.dirichletBoundaryCount(), 10u);
  int k = 0, kk = 0;
  for (size_t i = 0; i < qt.cellCount(); i++)
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0 && (o == SimulationDomain2::Orientation::LEFT ||
              o == SimulationDomain2::Orientation::RIGHT))
            k++;
          UNUSED_VARIABLE(p);
          UNUSED_VARIABLE(w);
        });
  for (size_t i = 0; i < qt.cellCount(); i++)
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, DomainBoundaryHandler2::BoundaryType bt,
               double v) {
          if (id < 0 && bt == DomainBoundaryHandler2::BoundaryType::DIRICHLET)
            kk++;
          UNUSED_VARIABLE(v);
          UNUSED_VARIABLE(w);
        });
  EXPECT_EQ(k, 10);
  EXPECT_EQ(k, kk);
}

TEST(SimRegularDomain2, gradientXAtCellCenteredScalar) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fx = [](Point2d p) {
    return -std::sin(p[0]) * std::cos(p[1]);
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  for (size_t i = 0; i < grid.cellCount(); i++) {
    ASSERT_NEAR(_fx(grid.cellCenterPosition(i)),
                grid.gradientXAtCellCenteredScalar(id, i), 1e-2);
  }
}

TEST(SimRegularDomain2, gradientYAtCellCenteredScalar) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fy = [](Point2d p) {
    return -std::cos(p[0]) * std::sin(p[1]);
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  for (size_t i = 0; i < grid.cellCount(); i++) {
    ASSERT_NEAR(_fy(grid.cellCenterPosition(i)),
                grid.gradientYAtCellCenteredScalar(id, i), 1e-2);
  }
}

TEST(SimRegularDomain2, gradientXYAtCellCenteredScalar) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fxy = [](Point2d p) {
    return std::sin(p[0]) * std::sin(p[1]);
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  for (size_t i = 0; i < grid.cellCount(); i++)
    ASSERT_NEAR(_fxy(grid.cellCenterPosition(i)),
                grid.gradientXYAtCellCenteredScalar(id, i), 1e-2);
}

TEST(SimRegularDomain2, faceDivergentAtCellCenter) {
  SimRegularDomain2 grid(100);
  // size_t id = grid.addCellCenteredScalarField(0.);
  size_t id = grid.addFaceCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fx = [](Point2d p) {
    return -std::sin(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fy = [](Point2d p) {
    return -std::cos(p[0]) * std::sin(p[1]);
  };
  grid.iterateDxScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  grid.iterateDyScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  // grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  for (size_t i = 0; i < grid.cellCount(); i++) {
    ASSERT_NEAR(_fx(grid.cellCenterPosition(i)) +
        _fy(grid.cellCenterPosition(i)),
                grid.faceDivergentAtCellCenter(id, i), 1e-5);
  }
}

TEST(SimRegularDomain2, cellGradientAtFaceCenter) {
  SimRegularDomain2 grid(30);
  size_t idcell = grid.addCellCenteredScalarField(0.);
  // size_t idface = grid.addFaceCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0] * PI) * std::cos(p[1] * PI);
  };
  std::function<double(Point2d)> _fx = [](Point2d p) {
    return -PI * std::sin(p[0] * PI) * std::cos(p[1] * PI);
  };
  std::function<double(Point2d)> _fy = [](Point2d p) {
    return -PI * std::cos(p[0] * PI) * std::sin(p[1] * PI);
  };

  // grid.iterateDxScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  // grid.iterateDyScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  grid.iterateCellScalar(idcell, [&](double &v, Point2d p) { v = _f(p); });
  for (size_t i = 0; i < grid.faceCount(); i++) {
    if (i < grid.faceCount() / 2)
      ASSERT_NEAR(_fx(grid.faceCenterPosition(i)),
                  grid.cellGradientAtFaceCenter(idcell, i), 1e-2);
    else
      ASSERT_NEAR(_fy(grid.faceCenterPosition(i)),
                  grid.cellGradientAtFaceCenter(idcell, i), 1e-2);
  }
}

TEST(SimRegularDomain2, CenteredGridCell) {
  {
    SimRegularDomain2 g(2);
    unsigned v[4];

    EXPECT_EQ(g.cellCenteredGridCell(Point2d(0, 0), &v[0]), Point2d(0, 0));

    EXPECT_EQ(g.cellCenteredGridCell(Point2d(0.5, 0.5), &v[0]), Point2d(0.5, 0.5));

    EXPECT_EQ(g.cellCenteredGridCell(Point2d(0.5, 0.), &v[0]), Point2d(0.5, 0.));
  }
  {
    SimRegularDomain2 g(2);
    auto id = g.addFaceCenteredScalarField(10.);
    unsigned v[4];
    double f[4];
    EXPECT_EQ(g.uFaceCenteredGridCell(Point2d(1, 1), &v[0], &f[0], id), Point2d(0, 0));
    
    EXPECT_EQ(g.uFaceCenteredGridCell(Point2d(0, 0), &v[0], &f[0], id), Point2d(0., 0.));
    
    EXPECT_EQ(g.uFaceCenteredGridCell(Point2d(0.5, 0.5), &v[0], &f[0], id), Point2d(0.5, 0.5));
    
    EXPECT_EQ(g.uFaceCenteredGridCell(Point2d(0.5, 0.), &v[0], &f[0], id),Point2d(0.5, 0.));
    
    EXPECT_EQ(g.vFaceCenteredGridCell(Point2d(0,0), &v[0], &f[0], id), Point2d(0., 0.));
    
    EXPECT_EQ(g.vFaceCenteredGridCell(Point2d(0.5, 0.5), &v[0], &f[0], id),Point2d(0.5, 0.5));

    EXPECT_EQ(g.vFaceCenteredGridCell(Point2d(0.5, 0.), &v[0], &f[0], id), Point2d(0.5, 0.));
  }
}

TEST(SimRegularDomain2, sampleCellCenteredScalar) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::BILINEAR);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      ASSERT_NEAR(grid.sampleCellCenteredScalar(id, target), _f(target), 1e-2);
    }
  }
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::BICUBIC);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      ASSERT_NEAR(grid.sampleCellCenteredScalar(id, target), _f(target), 1e-2);
    }
  }
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::CATMULL);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      ASSERT_NEAR(grid.sampleCellCenteredScalar(id, target), _f(target), 1e-2);
    }
  }
}

TEST(SimRegularDomain2, sampleScalarUsingCatmullRom) {
  size_t n = 5;
  SimRegularDomain2 grid(n);
  size_t id = grid.addCellCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [&](Point2d p) {
    // Get the index instead of the normalized float coordinates.
    return std::floor(p.x() / (1./n)) + std::floor(p.y() / (1./n)) + 1.;
  };
  grid.iterateCellScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::CATMULL);
      Point2d target(1.5, 2.8);
      double result = grid.sampleCellCenteredScalar(id, target);
      EXPECT_LT(4.0, result);
      EXPECT_GT(6.0, result);
  }
}

TEST(SimRegularDomain2, sampleFaceCenteredScalar) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addFaceCenteredScalarField(0.);
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  grid.iterateDxScalar(id, [&](double &v, Point2d p) { v = _f(p); });
  grid.iterateDyScalar(id, [&](double &v, Point2d p) { v = _f(p); });
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
  {
    Point2d p0(0, 0);
    Point2d p1(0.5, 0.5);
    Point2d p2(1, 1);
    Point2d p3(0.25, 0.25);

    Vector2d s = grid.sampleFaceCenteredScalar(id, p0);
    ASSERT_NEAR(s[0], _f(p0), 1e-2);
    ASSERT_NEAR(s[1], _f(p0), 1e-2);

    s = grid.sampleFaceCenteredScalar(id, p1);
    ASSERT_NEAR(s[0], _f(p1), 1e-2);
    ASSERT_NEAR(s[1], _f(p1), 1e-2);

    s = grid.sampleFaceCenteredScalar(id, p2);
    ASSERT_NEAR(s[0], _f(p2), 1e-2);
    ASSERT_NEAR(s[1], _f(p2), 1e-2);

    s = grid.sampleFaceCenteredScalar(id, p3);
    ASSERT_NEAR(s[0], _f(p3), 1e-2);
    ASSERT_NEAR(s[1], _f(p3), 1e-2);
  }
  {
    grid.setInterpolationMethod(
        SimulationDomain2::InterpolationMethod::BICUBIC);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      Point2d target = sampler.sample(BBox2D::make_unit_bbox());
      Vector2d s = grid.sampleFaceCenteredScalar(id, target);
      ASSERT_NEAR(s[0], _f(target), 1e-2);
      ASSERT_NEAR(s[1], _f(target), 1e-2);
    }
  }
}

TEST(SimRegularDomain2, extrapolateStaggeredField) {
  SimRegularDomain2 grid(100);
  size_t id = grid.addFaceCenteredScalarField(0);
  DomainBoundaryHandler2 bh =
      DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  grid.setBoundaryHandler(&bh);
  for (int i = 0; i < 100; ++i)
    for (int j = 0; j < 100; ++j) {
      unsigned cell = j * 100 + i;

      if (j > 50)
        bh.setCellType(cell, furoo::DomainBoundaryHandler2::MaterialType::AIR);
      else
        bh.setCellType(cell, furoo::DomainBoundaryHandler2::MaterialType::FLUID);
    }
  grid.setInterpolationMethod(
      SimulationDomain2::InterpolationMethod::BILINEAR);

  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };

  grid.iterateDxScalar(
      id, [&](double &v, Point2d p) { v = (p[1] > 0.5 ? 0 : _f(p)); });
  grid.iterateDyScalar(
      id, [&](double &v, Point2d p) { v = (p[1] > 0.5 ? 0 : _f(p)); });

  for (size_t i = 0; i < 10000; i++) {
    BoxSampler sampler;
    Point2d target = sampler.sample(BBox2D::make_unit_bbox());
    Vector2d s = grid.sampleFaceCenteredScalar(id, target);
    ASSERT_NEAR(s[0], _f(target), 1e-2);
    ASSERT_NEAR(s[1], _f(target), 1e-2);
  }

  Point2d pp(0.55, 0.65);
  Vector2d s = grid.sampleFaceCenteredScalar(id, pp);
  ASSERT_NEAR(s[0], 0, 1e-8);
  ASSERT_NEAR(s[1], 0, 1e-8);

  grid.extrapolateToDomain(id);

  for (size_t i = 0; i < 10000; i++) {
    BoxSampler sampler;
    Point2d target = sampler.sample(BBox2D::make_unit_bbox());
    Vector2d s = grid.sampleFaceCenteredScalar(id, target);
    ASSERT_NEAR(s[0], _f(target), 1e-2);
    ASSERT_NEAR(s[1], _f(target), 1e-2);
  }
}

TEST(SimRegularDomain2, extrapolateStaggeredField2) {
  SimRegularDomain2 grid(5);
  size_t id = grid.addFaceCenteredScalarField(0);
  DomainBoundaryHandler2 bh =
      DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  grid.setBoundaryHandler(&bh);

  for (size_t cell = 0; cell < grid.cellCount(); cell++)
    bh.setCellType(cell, furoo::DomainBoundaryHandler2::MaterialType::AIR);
  bh.setCellType(0, furoo::DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(24, furoo::DomainBoundaryHandler2::MaterialType::FLUID);
  bh.setCellType(12, furoo::DomainBoundaryHandler2::MaterialType::FLUID);
  grid.scalarAtFace(0, 0) = 1.;
  grid.scalarAtFace(0, 1) = 1.;
  grid.scalarAtFace(0, 28) = -1.;
  grid.scalarAtFace(0, 29) = -1.;
  grid.scalarAtFace(0, 30) = 1.;
  grid.scalarAtFace(0, 35) = 1.;
  grid.scalarAtFace(0, 54) = -1.;
  grid.scalarAtFace(0, 59) = -1.;
  grid.scalarAtFace(0, 14) = 2.;
  grid.scalarAtFace(0, 15) = 2.;
  grid.scalarAtFace(0, 42) = 2.;
  grid.scalarAtFace(0, 47) = 2.;

  for (int i = 4; i >= 0; i--) {
    for (size_t j = 0; j < 6; j++)
      std::cout << grid.faceCenteredScalar(0, i * 6 + j) << " ";
    std::cout << std::endl;
  }
  for (int i = 5; i >= 0; i--) {
    for (size_t j = 0; j < 5; j++)
      std::cout << grid.faceCenteredScalar(0, i * 5 + j + 30) << " ";
    std::cout << std::endl;
  }
  grid.extrapolateToDomain(id);
  std::cout << std::endl;
  for (int i = 4; i >= 0; i--) {
    for (size_t j = 0; j < 6; j++)
      std::cout << grid.faceCenteredScalar(0, i * 6 + j) << " ";
    std::cout << std::endl;
  }
  for (int i = 5; i >= 0; i--) {
    for (size_t j = 0; j < 5; j++)
      std::cout << grid.faceCenteredScalar(0, i * 5 + j + 30) << " ";
    std::cout << std::endl;
  }
}

struct SillyFunctionOp
{
  SillyFunctionOp() = default;
  Vector2d operator()(const Vector2d& xyz) const
  {
    return Vector2d(3. * xyz.y() + 1., 5. * xyz.x() + 7.);
  }
};

TEST(SimRegularDomain2, SampleStaggeredUsingBilinearInterpolation)
{
  SimRegularDomain2 grid(8);
  grid.setInterpolationMethod(furoo::SimulationDomain2::InterpolationMethod::BILINEAR);
  size_t fieldId = grid.addFaceCenteredScalarField(0.);

  SillyFunctionOp func;

  for (int i = 0; i <= 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      grid.setScalarAtFace(fieldId, 
                           grid.uFaceIdFromIJ(Point2i(i, j)), 
                           func(Vector2d(i, j)).x());
    }
  }

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j <= 8; ++j)
    {
      grid.setScalarAtFace(fieldId, 
                           grid.vFaceIdFromIJ(Point2i(i, j)), 
                           func(Vector2d(i, j)).y());
    }
  }

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      EXPECT_NEAR(
        grid.sampleFaceFieldAtCellCenter(fieldId, 
                                         grid.cellIdFromIJ(Point2i(i, j))).x(),
        func(Vector2d(i, j)).x(),
        1e-6);

      EXPECT_NEAR(
        grid.sampleFaceFieldAtCellCenter(fieldId, 
                                         grid.cellIdFromIJ(Point2i(i, j))).y(),
        func(Vector2d(i, j)).y(),
        1e-6);
    }
  }
}

TEST(SimRegularDomain2, SampleScalarUsingBicubicInterpolation)
{
  SimRegularDomain2 grid(8);
  grid.setInterpolationMethod(furoo::SimulationDomain2::InterpolationMethod::BICUBIC);
  size_t fieldId = grid.addCellCenteredScalarField(0.);

  SillyFunctionOp func;

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      grid.setScalarAtCell(fieldId, 
                           grid.cellIdFromIJ(Point2i(i, j)), 
                           func(Vector2d(i, j)).y());
    }
  }

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      EXPECT_NEAR(
        grid.sampleCellCenteredScalar(fieldId, 
                                      Point2d(i, j)),
        func(Vector2d(i, j)).y(),
        1e-6);
    }
  }
}

TEST(SimRegularDomain2, SampleStaggeredUsingBicubicInterpolation)
{
  SimRegularDomain2 grid(8);
  grid.setInterpolationMethod(furoo::SimulationDomain2::InterpolationMethod::BICUBIC);
  size_t fieldId = grid.addFaceCenteredScalarField(0.);

  SillyFunctionOp func;

  for (int i = 0; i <= 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      grid.setScalarAtFace(fieldId, 
                           grid.uFaceIdFromIJ(Point2i(i, j)), 
                           func(Vector2d(i, j)).x());
    }
  }

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j <= 8; ++j)
    {
      grid.setScalarAtFace(fieldId, 
                           grid.vFaceIdFromIJ(Point2i(i, j)), 
                           func(Vector2d(i, j)).y());
    }
  }

  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      EXPECT_NEAR(
        grid.sampleUFaceField(fieldId, Point2d(i, j)),
        func(Vector2d(i, j)).x(),
        1e-6);

      EXPECT_NEAR(
        grid.sampleVFaceField(fieldId, Point2d(i, j)),
        func(Vector2d(i, j)).y(),
        1e-6);    
    }
  }
}
