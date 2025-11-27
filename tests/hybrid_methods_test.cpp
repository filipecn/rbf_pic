#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(HybridMethods, ApplyExternalForces) {
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 5);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(0.);
    Vector2d gravity(0., -9.8);
    double dt = 0.01;
    qt.iterateDxScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      v += dt * gravity.x();
    });
    qt.iterateDyScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      v += dt * gravity.y();
    });
    qt.iterateDxScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      EXPECT_EQ(IS_EQUAL(v, 0.), true);
    });
    qt.iterateDyScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      EXPECT_EQ(IS_EQUAL(v, dt * gravity.y()), true);
    });
  }
  {
    SimRegularDomain2 qt(10);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(0.);
    Vector2d gravity(0., -9.8);
    double dt = 0.01;
    qt.iterateDxScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      v += dt * gravity.x();
    });
    qt.iterateDyScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      v += dt * gravity.y();
    });
    qt.iterateDxScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      EXPECT_EQ(IS_EQUAL(v, 0.), true);
    });
    qt.iterateDyScalar(velocityField, [&](double &v, Point2d p) {
      UNUSED_VARIABLE(p);
      EXPECT_EQ(IS_EQUAL(v, dt * gravity.y()), true);
    });
  }
}

TEST(HybridMethods, CFL) {
  {
    SimRegularDomain2 qt(10);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(1.);
    Vector2d gravity(0., -9.8);
    double dt = 0.01;
    double maxVel = 0.;
    qt.iterateDxScalar(velocityField, [&](size_t i, double &v) {
      UNUSED_VARIABLE(i);
      maxVel = max(maxVel, std::fabs(v + dt * gravity.x()));
    });
    qt.iterateDyScalar(velocityField, [&](size_t i, double &v) {
      UNUSED_VARIABLE(i);
      maxVel = max(maxVel, std::fabs(v + dt * gravity.y()));
    });
    double cfl =
        maxVel / std::min(qt.smallestCellSize()[0], qt.smallestCellSize()[1]);
    EXPECT_EQ(cfl * dt < 1., true);
    ASSERT_NEAR(cfl * dt, dt / (1. / 10.), 1e-8);
  }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 6);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(1.);
    Vector2d gravity(0., -9.8);
    double dt = 0.01;
    double maxVel = 0.;
    qt.iterateDxScalar(velocityField, [&](size_t i, double &v) {
      UNUSED_VARIABLE(i);
      maxVel = max(maxVel, std::fabs(v + dt * gravity.x()));
    });
    qt.iterateDyScalar(velocityField, [&](size_t i, double &v) {
      UNUSED_VARIABLE(i);
      maxVel = max(maxVel, std::fabs(v + dt * gravity.y()));
    });
    double cfl =
        maxVel / std::min(qt.smallestCellSize()[0], qt.smallestCellSize()[1]);
    ASSERT_LT(cfl * dt, 1.);
    ASSERT_NEAR(cfl * dt, dt / (1. / std::pow(2., 6.)), 1e-8);
  }
}

TEST(HybridMethods, ApplyBoundaryVelocities) {
  {
    SimRegularDomain2 qt(20);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(1.);
    bh.iterateBoundaryFaces(
        [&](size_t i, DomainBoundaryHandler2::BoundaryType t) {
          UNUSED_VARIABLE(t);
          qt.faceCenteredScalar(velocityField, i) = 0.;
        });
    qt.iterateDxScalar(velocityField, [&](size_t i, double v) {
      if (bh.faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
        ASSERT_NEAR(v, 0., 1e-8);
      else
        ASSERT_NEAR(v, 1., 1e-8);
    });
    qt.iterateDyScalar(velocityField, [&](size_t i, double v) {
      if (bh.faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
        ASSERT_NEAR(v, 0., 1e-8);
      else
        ASSERT_NEAR(v, 1., 1e-8);
    });
  }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 6);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    size_t velocityField = qt.addFaceCenteredScalarField(1.);
    bh.iterateBoundaryFaces(
        [&](size_t i, DomainBoundaryHandler2::BoundaryType t) {
          UNUSED_VARIABLE(t);
          qt.faceCenteredScalar(velocityField, i) = 0.;
        });
    qt.iterateDxScalar(velocityField, [&](size_t i, double v) {
      if (bh.faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
        ASSERT_NEAR(v, 0., 1e-8);
      else
        ASSERT_NEAR(v, 1., 1e-8);
    });
    qt.iterateDyScalar(velocityField, [&](size_t i, double v) {
      if (bh.faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
        ASSERT_NEAR(v, 0., 1e-8);
      else
        ASSERT_NEAR(v, 1., 1e-8);
    });
  }
}

TEST(HybridMethods, AdvectParticles) {
  {
    SimRegularDomain2 qt(20);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    ParticleSystem2 ps(new PointZGrid2(64, 64, BBox2D::make_unit_bbox()),
                       new RBFInterpolant<Point2d>());
    double dt = 0.01;
    size_t vx = ps.addScalarProperty();
    size_t vy = ps.addScalarProperty();
    ps.addParticle(Point2d(.5, .5));
    ps.setScalarProperty(vx, 0, 1.);
    ps.setScalarProperty(vy, 0, -1.);
    EulerTemporalIntegrator2 eti;
    auto p = eti.integrate(
        ps[0],
        Vector2d(ps.getScalarProperty(vx, 0), ps.getScalarProperty(vy, 0)), dt);
    EXPECT_EQ(p, Point2d(.5 + dt * 1., .5 - dt * 1.));
    dt = 0.5;
    ps.addParticle(Point2d(0.9, 0.5));
    ps.setScalarProperty(vx, 1, 1.);
    p = eti.integrate(
        ps[1],
        Vector2d(ps.getScalarProperty(vx, 1), ps.getScalarProperty(vy, 1)), dt);
    EXPECT_EQ(p, Point2d(1.4, 0.5));
    auto i = bh.intersect(&qt, ps[1], p);
    EXPECT_EQ(i.isValid, true);
    auto r = Transform2D::applyReflection(p - i.point, i.normal.normalized());
    EXPECT_EQ(i.point + r, Point2d(0.6, 0.5));
    ps.setScalarProperty(vx, 1, 1.);
    ps.setScalarProperty(vy, 1, 1.);
    std::cout << ps.getScalarProperty(vy, 1) << std::endl;
    dt = 0.2;
    p = eti.integrate(
        ps[1],
        Vector2d(ps.getScalarProperty(vx, 1), ps.getScalarProperty(vy, 1)), dt);
    EXPECT_EQ(p, Point2d(1.1, 0.7));
    i = bh.intersect(&qt, ps[1], p);
    EXPECT_EQ(i.isValid, true);
    EXPECT_EQ(i.point, Point2d(1., 0.6));
    r = Transform2D::applyReflection(p - i.point, i.normal.normalized());
    EXPECT_EQ(i.point + r, Point2d(0.9, 0.7));
  }
}

TEST(HybridMethods, MarkCells) {
  {
    SimRegularDomain2 qt(10);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    ParticleSystem2 ps(new PointZGrid2(64, 64, BBox2D::make_unit_bbox()),
                       new RBFInterpolant<Point2d>());
    for (size_t i = 0; i < 5; i++)
      for (size_t j = 0; j < 10; j++)
        ps.addParticle(Point2d(j * 0.1 + 0.05, i * 0.1 + 0.05));
    bh.markCells(&qt, &ps);
    size_t airCount = 0, fluidCount = 0;
    for (size_t i = 0; i < qt.cellCount(); i++)
      if (bh.cellType(i) == DomainBoundaryHandler2::MaterialType::AIR)
        airCount++;
      else if (bh.cellType(i) == DomainBoundaryHandler2::MaterialType::FLUID)
        fluidCount++;
      else
        EXPECT_EQ(true, false);
    EXPECT_EQ(fluidCount, 50u);
    EXPECT_EQ(airCount, 50u);
  }
  {
    SimQuadTree qt(BBox2D::make_unit_bbox(), 6);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    ParticleSystem2 ps(new PointZGrid2(64, 64, BBox2D::make_unit_bbox()),
                       new RBFInterpolant<Point2d>());
    double h = 1. / 64.;
    for (size_t i = 0; i < 32; i++)
      for (size_t j = 0; j < 64; j++)
        ps.addParticle(Point2d(j * h + h * .5, i * h + h * .5));
    bh.markCells(&qt, &ps);
    size_t airCount = 0, fluidCount = 0;
    for (size_t i = 0; i < qt.cellCount(); i++)
      if (bh.cellType(i) == DomainBoundaryHandler2::MaterialType::AIR)
        airCount++;
      else if (bh.cellType(i) == DomainBoundaryHandler2::MaterialType::FLUID)
        fluidCount++;
      else
        EXPECT_EQ(true, false);
    EXPECT_EQ(fluidCount, 64u * 32u);
    EXPECT_EQ(airCount, 64u * 32u);
  }
}
