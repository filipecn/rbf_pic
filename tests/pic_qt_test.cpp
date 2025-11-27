#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

class PICSolverTestClass : public PICSolver2 {
public:
  PICSolverTestClass() {}

  PICSolverTestClass(DomainBoundaryHandler2 *b, SimulationDomain2 *d,
                     ParticleSystem2 *ps)
      : PICSolver2(b, d, ps) {}

  void runOnInitialize() { this->onInitialize(); }

  void runApplyBoundaryCondition() { this->applyBoundaryCondition(); }

  void runBeginAdvanceTimeStep(double dt) { this->beginAdvanceTimeStep(dt); }

  void runComputeExternalForces(double dt) { this->computeExternalForces(dt); }

  void runComputeViscosity(double dt) { this->computeViscosity(dt); }

  void runComputePressure(double dt) { this->computePressure(dt); }

  void runCorrectVelocities(double dt) { this->correctVelocities(dt); }

  void runTransferFromParticlesToGrid() { this->transferFromParticlesToGrid(); }

  void runTransferFromGridToParticles() { this->transferFromGridToParticles(); }

  void runComputeAdvection(double dt) { this->computeAdvection(dt); }

  void runEndAdvanceTimeStep(double dt) { this->endAdvanceTimeStep(dt); }

  void runExtrapolateToDomain() { this->propagateVelocitiesToDomain(); }
};

TEST(PICSolverQT, Constructors) {
  {
    auto domain = new SimQuadTree(BBox2D::make_unit_bbox(), 3);
    auto ps = new ParticleSystem2(new PointZGrid2(16));
    domain->setBoundaryHandler(new DomainBoundaryHandler2(
        DomainBoundaryHandler2::BoundaryType::NEUMANN));
    PICSolver2 pic;
    pic.setSimulationDomain(domain);
    pic.setParticleSystem(ps);
  }
  {
    auto domain = new SimRegularDomain2(10);
    auto ps = new ParticleSystem2(new PointZGrid2(16));
    PICSolver2 pic(new DomainBoundaryHandler2(
        DomainBoundaryHandler2::BoundaryType::NEUMANN),
                   domain, ps);
  }
}

TEST(PICSolverQT, Frame) {
  auto domain = new SimQuadTree(BBox2D::make_unit_bbox(), 1);
  auto ps = new ParticleSystem2(new PointZGrid2(16));
  PICSolver2 pic(
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN),
      domain, ps);
  auto f = pic.currentFrame();
  EXPECT_EQ(f.index(), 0u);
  ASSERT_NEAR(f.ellapsedTime(), 0., 1e-8);
  ASSERT_NEAR(pic.currentTime(), 0., 1e-8);
  pic.advanceFrame();
  f = pic.currentFrame();
  EXPECT_EQ(f.index(), 1u);
  ASSERT_NEAR(f.ellapsedTime(), 1. / 60., 1e-8);
}

TEST(PICSolverQT, ComputeExternalForces) {
  auto domain = new SimQuadTree(BBox2D::make_unit_bbox(), 4);
  auto ps = new ParticleSystem2(new PointZGrid2(16));
  PICSolverTestClass pic(
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN),
      domain, ps);
  double dt = 0.01;
  pic.runOnInitialize();
  pic.runComputeExternalForces(dt);
  domain->iterateDxScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    ASSERT_NEAR(v, 0., 1e-8);
  });
  domain->iterateDyScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    ASSERT_NEAR(v, -9.8 * dt, 1e-8);
  });
}

TEST(PICSolverQT, ApplyBoundaryCondition) {
  auto domain = new SimQuadTree(BBox2D::make_unit_bbox(), 3);
  auto ps = new ParticleSystem2(new PointZGrid2(16));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  domain->setBoundaryHandler(bh);
  PICSolverTestClass pic(bh, domain, ps);
  pic.runOnInitialize();
  domain->iterateDxScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = 1.;
  });
  domain->iterateDyScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = 1.;
  });
  pic.runApplyBoundaryCondition();
  domain->iterateDxScalar(0, [&](size_t i, double v) {
    if (bh->faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
      ASSERT_NEAR(v, 0., 1e-8);
    else
      ASSERT_NEAR(v, 1., 1e-8);
  });
  domain->iterateDyScalar(0, [&](size_t i, double v) {
    if (bh->faceBoundaryType(i) != DomainBoundaryHandler2::BoundaryType::NONE)
      ASSERT_NEAR(v, 0., 1e-8);
    else
      ASSERT_NEAR(v, 1., 1e-8);
  });
}

TEST(PICSolverQT, MarkCells) {
  auto region = BBox2D::make_unit_bbox();
  auto domain = new SimQuadTree(BBox2D::make_unit_bbox(), 4);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh = new DomainBoundaryHandler2();
  PICSolverTestClass pic(bh, domain, ps);
  pic.runOnInitialize();
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 10; j++)
      ps->addParticle(Point2d(j * 0.1 + 0.05, i * 0.1 + 0.05));
  bh->markCells(domain, ps);
  size_t airCount = 0, fluidCount = 0;
  for (size_t i = 0; i < domain->cellCount(); i++)
    if (bh->cellType(i) == DomainBoundaryHandler2::MaterialType::AIR)
      airCount++;
    else if (bh->cellType(i) == DomainBoundaryHandler2::MaterialType::FLUID)
      fluidCount++;
    else
      EXPECT_EQ(true, false);
  EXPECT_EQ(fluidCount, 50u);
  EXPECT_EQ(airCount, 206u);
}

TEST(PICSolverQT, ComputeAdvection) {
  auto region = BBox2D::make_unit_bbox();
  auto domain = new SimQuadTree(region, 4);
  domain->setPointSetStructures(new PointZGrid2(16, 16, region),
                                new PointZGrid2(16, 16, region),
                                new PointZGrid2(16, 16, region));
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh = new DomainBoundaryHandler2();
  PICSolverTestClass pic(bh, domain, ps);
  pic.runOnInitialize();
  BoxSampler s;
  std::vector<Point2d> points;
  std::vector<Vector2d> velocities;
  for (size_t i = 0; i < 50; i++) {
    auto p = s.sample(region);
    auto pv = s.sample();
    auto v = Vector2d(pv.x(), pv.y());
    ps->addParticle(p);
    ps->setScalarProperty(0, i, v.x());
    ps->setScalarProperty(1, i, v.y());
    points.emplace_back(p);
    velocities.emplace_back(v);
  }
  double dt = 0.5;
  ps->iterateParticles([&](unsigned int i, Point2d p) {
    EXPECT_EQ(points[i], p);
    EXPECT_EQ(velocities[i], Vector2d(ps->getScalarProperty(0, i),
                                      ps->getScalarProperty(1, i)));
  });
  // std::cout << "Advecting" << '\n';
  pic.runComputeAdvection(dt);
  // advect here
  // std::cout << "NOW\n";
  for (size_t i = 0; i < 50; i++) {
    EulerTemporalIntegrator2 eti;
    // std::cout << "- integrate " << i << std::endl
    //   << points[i] << " with "
    //   << Vector2d(velocities[i].x(), velocities[i].y());
    auto p = eti.integrate(points[i],
                           Vector2d(velocities[i].x(), velocities[i].y()), dt);
    // std::cout << p << std::endl;
    auto in = bh->intersect(domain, points[i], p);
    if (in.isValid) {
      //   std::cout << "reflect" << '\n';
      auto r =
          Transform2D::applyReflection(p - in.point, in.normal.normalized());
      p = in.point + r;
    }
    points[i] = p;
  }
  // ps->size();
  ps->iterateParticles([&](unsigned int i, Point2d p) {
    // std::cout << i << std::endl;
    UNUSED_VARIABLE(p);
    EXPECT_EQ(points[i], p);
  });
}
