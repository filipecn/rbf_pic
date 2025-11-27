#include <furoo.h>
#include <gtest/gtest.h>
#include <fstream>

using namespace furoo;

void saveMatrixToMatlab(const LinearMatrix& matrix, const std::string& filename)
{
  int sizeX = (int)matrix.rowCount();
  int sizeY = (int)matrix.columnCount();

  std::cout << "matrix size here:: " << matrix.rowCount() << " " << matrix.columnCount() << std::endl;
  /// build string
  std::ostringstream oss;
  oss << "furoo_matrix =[";

  for (int i = 0; i < sizeX; ++i) {
    for (int j = 0; j < sizeY; ++j) {
      oss << " " << matrix(i, j);
    }
    oss << ";";
  }
  oss << "];" << std::endl;

  std::ofstream file;
  file.open(filename);

  file << oss.str();

  file.close();
}

void saveVectorToMatlab(const LinearVector& vector, const std::string& filename, const std::string& vectorName)
{
    int size = (int)vector.size();

    /// build string
    std::ostringstream oss;
    oss << vectorName;
    oss << " = [";

    for (int i = 0; i < size; ++i) {
            oss << " " << vector[i];
        }

    oss << "];" << std::endl;

    std::ofstream file;
    file.open(filename);

    file << oss.str();
    file.close();
}

void saveLinearSystem(const LinearMatrix& matrix, const LinearVector& vector,
                      const LinearVector& solution,
                      const LinearVector& divergent,
                      int frame_index)
{
  std::string file = "matrix." + std::to_string(frame_index) + ".m";
  saveMatrixToMatlab(matrix, file);

  file = "rhs." + std::to_string(frame_index) + ".m";
  saveVectorToMatlab(vector, file, "furoo_rhs");

  file = "solution." + std::to_string(frame_index) + ".m";
  saveVectorToMatlab(solution, file, "furoo_sol");

  file = "divergent." + std::to_string(frame_index) + ".m";
  saveVectorToMatlab(divergent, file, std::string("furoo_div" + std::to_string(frame_index)));
}


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

class PoissonDomain : public SimRegularDomain2 {
public:
  PoissonDomain(size_t n) : SimRegularDomain2(n) {}

  double faceDivergentAtCellCenter(unsigned int vf,
                                   unsigned int fid) const override {
    UNUSED_VARIABLE(vf);
    Point2d cp = this->cellCenterPosition(fid);
    return -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
  }
};

class PoissonQTDomain : public SimQuadTree {
public:
  PoissonQTDomain(BBox2D b, size_t n) : SimQuadTree(b, n) {}

  double faceDivergentAtCellCenter(unsigned int vf,
                                   unsigned int fid) const override {
    UNUSED_VARIABLE(vf);
    Point2d cp = this->cellCenterPosition(fid);
    return -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
  }
};

TEST(PICSolver2, Constructors) {
  {
    auto domain = new SimRegularDomain2(10);
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

TEST(PICSolver2, Frame) {
  auto domain = new SimRegularDomain2(10);
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

TEST(PICSolver2, ComputeExternalForces) {
  auto domain = new SimRegularDomain2(10);
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

TEST(PICSolver2, ApplyBoundaryCondition) {
  auto domain = new SimRegularDomain2(10);
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

TEST(PICSolver2, MarkCells) {
  auto region = BBox2D::make_unit_bbox();
  auto domain = new SimRegularDomain2(10);
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
  EXPECT_EQ(airCount, 50u);
}

TEST(PICSolver2, ComputePressure_poisson) {
  return;
  {
    auto domain = new PoissonDomain(16);
    auto ps = new ParticleSystem2(
        new PointZGrid2(16u, 16u, BBox2D::make_unit_bbox()));
    auto bh = new DomainBoundaryHandler2(
        DomainBoundaryHandler2::BoundaryType::DIRICHLET,
        DomainBoundaryHandler2::BoundaryType::NEUMANN,
        DomainBoundaryHandler2::BoundaryType::DIRICHLET,
        DomainBoundaryHandler2::BoundaryType::NEUMANN);
    PICSolverTestClass pic(bh, domain, ps);
    pic.runOnInitialize();
    bh->markFaces(domain);
    for (size_t i = 0; i < domain->cellCount(); i++)
      bh->setCellType(i, DomainBoundaryHandler2::MaterialType::FLUID);
    bh->iterateBoundaryFaces(
        [&](size_t i, DomainBoundaryHandler2::BoundaryType bt) {
          if (bt == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
            Point2d p = domain->faceCenterPosition(i);
            bh->setBoundaryFaceValue(i, cos(PI * p.x()) * cos(PI * p.y()));
          }
        });
    pic.runComputePressure(1);
    LinearVector error;
    domain->iterateCellScalar(0, [&](size_t i, double v) {
      Point2d cp = domain->cellCenterPosition(i);
      double s = cos(PI * cp.x()) * cos(PI * cp.y());
      ASSERT_NEAR(s, v, 1e-2);
      error.emplace_back(v - s);
    });
    std::cout << "ERROR " << Blas::norm(&error) / domain->cellCount()
              << std::endl;
    EXPECT_EQ(Blas::norm(&error) / domain->cellCount() < 1e-3, true);
  }
  {
    auto domain = new PoissonQTDomain(BBox2D::make_unit_bbox(), 4);
    auto ps = new ParticleSystem2(
        new PointZGrid2(16u, 16u, BBox2D::make_unit_bbox()));
    auto bh = new DomainBoundaryHandler2(
        DomainBoundaryHandler2::BoundaryType::DIRICHLET,
        DomainBoundaryHandler2::BoundaryType::NEUMANN,
        DomainBoundaryHandler2::BoundaryType::DIRICHLET,
        DomainBoundaryHandler2::BoundaryType::NEUMANN);
    PICSolverTestClass pic(bh, domain, ps);
    pic.runOnInitialize();
    bh->markFaces(domain);
    for (size_t i = 0; i < domain->cellCount(); i++)
      bh->setCellType(i, DomainBoundaryHandler2::MaterialType::FLUID);
    bh->iterateBoundaryFaces(
        [&](size_t i, DomainBoundaryHandler2::BoundaryType bt) {
          if (bt == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
            Point2d p = domain->faceCenterPosition(i);
            bh->setBoundaryFaceValue(i, cos(PI * p.x()) * cos(PI * p.y()));
          }
        });
    pic.runComputePressure(1);
    LinearVector error;
    domain->iterateCellScalar(0, [&](size_t i, double v) {
      Point2d cp = domain->cellCenterPosition(i);
      double s = cos(PI * cp.x()) * cos(PI * cp.y());
      ASSERT_NEAR(s, v, 1e-1);
      error.emplace_back(v - s);
    });
    std::cout << "ERROR " << Blas::norm(&error) / domain->cellCount()
              << std::endl;
    EXPECT_EQ(Blas::norm(&error) / domain->cellCount() < 1e-3, true);
  }
}

TEST(PICSolver2, ComputeAdvection) {
  auto region = BBox2D::make_unit_bbox();
  auto domain = new SimQuadTree(region, 4);
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

TEST(PICSolver2, SimulationStep) {
  return;
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(16);
  // dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BICUBIC);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 16.;
  for (size_t i = 0; i < 8; i++)
    for (size_t j = 0; j < 8; j++) {
      Point2d cellLower(i * h, j * h);
      // ps->addParticle(cellLower + Vector2d(h / 2, h / 2));
      // continue;
      double subh = h / 4.;
      for (int k = 1; k <= 4; k++)
        for (int l = 1; l <= 4; l++) {
          if ((i == 10 && k == 4) || (j == 10 && l == 4))
            continue;
          ps->addParticle(cellLower + Vector2d(k * subh, l * subh));
        }
    }
  pic.runOnInitialize();
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 100; i++) {
    sprintf(filename, "simulation/kim%06d", i);
    IO::saveToKim(filename, ps);

    pic.runBeginAdvanceTimeStep(dt);
    std::cout << "Velocities before pressure" << '\n';
    dm->iterateDxScalar(0, [](size_t id, double &v) {
      std::cout << id << " " << v << std::endl;
    });
    std::cout << std::endl;
    dm->iterateDyScalar(0, [](size_t id, double &v) {
      std::cout << id << " " << v << std::endl;
    });
    std::cout << '\n';
    pic.runComputeExternalForces(dt);
    pic.runApplyBoundaryCondition();
    bh->markCells(dm, ps);
    pic.runComputePressure(dt);
    pic.runCorrectVelocities(dt);
    std::cout << "Corrected velocities components" << '\n';
    dm->iterateDxScalar(0, [](size_t id, double &v) {
      std::cout << id << " " << v << std::endl;
    });
    std::cout << std::endl;
    dm->iterateDyScalar(0, [](size_t id, double &v) {
      std::cout << id << " " << v << std::endl;
    });
    std::cout << std::endl;
    pic.runTransferFromGridToParticles();
    std::cout << "Particle velocity" << '\n';
    for (size_t i = 0; i < ps->size(); i++) {
      std::cout << ps->getScalarProperty(0, i) << " "
                << ps->getScalarProperty(1, i) << '\n';
    }
    std::cout << '\n';
    std::cout << "Particle positions" << '\n';
    ps->iterateParticles(
        [&](unsigned int i, Point2d p) { std::cout << i << " " << p; });
    std::cout << std::endl;
    pic.runComputeAdvection(dt);
    std::cout << "Updated positions" << '\n';
    ps->iterateParticles(
        [&](unsigned int i, Point2d p) { std::cout << i << " " << p; });
    std::cout << std::endl;
  }
}

TEST(PICSolver2, Simulation) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
          new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolver2 pic(bh, dm, ps);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BILINEAR);
  BoxSampler sampler;
  auto h = 1. / 100.;
  for (size_t i = 45; i < 55; i++)
    for (size_t j = 5; j < 67; j++) {
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h + h / 4. + h
      // / 2., j * h + h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h - h / 4. + h
      // / 2., j * h + h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h + h / 4. + h
      // / 2., j * h - h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h - h / 4. + h
      // / 2., j * h - h / 4. + h / 2.), h / 4.)));
      ps->addParticle(
              Point2d(i * h + h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h - h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h + h / 4. + h / 2., j * h - h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h - h / 4. + h / 2., j * h - h / 4. + h / 2.));
    }
  char filename[100];
  for (int i = 0; i < 200; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename, "simulation/kim%06d", i);
    IO::saveToKim(filename, ps);
    Timer timer;
//    sprintf(filename, "simulation2/%06d", i);
//    IO::saveToPV(filename, ps, true);
//    IO::saveToCF(filename, dm, 0);
    pic.advanceFrame();
    std::cerr << timer.ellapsedTimeInSeconds() << "s" << std::endl;
  }
}

TEST(PICSolver2, DamBreaking) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolver2 pic(bh, dm, ps);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BILINEAR);

  auto h = 1. / 100.;
  for (size_t i = 0; i < 20; i++)
    for (size_t j = 0; j < 80; j++) {
      ps->addParticle(
          Point2d(i * h + h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h - h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h + h / 4. + h / 2., j * h - h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h - h / 4. + h / 2., j * h - h / 4. + h / 2.));
    }
  char filename[100];
  for (int i = 0; i < 320; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename, "simulation/kim%06d", i);
    IO::saveToKim(filename, ps);
    pic.advanceFrame();
  }
}

TEST(PICSolver2, PressureSolve) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(3);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
          new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolver2 pic(bh, dm, ps);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BILINEAR);
  BoxSampler sampler;
  std::cout << "CELL COUNT: "<< dm->grid().cellCount() << std::endl;
  std::cout << "GRID SIZE: " << dm->grid().gridSize() << std::endl;
  auto h = 1. / 3.;
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++) {
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h + h / 4. + h
      // / 2., j * h + h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h - h / 4. + h
      // / 2., j * h + h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h + h / 4. + h
      // / 2., j * h - h / 4. + h / 2.), h / 4.)));
      // ps->addParticle(sampler.sample(BBox2D(Point2d(i * h - h / 4. + h
      // / 2., j * h - h / 4. + h / 2.), h / 4.)));
      ps->addParticle(
              Point2d(i * h + h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h - h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h + h / 4. + h / 2., j * h - h / 4. + h / 2.));
      ps->addParticle(
              Point2d(i * h - h / 4. + h / 2., j * h - h / 4. + h / 2.));
    }
  char filename[100];
  for (int i = 0; i < 3; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename, "simulation/kim%06d", i);
//    IO::saveToKim(filename, ps);
    Timer timer;
//    sprintf(filename, "simulation2/%06d", i);
//    IO::saveToPV(filename, ps, true);
//    IO::saveToCF(filename, dm, 0);
    pic.advanceFrame();
    saveLinearSystem(pic.pressureMatrix,
                     pic.pressureRhs,
                     pic.pressureSolution,
                     pic.pressureDivergent,
                     i);
    std::cerr << timer.ellapsedTimeInSeconds() << "s" << std::endl;
  }
}

TEST(PICSolver2, SimulationQT) {
//  auto region = BBox2D::make_unit_bbox();
//  auto dm = new SimQuadTree(region, 3);
//  dm->setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
//                            new PointZGrid2(16));
//  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
//  auto bh =
//      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
//  PICSolverTestClass pic(bh, dm, ps);
//  auto h = 1. / 100.;
//  for (size_t i = 29; i < 42; i++)
//    for (size_t j = 5; j < 86; j++) {
//      if ((i == 34 || i == 35) && (j == 53 || j == 54))
//        continue;
//      Point2d cellLower(i * h, j * h);
//      double subh = h / 2.;
//      for (int k = 1; k <= 2; k++)
//        for (int l = 1; l <= 2; l++) {
//          // if ((i == 11 && k == 4) || (j == 11 && l == 4))
//          //  continue;
//          ps->addParticle(cellLower + Vector2d(k * subh, l * subh));
//        }
//    }
//  char filename[100];
//  double dt = 1. / 60.;
//  pic.runOnInitialize();
//  for (int i = 0; i < 20; i++) {
//    std::cerr << "Frame: " << i << '\n';
//    sprintf(filename, "simulation/kim%06d", i);
//    Timer timer;
//    pic.runBeginAdvanceTimeStep(dt);
//    pic.runComputeExternalForces(dt);
//    pic.runApplyBoundaryCondition();
//    sprintf(filename, "simulation/%06d", i);
//    IO::saveToPV(filename, ps, true);
//    IO::saveToCF(filename, dm, 0);
//    IO::saveToQT(filename, dm);
//    pic.runComputeAdvection(dt);
//    std::cerr << timer.ellapsedTimeInSeconds() << "s" << std::endl;
//  }
}

TEST(PICSolver2, NoGravityTest) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(16);
  // dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BICUBIC);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  BoxSampler sampler;
  for (unsigned i = 0; i < 10; i++)
    ps->addParticle(sampler.sample());
  ps->addParticle(Point2d(0.5, 0.));
  pic.runOnInitialize();
  double dt = 1. / 60.;
  char filename[100];
  pic.setGravity(Vector2d(0., 0.));
  for (int i = 0; i < 100; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runBeginAdvanceTimeStep(dt);
    pic.runComputeExternalForces(dt);
    pic.runApplyBoundaryCondition();
    // bh->markCells(dm, ps);
    // pic.runComputePressure(dt);
    // pic.runCorrectVelocities(dt);
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, ConstantVelocity1) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 20; i < 30; i++)
    for (size_t j = 20; j < 30; j++) {
      if (i >= 23 && i <= 27 && j >= 23 && j <= 27)
        continue;
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++)
          ps->addParticle(cellLower * 2. + Vector2d(k * subh, l * subh));
    }
  pic.runOnInitialize();
  double dt = 1. / 60.;
  char filename[100];
  for (int f = 0; f < 200; f++) {
    std::cerr << "Frame: " << f << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            f);
    IO::saveToKim(filename, ps);
    for (size_t i = 0; i < ps->size(); i++) {
      auto p = (*ps)[i];
      ps->setScalarProperty(0, i, -(p.y() - 0.5) * 0.5);
      ps->setScalarProperty(1, i, (p.x() - 0.5) * 0.5);
    }

    pic.runBeginAdvanceTimeStep(dt);
    bh->markCells(dm, ps);
    pic.runExtrapolateToDomain();
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
  for (int f = 200; f < 400; f++) {
    std::cerr << "Frame: " << f << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            f);
    IO::saveToKim(filename, ps);
    for (size_t i = 0; i < ps->size(); i++) {
      auto p = (*ps)[i];
      ps->setScalarProperty(0, i, -(p.y() - 0.5) * 0.5);
      ps->setScalarProperty(1, i, (p.x() - 0.5) * 0.5);
    }
    pic.runBeginAdvanceTimeStep(dt);
    bh->markCells(dm, ps);
    pic.runExtrapolateToDomain();
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, ConstantVelocity2) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 2; i < 4; i++)
    for (size_t j = 5; j < 10; j++) {
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++)
          ps->addParticle(cellLower + Vector2d(k * subh, l * subh));
    }
  pic.runOnInitialize();
  for (size_t i = 0; i < ps->size(); i++) {
    ps->setScalarProperty(0, i, 0.);
    ps->setScalarProperty(1, i, -0.3);
  }
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 100; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runBeginAdvanceTimeStep(dt);
    // pic.runComputeExternalForces(dt);
    // pic.runApplyBoundaryCondition();
    // bh->markCells(dm, ps);
    // pic.runComputePressure(dt);
    // pic.runCorrectVelocities(dt);
    // pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
    for (int l = 15; l >= 0; l--) {
      for (size_t c = 0; c < 15; c++)
        std::cerr << dm->scalarAtFace(0, 101 * 100 + l * 100 + c) << " ";
      std::cerr << std::endl;
    }
  }
}

TEST(PICSolver2, Zalesak) {
  /* {
     auto region = BBox2D::make_unit_bbox();
     auto dm = new SimRegularDomain2(100);
     auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
     auto bh = new DomainBoundaryHandler2(
         DomainBoundaryHandler2::BoundaryType::NEUMANN);
     PICSolverTestClass pic(bh, dm, ps);
     auto h = 1. / 100.;
     pic.runOnInitialize();
     double dt = 1. / 60.;
     char filename[100];
     for (size_t i = 0; i < 30; i++)
       for (size_t j = 0; j < 30; j++) {
         if (i >= 11 && i <= 18 && j <= 25)
           continue;
         Point2d cellLower(i * h, j * h);
         double subh = h / 1.;
         for (int k = 1; k < 2; k++)
           for (int l = 1; l < 2; l++) {
             auto p = cellLower + Vector2d(k * subh, l * subh);
             if (p.distance(Point2d(0.15, 0.15)) >= 0.15)
               continue;
             std::cout << p << std::endl;
             ps->addParticle(p + Vector2d(0.5 - 0.15, 0.75 - 0.15));
           }
       }
     for (int f = 0; f < 200; f++) {
       std::cerr << "Frame: " << f << '\n';
       sprintf(filename, "simulation/"
                         "kim%06d",
               f);
       IO::saveToKim(filename, ps);
       for (size_t i = 0; i < ps->size(); i++) {
         auto p = (*ps)[i];
         ps->setScalarProperty(0, i, (PI / 2) * (0.5 - p.y()));
         ps->setScalarProperty(1, i, (PI / 2) * (p.x() - 0.5));
       }

       pic.runBeginAdvanceTimeStep(dt);
       pic.runTransferFromGridToParticles();
       pic.runComputeAdvection(dt);
     }
     for (int f = 200; f < 400; f++) {
       std::cerr << "Frame: " << f << '\n';
       sprintf(filename, "simulation/"
                         "kim%06d",
               f);
       IO::saveToKim(filename, ps);
       for (size_t i = 0; i < ps->size(); i++) {
         auto p = (*ps)[i];
         ps->setScalarProperty(0, i, -(PI / 2) * (0.5 - p.y()));
         ps->setScalarProperty(1, i, -(PI / 2) * (p.x() - 0.5));
       }
       pic.runBeginAdvanceTimeStep(dt);
       pic.runTransferFromGridToParticles();
       pic.runComputeAdvection(dt);
     }
   }
   return;*/
  // ZALESAK's TEST
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BICUBIC);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 0; i < 30; i++)
    for (size_t j = 0; j < 30; j++) {
      if (i >= 11 && i <= 18 && j <= 25)
        continue;
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++) {
          auto p = cellLower + Vector2d(k * subh, l * subh);
          if (p.distance(Point2d(0.15, 0.15)) >= 0.15)
            continue;
          ps->addParticle(p + Vector2d(0.5 - 0.15, 0.75 - 0.15));
        }
    }

  pic.runOnInitialize();
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 400; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    dm->iterateDxScalar(
        0, [](double &v, Point2d p) { v = -(PI / 2) * (0.5 - p.y()); });
    dm->iterateDyScalar(
        0, [](double &v, Point2d p) { v = -(PI / 2) * (p.x() - 0.5); });
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    bh->markCells(dm, ps);
    pic.runExtrapolateToDomain();
    pic.runComputeAdvection(dt);
  }
  for (int i = 400; i < 800; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    dm->iterateDxScalar(
        0, [](double &v, Point2d p) { v = (PI / 2) * (0.5 - p.y()); });
    dm->iterateDyScalar(
        0, [](double &v, Point2d p) { v = (PI / 2) * (p.x() - 0.5); });
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    bh->markCells(dm, ps);
    pic.runExtrapolateToDomain();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, ConstantGridVelocity1) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 2; i < 40; i++)
    for (size_t j = 5; j < 10; j++) {
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++)
          ps->addParticle(cellLower * 2. + Vector2d(k * subh, l * subh));
    }
  pic.runOnInitialize();
  dm->iterateDxScalar(0, [](double &v, Point2d p) {
    v = 0.3;
    UNUSED_VARIABLE(p);
  });
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 100; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, ConstantGridVelocity2) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 20; i < 30; i++)
    for (size_t j = 20; j < 30; j++) {
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++)
          ps->addParticle(cellLower * 2. + Vector2d(k * subh, l * subh));
    }
  pic.runOnInitialize();
  dm->iterateDxScalar(0, [](double &v, Point2d p) { v = -(p.y() - 0.5) * 1.; });
  dm->iterateDyScalar(0, [](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = (p.x() - 0.5) * 1.;
  });
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 200; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
  dm->iterateDxScalar(0, [](double &v, Point2d p) { v = (p.y() - 0.5) * 1.; });
  dm->iterateDyScalar(0, [](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = -(p.x() - 0.5) * 1.;
  });
  for (int i = 200; i < 400; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, ConstantGridVelocity3) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(64u, 64u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverTestClass pic(bh, dm, ps);
  auto h = 1. / 100.;
  for (size_t i = 2; i < 40; i++)
    for (size_t j = 5; j < 10; j++) {
      Point2d cellLower(i * h, j * h);
      double subh = h / 1.;
      for (int k = 1; k < 2; k++)
        for (int l = 1; l < 2; l++)
          ps->addParticle(cellLower * 2. + Vector2d(k * subh, l * subh));
    }
  pic.runOnInitialize();
  dm->iterateDxScalar(0, [](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = 0.3;
  });
  dm->iterateDyScalar(0, [](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v = 0.3;
  });
  double dt = 1. / 60.;
  char filename[100];
  for (int i = 0; i < 200; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename,
            "simulation/"
                "kim%06d",
            i);
    IO::saveToKim(filename, ps);
    pic.runTransferFromGridToParticles();
    pic.runComputeAdvection(dt);
  }
}

TEST(PICSolver2, WallSplash) {
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolver2 pic(bh, dm, ps);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BICUBIC);
  auto h = 1. / 100.;
  for (size_t i = 29; i < 42; i++)
    for (size_t j = 40; j < 86; j++) {
      Point2d cellLower(i * h, j * h);
      double subh = h / 2.;
      for (int k = 1; k <= 2; k++)
        for (int l = 1; l <= 2; l++) {
          // if ((i == 11 && k == 4) || (j == 11 && l == 4))
          //  continue;
          ps->addParticle(cellLower + Vector2d(k * subh, l * subh));
        }
    }
  char filename[100];
  pic.setGravity(Vector2d(-9.8, 0.));
  for (int i = 0; i < 400; i++) {
    std::cerr << "Frame: " << i << '\n';
    sprintf(filename, "simulation/kim%06d", i);
    IO::saveToKim(filename, ps);
    Timer timer;
    pic.advanceFrame();
    if (i == 0) {
      for (size_t j = 0; j < ps->size(); j++) {
        ps->setScalarProperty(0, j, -4.5);
        ps->setScalarProperty(1, j, 0.);
      }
    }
    std::cerr << timer.ellapsedTimeInSeconds() << "s" << std::endl;
  }
}
