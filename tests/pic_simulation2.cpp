#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

Point2d advectExact(Point2d p,
                    double dt,
                    std::function<Vector2d(Point2d)> velocity) {
  unsigned int numSubSteps = 2;
  double subdt = dt / numSubSteps;
  Point2d pt = p;
  for (unsigned int t = 0; t < numSubSteps; t++) {
    Vector2d v0 = velocity(pt);
    // Mid-point rule
    Point2d midPt = pt + 0.5 * subdt * v0;
    Vector2d midVel = velocity(midPt);
    pt += subdt * midVel;
  }
  return pt;
}

Vector2d enright(Point2d p) {
  return {
      2. * SQR(sin(PI * p.x())) * sin(PI * p.y())
          * cos(PI * p.y()),
      -2. * SQR(sin(PI * p.y())) * sin(PI * p.x())
          * cos(PI * p.x())};
}

TEST(AdaptiveCenteredPic2, ParticleGridTransfer) {
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 5, 6);
  auto particles = simulation.particles();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &velocityXField = *structure->field<double>(8);
  auto &velocityYField = *structure->field<double>(9);
  auto &cellMaterialField = *structure->field<Definitions::Material>(2);
  Injector2 injector(particles, 0.013);
  injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  simulation.buildTree();
  for (size_t i = 0; i < 50; i++) {
    simulation.markCells();
    simulation.markFaces();
    // set particles velocities from analytical function
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setScalarProperty(0, id, enright(p).x());
      particles->setScalarProperty(1, id, enright(p).y());
    });
    // transfer particles -> grid
    simulation.transferFromParticlesToGrid();
    // compute grid error
    double maxError = 0.;
    structure->iterateCells([&](size_t id) {
      if (cellMaterialField[id] == Definitions::Material::FLUID) {
        auto v = enright(structure->cellCenterPosition(id));
        auto vx = velocityXField[id];
        auto vy = velocityYField[id];
        maxError = std::max(maxError, std::fabs(v.x() - vx));
        maxError = std::max(maxError, std::fabs(v.y() - vy));
      }
    });
    std::cerr << maxError << std::endl;
    EXPECT_LT(maxError, 1e-2);
    // advect particles
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setPosition(id, advectExact(p, 0.001, enright));
    });
    simulation.updateTree();
  }

}

TEST(AdaptiveCenteredPic2, GridParticleTransfer) {
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 5, 6);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &velocityField = *structure->field<double>(2);
  auto &faceMaterialField = *structure->field<Definitions::Material>(1);
  Injector2 injector(particles, 0.013);
  injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  simulation.buildTree();
  for (size_t i = 0; i < 50; i++) {
    simulation.markCells();
    simulation.markFaces();
    // set grid velocities
    structure->iterateFaces([&](size_t id) {
      if (faceMaterialField[id] == Definitions::Material::FLUID) {
        double v = (structure->faceOrientation(id)
            == Definitions::Orientation::HORIZONTAL) ?
                   enright(structure->faceCenterPosition(id)).y() :
                   enright(structure->faceCenterPosition(id)).x();
        velocityField[id] = v;
      }
    });
    // transfer grid -> particles
    simulation.transferFromGridToParticles();
    // compute grid error
    Vector2d maxError;
    particles->iterateParticles([&](size_t id, Point2d p) {
      Vector2d v(particles->getScalarProperty(0, id),
                 particles->getScalarProperty(1, id));
      auto av = enright(p);
      maxError = Vector2d(std::max(maxError.x(), std::fabs(v.x() - av.x())),
                          std::max(maxError.y(), std::fabs(v.y() - av.y())));
      particles->setPosition(id, advectExact(p, 0.001, enright));
    });
    std::cerr << maxError << std::endl;
    // advect particles
    simulation.updateTree();
  }

}

TEST(AdaptivePic2, ParticleGridTransfer) {
  auto domainRegion = BBox2d::squareBox();
  AdaptivePic2 simulation(domainRegion, 5, 6);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &velocityField = *structure->field<double>(2);
  auto &faceMaterialField = *structure->field<Definitions::Material>(1);
  Injector2 injector(particles, 0.013);
  injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  simulation.buildTree();
  for (size_t i = 0; i < 50; i++) {
    simulation.markCells();
    simulation.markFaces();
    // set particles velocities from analytical function
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setScalarProperty(0, id, enright(p).x());
      particles->setScalarProperty(1, id, enright(p).y());
    });
    // transfer particles -> grid
    simulation.transferFromParticlesToGrid();
    // compute grid error
    double maxError = 0.;
    structure->iterateFaces([&](size_t id) {
      if (faceMaterialField[id] == Definitions::Material::FLUID) {
        double v = (structure->faceOrientation(id)
            == Definitions::Orientation::HORIZONTAL) ?
                   enright(structure->faceCenterPosition(id)).y() :
                   enright(structure->faceCenterPosition(id)).x();
        maxError = std::max(maxError, std::fabs(v - velocityField[id]));
      }
    });
    std::cerr << maxError << std::endl;
    // advect particles
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setPosition(id, advectExact(p, 0.001, enright));
    });
    simulation.updateTree();
  }

}

TEST(AdaptivePic2, GridParticleTransfer) {
  auto domainRegion = BBox2d::squareBox();
  AdaptivePic2 simulation(domainRegion, 5, 6);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &velocityField = *structure->field<double>(2);
  auto &faceMaterialField = *structure->field<Definitions::Material>(1);
  Injector2 injector(particles, 0.013);
  injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  simulation.buildTree();
  for (size_t i = 0; i < 50; i++) {
    simulation.markCells();
    simulation.markFaces();
    // set grid velocities
    structure->iterateFaces([&](size_t id) {
      if (faceMaterialField[id] == Definitions::Material::FLUID) {
        double v = (structure->faceOrientation(id)
            == Definitions::Orientation::HORIZONTAL) ?
                   enright(structure->faceCenterPosition(id)).y() :
                   enright(structure->faceCenterPosition(id)).x();
        velocityField[id] = v;
      }
    });
    // transfer grid -> particles
    simulation.transferFromGridToParticles();
    // compute grid error
    Vector2d maxError;
    particles->iterateParticles([&](size_t id, Point2d p) {
      Vector2d v(particles->getScalarProperty(0, id),
                 particles->getScalarProperty(1, id));
      auto av = enright(p);
      maxError = Vector2d(std::max(maxError.x(), std::fabs(v.x() - av.x())),
                          std::max(maxError.y(), std::fabs(v.y() - av.y())));
      particles->setPosition(id, advectExact(p, 0.001, enright));
    });
    std::cerr << maxError << std::endl;
    // advect particles
    simulation.updateTree();
  }

}

TEST(AdaptivePic2, SolvePressure) {
  auto domainRegion = BBox2d::squareBox();
  AdaptivePic2 simulation(domainRegion, 5, 6);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  Injector2 injector(particles, 0.013);
  injector.setupBoxShape(BBox2d(Point2d(0.3, 0.01), Point2d(0.6, 0.3)));
  // set particles velocities from analytical function
  simulation.buildTree();
  simulation.markCells();
  simulation.markFaces();
  simulation.computeExternalForces(0.001);
  simulation.applyBoundaryCondition();
  simulation.solvePressure(0.001);
  simulation.updateTree();
}

TEST(AdaptiveCenteredPic2, CornerDivergence) {
  DifferentialRBF<2> rbf;
  std::vector<Point2d> points;
  std::vector<double> values;

  points.emplace_back(0.5, 0.5);
  points.emplace_back(0.4, .5);
  points.emplace_back(0.55, 0.5);
  points.emplace_back(0.5, 0.6);
  points.emplace_back(0.5, 0.45);

  values.emplace_back(-9.81);
  values.emplace_back(-9.81);
  values.emplace_back(-9.81);
  values.emplace_back(-9.81);
  values.emplace_back(0);

  std::cerr << rbf.gradientAt(Point2d(0.5, 0.5), 1, points, values)
            << std::endl;

  points[2] = Point2d(0.6, 0.5);
  std::cerr << rbf.gradientAt(Point2d(0.5, 0.5), 1, points, values)
            << std::endl;

  points[2] = Point2d(0.55, 0.5);
  values[2] = 0;

  std::cerr << rbf.gradientAt(Point2d(0.5, 0.5), 1, points, values)
            << std::endl;
}