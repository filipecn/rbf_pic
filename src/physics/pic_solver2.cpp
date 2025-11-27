#include "physics/pic_solver2.h"
#include <fstream>
#include <furoo.h>
#include <geometry/transform.h>
#include <map>
#include <sstream>

namespace furoo {

PICSolver2::PICSolver2() {
  _temporalIntegrator.reset(new EulerTemporalIntegrator2());
}

PICSolver2::PICSolver2(DomainBoundaryHandler2 *b, SimulationDomain2 *d,
                       ParticleSystem2 *p, EulerTemporalIntegrator2 *e) {
  _domain.reset(d);
  _particles = p;
  _boundaryHandler = b;
  _domain->setBoundaryHandler(b);
  _temporalIntegrator.reset(e);
  _pressureField = _domain->addCellCenteredScalarField(0.);
  _velocityField = _domain->addFaceCenteredScalarField(0.);
  _particleVx = _particles->addScalarProperty();
  _particleVy = _particles->addScalarProperty();
}

PICSolver2::~PICSolver2() {}

void PICSolver2::setSimulationDomain(SimulationDomain2 *d) {
  _domain.reset(d);
  _pressureField = _domain->addCellCenteredScalarField(0.);
  _velocityField = _domain->addFaceCenteredScalarField(0.);
}

void PICSolver2::setBoundaryHandler(DomainBoundaryHandler2 *bh) {
  _boundaryHandler = bh;
  if (_domain)
    _domain->setBoundaryHandler(_boundaryHandler);
  else
    throw("Domain has not been set yet!");
}

DomainBoundaryHandler2 *PICSolver2::boundaryHandler() {
  return _boundaryHandler;
}

void PICSolver2::setParticleSystem(ParticleSystem2 *ps) {
  _particles = ps;
  _particleVx = _particles->addScalarProperty();
  _particleVy = _particles->addScalarProperty();
}

SimulationDomain2 *PICSolver2::domain() { return _domain.get(); }

ParticleSystem2 *PICSolver2::particles() { return _particles; }

void PICSolver2::setGravity(Vector2d g) { _gravity = g; }

void PICSolver2::onInitialize() {
  _boundaryHandler->markCells(_domain.get(), _particles);
  _boundaryHandler->markFaces(_domain.get());
}

void PICSolver2::onAdvanceTimeStep(double dt) {
  ASSERT_FATAL(_domain);
  ASSERT_FATAL(_particles);
  unsigned int subdt = numberOfSubTimeSteps(dt);
  for (size_t i = 0; i < subdt; i++) {
    beginAdvanceTimeStep(dt / subdt);
    computeExternalForces(dt / subdt);
    applyBoundaryCondition();
    computeViscosity(dt / subdt);
    _boundaryHandler->markCells(_domain.get(), _particles);
    _boundaryHandler->markFaces(_domain.get());
    propagateVelocitiesToDomain();
    std::cerr << "computePressure .. ";
    computePressure(dt / subdt);
    std::cerr << " finished.\n";
    correctVelocities(dt / subdt);
    propagateVelocitiesToDomain();
    applyBoundaryCondition();
    transferFromGridToParticles();
    computeAdvection(dt / subdt);
    endAdvanceTimeStep(dt / subdt);
  }
}

double PICSolver2::cfl(double dt) const {
  double maxVel = 0.;
  _domain->iterateDxScalar(_velocityField, [&](size_t i, double &v) {
    UNUSED_VARIABLE(i);
    maxVel = max(maxVel, std::fabs(v + dt * _gravity.x()));
  });
  _domain->iterateDyScalar(_velocityField, [&](size_t i, double &v) {
    UNUSED_VARIABLE(i);
    maxVel = max(maxVel, std::fabs(v + dt * _gravity.y()));
  });
  return (maxVel * dt) / std::min(_domain->smallestCellSize()[0],
                                  _domain->smallestCellSize()[1]);
}

unsigned int PICSolver2::numberOfSubTimeSteps(double dt) const {
  double currentCfl = cfl(dt);
  return static_cast<unsigned int>(std::max(std::ceil(currentCfl * 1.1), 1.0));
}

void PICSolver2::applyBoundaryCondition() {
  _boundaryHandler->iterateBoundaryFaces(
      [&](size_t i, DomainBoundaryHandler2::BoundaryType t) {
        if (t == DomainBoundaryHandler2::BoundaryType::NEUMANN)
          _domain->faceCenteredScalar(_velocityField, i) = 0.;
      });
}

void PICSolver2::beginAdvanceTimeStep(double dt) {
  UNUSED_VARIABLE(dt);
  transferFromParticlesToGrid();
}

void PICSolver2::computeExternalForces(double dt) {
  _domain->iterateDxScalar(_velocityField, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v += dt * _gravity.x();
  });
  _domain->iterateDyScalar(_velocityField, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    v += dt * _gravity.y();
  });
}

void PICSolver2::computeViscosity(double dt) {
  UNUSED_VARIABLE(dt);
  NOT_IMPLEMENTED();
}

void PICSolver2::computePressure(double dt) {
  UNUSED_VARIABLE(dt);
  _domain->iterateCellScalar(_pressureField, [](double &v, Point2d id) {
    UNUSED_VARIABLE(id);
    v = 0.;
  });
  std::map<unsigned, unsigned> indexMap;
  unsigned fluidCellCount = 0;
  for (size_t i = 0; i < _domain->cellCount(); i++)
    if (_boundaryHandler->cellType(i) ==
        DomainBoundaryHandler2::MaterialType::FLUID)
      indexMap[i] = fluidCellCount++;
  LinearVector b(fluidCellCount + 1, 0.);
  LinearVector x(fluidCellCount + 1, 0.);
  LinearMatrix A(fluidCellCount + 1);
  // double dx = _domain->smallestCellSize()[0];
  for (auto c : indexMap) {
    _domain->iterateStencilAtCellCenter(
        c.first, [&](int id, double w, DomainBoundaryHandler2::BoundaryType bt,
                     double bv) {
          UNUSED_VARIABLE(bt);
          UNUSED_VARIABLE(bv);
          if (id >= 0) {
            auto it = indexMap.find(id);
            if (it != indexMap.end())
              A(c.second, indexMap[id]) = w;
          }
        });
    double nsum = 0.;
    double dx = _domain->cellRegion(c.first).size(0);
    auto neighbors = _domain->cellFaces(c.first);
    for (auto n : neighbors) {
      if (_boundaryHandler->faceBoundaryType(n.id) !=
          DomainBoundaryHandler2::BoundaryType::NEUMANN)
        continue;
      switch (n.orientation) {
      case SimulationDomain2::Orientation::BOTTOM:
        nsum -= _domain->faceCenteredScalar(_velocityField, n.id);
        break;
      case SimulationDomain2::Orientation::LEFT:
        nsum -= _domain->faceCenteredScalar(_velocityField, n.id);
        break;
      case SimulationDomain2::Orientation::TOP:
        nsum += _domain->faceCenteredScalar(_velocityField, n.id);
        break;
      case SimulationDomain2::Orientation::RIGHT:
        nsum += _domain->faceCenteredScalar(_velocityField, n.id);
        break;
      default:
        ASSERT_FATAL(false);
      }
    }
    // double d = -_domain->faceDivergentAtCellCenter(_velocityField, c.first);
    // std::cerr << "divergence at cell " << c.first << " " << d << std::endl;
    b[c.second] =
        -_domain->faceDivergentAtCellCenter(_velocityField, c.first) / dt +
        nsum / (dt * dx);
  }
  //  std::cout << dt << std::endl;
  //  A(fluidCellCount, fluidCellCount) = 1.;
  //  this->pressureMatrix = A;
  //  b[fluidCellCount] = 0.;
  //  this->pressureRhs = b;
  //  this->pressureDivergent = div;
  ConjugateGradientLinearSolver bicg;
  bicg.solve(&A, &b, &x);
  //  this->pressureSolution = x;

  for (auto c : indexMap)
    _domain->scalarAtCell(_pressureField, c.first) = x[c.second];
}

void PICSolver2::correctVelocities(double dt) {
  UNUSED_VARIABLE(dt);
  _domain->iterateDxScalar(_velocityField, [&](size_t id, double &v) {
    if (_boundaryHandler->faceType(id) ==
        DomainBoundaryHandler2::MaterialType::FLUID)
      v -= _domain->cellGradientAtFaceCenter(_pressureField, id) * dt;
  });
  _domain->iterateDyScalar(_velocityField, [&](size_t id, double &v) {
    if (_boundaryHandler->faceType(id) ==
        DomainBoundaryHandler2::MaterialType::FLUID)
      v -= _domain->cellGradientAtFaceCenter(_pressureField, id) * dt;
  });
}

void PICSolver2::propagateVelocitiesToDomain() {
  _domain->extrapolateToDomain(_velocityField);
}

void PICSolver2::computeAdvection(double dt) {
  _particles->iterateParticles([&](unsigned int i, Point2d p) {
    // Runge kutta 2-order
    // unsigned int numSubSteps =
    //     static_cast<unsigned int>(std::max(_maxCfl, 1.0));
    // numSubSteps = 2;
    // double subdt = dt / numSubSteps;
    // Point2d pt = p;
    // try {
    // for (unsigned int t = 0; t < numSubSteps; t++) {
    //   Vector2d v0 = _domain->sampleFaceCenteredScalar(_velocityField, pt);
    //   // Mid-point rule
    //   Point2d midPt = pt + 0.5 * subdt * v0;
    //   Vector2d midVel =
    //       _domain->sampleFaceCenteredScalar(_velocityField, midPt);
    //   pt += subdt * midVel;
    // }
    // } catch (const char *e) {
    //   std::cerr << "Thrown from " << e << '\n';
    //   throw("computeAdvection");
    // }
    // auto np = pt;

    // EULER
    auto v = Vector2d(_particles->getScalarProperty(_particleVx, i),
                      _particles->getScalarProperty(_particleVy, i));
    auto np = _temporalIntegrator->integrate(p, v, dt);
    // Vector2d v = _domain->sampleFaceCenteredScalar(_velocityField, p);

    // Check if point is still inside domain
    if (p != np) {
      auto in = _boundaryHandler->intersect(_domain.get(), p, np);
      if (in.isValid) {
        // auto r =
        //    Transform2D::applyReflection(np - in.point,
        //    in.normal.normalized());
        // np = in.point + r;
        // auto nv = (in.normal.x() != 0.) ? Vector2d(-v.x(), v.y())
        //                                : Vector2d(v.x(), -v.y());
        if (in.normal.x() != 0.) {
          _particles->setScalarProperty(_particleVx, i, 0.);
          np = Point2d(in.point.x(), np.y());
        }
        if (in.normal.y() != 0.) {
          _particles->setScalarProperty(_particleVy, i, 0.);
          np = Point2d(np.x(), in.point.y());
        }
        if (!_domain->region().inside(np))
          np = in.point; // + 0.001 * (p - in.point);
      }
    }
    np = Point2d(clamp(np.x(), 0.0, 1.0), clamp(np.y(), 0.0, 1.0));
    _particles->setPosition(i, np);
    if (!_domain->region().inside(np)) {
      std::cerr << "Point: " << p << '\n';
      std::cerr << "NewPoint: " << np << '\n';
      throw("computeAdvection: new point outside domain");
    }
  });
}

void PICSolver2::endAdvanceTimeStep(double dt) {
  UNUSED_VARIABLE(dt);
  NOT_IMPLEMENTED();
}

void PICSolver2::transferFromParticlesToGrid() {
  // double radius = _domain->smallestCellSize()[0] * 2.5;
  /* WAAAAAAAAT???
  _domain->iterateDxScalar(_velocityField, [&](size_t i, double &v) {
    auto neighborhood = _domain->faceNeighborhood(i);
    auto p = _domain->faceCenterPosition(i);
    if (_boundaryHandler->faceBoundaryType(i) ==
  DomainBoundaryHandler2::BoundaryType::DIRICHLET) { v =
  _particles->sampleScalarProperty(_particleVx, radius, p, Vector2d(0., 1.)); }
  else if ((neighborhood[0].id >= 0
        && _boundaryHandler->cellType(neighborhood[0].id) ==
  DomainBoundaryHandler2::MaterialType::FLUID) && (neighborhood[1].id >= 0
            && _boundaryHandler->cellType(neighborhood[1].id) ==
  DomainBoundaryHandler2::MaterialType::FLUID)) { v =
  _particles->sampleScalarProperty(_particleVx, radius, p);
    }
  });
  _domain->iterateDyScalar(_velocityField, [&](size_t i, double &v) {
    auto neighborhood = _domain->faceNeighborhood(i);
    auto p = _domain->faceCenterPosition(i);
    if (_boundaryHandler->faceBoundaryType(i) ==
  DomainBoundaryHandler2::BoundaryType::DIRICHLET) { v =
  _particles->sampleScalarProperty(_particleVy, radius, p, Vector2d(1., 0.)); }
  else if ((neighborhood[0].id >= 0
        && _boundaryHandler->cellType(neighborhood[0].id) ==
  DomainBoundaryHandler2::MaterialType::FLUID) && (neighborhood[1].id >= 0
            && _boundaryHandler->cellType(neighborhood[1].id) ==
  DomainBoundaryHandler2::MaterialType::FLUID)) { v =
  _particles->sampleScalarProperty(_particleVy, radius, p);
    }
  });*/
  /*_domain->iterateDxScalar(_velocityField, [&](double &v, Point2d p) {
    v = _particles->sampleScalarProperty(_particleVx, radius, p);
  });
  _domain->iterateDyScalar(_velocityField, [&](double &v, Point2d p) {
    v = _particles->sampleScalarProperty(_particleVy, radius, p);
  });*/

  _domain->iterateDxScalar(_velocityField, [&](size_t i, double &v) {
    auto p = _domain->faceCenterPosition(i);
    if (_boundaryHandler->faceType(i) ==
        DomainBoundaryHandler2::MaterialType::FLUID)
      v = _particles->sampleScalarProperty(_particleVx, p, 13u);
  });
  _domain->iterateDyScalar(_velocityField, [&](size_t i, double &v) {
    auto p = _domain->faceCenterPosition(i);
    if (_boundaryHandler->faceType(i) ==
        DomainBoundaryHandler2::MaterialType::FLUID)
      v = _particles->sampleScalarProperty(_particleVy, p, 13u);
  });
}

void PICSolver2::transferFromGridToParticles() {
  _particles->iterateParticles([&](unsigned int i, Point2d p) {
    auto v = _domain->sampleFaceCenteredScalar(_velocityField, p);
    _particles->setScalarProperty(_particleVx, i, v.x());
    _particles->setScalarProperty(_particleVy, i, v.y());
  });
}

void PICSolver2::exportVelocityField(int frame) const {
  std::ostringstream filename;
  filename << "velocity." << std::to_string(frame) << ".csv";
  std::ofstream outFile(filename.str() + ".temp");

  _domain->iterateDxScalar(_velocityField, [&](size_t i, double &v) {
    UNUSED_VARIABLE(i);
    //    std::cout << "DxScalar = " << i << std::endl;
    outFile << v << "," << std::endl;
  });
  outFile.close();

  std::ifstream inFile(filename.str() + ".temp");
  std::ofstream finalFile(filename.str());

  std::string str;
  _domain->iterateDyScalar(_velocityField, [&](size_t i, double &v) {
    UNUSED_VARIABLE(i);
    //    std::cout << "DyScalar = " << i << std::endl;
    std::getline(inFile, str);
    finalFile << str << v << std::endl;
  });

  inFile.close();
  finalFile.close();
}

} // namespace furoo
