// #include "ponos_wrapper.h"
#include "test_utils.h"
#include <furoo.h>
#include <physics/particle_system2.h>
#include <physics/simulation_regular_domain.h>

using namespace furoo;

#define WIDTH 1000
#define HEIGHT 1000

Vector2d velocity(Point2d p) {
  return Vector2d(2. * SQR(sin(PI * p.x())) * sin(PI * p.y()) * cos(PI * p.y()),
                  -2. * SQR(sin(PI * p.y())) * sin(PI * p.x()) *
                      cos(PI * p.x()));
  // return Vector2d(-(PI / 2.) * (0.5 - p.y()), -(PI / 2.) * (p.x() - 0.5));
}

int main(int argc, char const *argv[]) {

  double injectorSpacing = 0.01;
  size_t quadtreeLevel = 4;
  char folderName[100];
  sprintf(folderName, "simulation/");
  if (argc == 2) {
    sscanf(argv[1], "%zu", &quadtreeLevel);
  }
  if (argc == 3) {
    sscanf(argv[2], "%lf", &injectorSpacing);
  }
  if (argc == 4) {
    sprintf(folderName, "%s", argv[3]);
  }

  ParticleSystem2 *numericalParticles;
  std::shared_ptr<DomainBoundaryHandler2> bh;
  std::shared_ptr<SimQuadTree> domain;
  // exact
  ParticleSystem2 *analyticalParticles;

  double dt = 1 / 200.;
  auto region = BBox2D::make_unit_bbox();

  numericalParticles = new ParticleSystem2(
      new PointZGrid2(32u, 32u, region),
      new RBFInterpolant<Point2d>(new QuinticKernel<Point2d>()));
  Injector injector(numericalParticles, injectorSpacing);
  injector.setupCircleShape(Point2d(0.5, 0.75), 0.15);
  std::cerr << "injectorSpacing: " << injectorSpacing << '\n';
  std::cerr << "Psize: " << numericalParticles->size() << '\n';
  analyticalParticles = new ParticleSystem2(
      new PointZGrid2(32u, 32u, region),
      new RBFInterpolant<Point2d>(new QuinticKernel<Point2d>()));
  analyticalParticles->addScalarProperty();
  analyticalParticles->addScalarProperty();
  numericalParticles->iterateParticles([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(i);
    analyticalParticles->addParticle(p);
  });

  // Setting up domain
  HaltonSequence rngX(3), rngY(5);
  std::vector<Point2d> samples;
  for (int i = 0; i < 47; i++)
    samples.emplace_back(rngX.randomDouble(), rngY.randomDouble());

  domain.reset(new SimQuadTree(
      BBox2D::make_unit_bbox(), [&](QuadTree::Node &node) -> bool {
        if (node.level() > quadtreeLevel)
          return false;
        return true;
        for (auto &p : samples)
          if (node.region().inside(p))
            return true;
        return false;
        int count = 0;
        numericalParticles->iterateParticles(node.region(),
                                             [&](unsigned int i, Point2d p) {
                                               UNUSED_VARIABLE(p);
                                               UNUSED_VARIABLE(i);
                                               count++;
                                             });
        return count > 0 && node.level() < 7;
      }));
  domain->setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                                new PointZGrid2(16));
  bh.reset(new DomainBoundaryHandler2(
      DomainBoundaryHandler2::BoundaryType::NEUMANN));

  PICSolverDebuggerClass pic(bh.get(), domain.get(), numericalParticles);

  // Simulation steps
  for (int curFrame = 0; curFrame < 800; curFrame++) {
    std::cerr << "Frame: " << curFrame << std::endl;

    // advect analytical particles
    analyticalParticles->iterateParticles([&](int i, Point2d p) {
      double tdt = dt;
      if (curFrame >= 400)
        tdt = -dt;
      auto np = advectExact(p, tdt, velocity);
      analyticalParticles->setPosition(i, np);
    });
    // mark grid with numerical particles
    bh->markCells(domain.get(), numericalParticles);
    bh->markFaces(domain.get());

    // update analytical particles velocities
    for (size_t i = 0; i < numericalParticles->size(); i++) {
      auto p = (*numericalParticles)[i];
      numericalParticles->setScalarProperty(
          0, i, velocity(p).x() * ((curFrame >= 400) ? -1. : 1.));
      numericalParticles->setScalarProperty(
          1, i, velocity(p).y() * ((curFrame >= 400) ? -1. : 1.));
    }

    pic.runTransferFromParticlesToGrid();
    pic.runComputeAdvection(dt);

    for (size_t i = 0; i < numericalParticles->size(); i++) {
      auto p = (*numericalParticles)[i];
      numericalParticles->setScalarProperty(0, i, velocity(p).x());
      numericalParticles->setScalarProperty(1, i, velocity(p).y());
    }
    bh->markCells(domain.get(), numericalParticles);
    bh->markFaces(domain.get());
    pic.runTransferFromParticlesToGrid();

    char filename[100];
    sprintf(filename, "%s/particles_%06d", folderName, curFrame);
    IO::saveToPV(filename, numericalParticles, true);
    sprintf(filename, "%s/particles_exact%06d", folderName, curFrame);
    IO::saveToPV(filename, analyticalParticles, true);
  }

  delete numericalParticles;
  delete analyticalParticles;

  return 0;
}
