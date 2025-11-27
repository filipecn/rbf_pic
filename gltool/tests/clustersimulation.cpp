// #include "ponos_wrapper.h"
#include "test_utils.h"
#include <furoo.h>
#include <physics/particle_system2.h>
#include <physics/simulation_regular_domain.h>

using namespace furoo;

#define WIDTH 1000
#define HEIGHT 1000

int main(int argc, char const *argv[]) {
  // Configuration
  if (argc < 6) {
    std::cout << "Usage: <quadtree level> <particle spacing> <# of frames> "
                 "<dt> <folder name>\n";
    return 0;
  }
  size_t quadtreeLevel = 0;
  double particleSpacing = 0.;
  int nframes = 0;
  double dt = 0.;
  char folderName[100];
  sscanf(argv[1], "%zu", &quadtreeLevel);
  sscanf(argv[2], "%lf", &particleSpacing);
  sscanf(argv[3], "%d", &nframes);
  sscanf(argv[4], "%lf", &dt);
  sprintf(folderName, "%s", argv[5]);

  // Initialization
  BBox2D region = BBox2D::make_unit_bbox();
  std::shared_ptr<SimQuadTree> domain(
      new SimQuadTree(region, [quadtreeLevel](QuadTree::Node &node) -> bool {
        return node.level() < quadtreeLevel;
      }));
  domain->setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                                new PointZGrid2(16));
  std::shared_ptr<DomainBoundaryHandler2> boundaryHandler(
      new DomainBoundaryHandler2(
          DomainBoundaryHandler2::BoundaryType::NEUMANN));
  std::shared_ptr<ParticleSystem2> particles(new ParticleSystem2(
      new PointZGrid2(32u, 32u, region), new ShepardInterpolant2()));
  // new RBFInterpolant<Point2d>(new QuinticKernel<Point2d>())));
  domain->setInterpolationMethod(
      SimulationDomain2::InterpolationMethod::BILINEAR);

  // Inject Particles
  Injector injector(particles.get(), particleSpacing);
  // DROP PARTICLES'
  // injector.setupCircleShape(Point2d(0.5, 0.8), 0.08);
  // injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.999, 0.2));
  // DAM PARTICLES
  injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.3, 0.999));

  // PIC
  PICSolverDebuggerClass pic(boundaryHandler.get(), domain.get(),
                             particles.get());

  // Simulation steps
  for (int curFrame = 0; curFrame < nframes; curFrame++) {
    std::cerr << "Frame: " << curFrame << std::endl;

    unsigned int subdt = pic.getNumberOfSubTimeSteps(dt);
    for (size_t i = 0; i < subdt; i++) {
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runTransferFromParticlesToGrid();
      pic.runComputeExternalForces(dt);
      pic.runApplyBoundaryCondition();
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runExtrapolateToDomain();
      pic.runComputePressure(dt);
      pic.runCorrectVelocities(dt);
      pic.runExtrapolateToDomain();
      pic.runApplyBoundaryCondition();
      pic.runTransferFromGridToParticles();
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runComputeAdvection(dt);
    }

    char filename[100];
    sprintf(filename, "%s/particles_%06d", folderName, curFrame);
    pic.saveParticleSystem(filename);
  }

  return 0;
}
