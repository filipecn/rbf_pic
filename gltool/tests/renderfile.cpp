#include "SearchRenderer.h"
#include "boundary_handler2_model.h"
#include "particle_system2_model.h"
#include "simulation_quad_tree_domain_model.h"
#include "test_utils.h"
#include <aergia/aergia.h>
#include <furoo.h>
#include <ponos_wrapper.h>

using namespace furoo;

#define WIDTH 800
#define HEIGHT 800

int main(int argc, char **argv) {
  // Configuration
  if (argc < 4) {
    std::cout << "Usage: <quadtree level> <folder name> <# of frames>\n";
    return 0;
  }
  size_t numberOfFrames;
  size_t quadtreeLevel;
  char filename[100];

  sscanf(argv[1], "%zu", &quadtreeLevel);
  sscanf(argv[3], "%zu", &numberOfFrames);
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
  Injector injector(particles.get(), 0.002);
  // DAM BREAK
  injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.3, 0.999));
  // DROP
  // injector.setupCircleShape(Point2d(0.5, 0.8), 0.08);
  // injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.999, 0.2));
  // PIC
  PICSolverDebuggerClass pic(boundaryHandler.get(), domain.get(),
                             particles.get());
  // VISUALIZATION
  aergia::SceneApp<> app(WIDTH, HEIGHT, "PIC 2D", false);
  app.init();
  app.addViewport2D(0, 0, WIDTH, HEIGHT);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(region), 1.1f);
  // HELPERS
  ParticleSystem2Model particlesModel(particles.get(), 0.003);
  SimulationQuadTreeDomainModel domainModel(domain.get(),
                                            boundaryHandler.get());
  app.scene.add(&particlesModel);
  app.scene.add(&domainModel);

  // sprintf(filename, "./%s/particles_000000", argv[2]);
  // printf("%s\n", filename);
  // injector.loadFromFile(filename);
  // std::cout << "loaded; " << particles->size() << '\n';

  app.renderCallback = [&]() {
    static int curFrame = 0;
    if (curFrame > static_cast<int>(numberOfFrames))
      return;
    // return;
    // pic.runAdvanceTimeStep(dt);
    sprintf(filename, "./%s/particles_%06d", argv[2], curFrame);
    printf("%s\n", filename);
    injector.loadFromFile(filename);

    // unsigned int subdt = pic.getNumberOfSubTimeSteps(dt);
    // for (size_t i = 0; i < subdt; i++) {
    //   boundaryHandler->markCells(domain.get(), particles.get());
    //   boundaryHandler->markFaces(domain.get());
    //   pic.runTransferFromParticlesToGrid();
    //   pic.runComputeExternalForces(dt);
    //   pic.runApplyBoundaryCondition();
    //   boundaryHandler->markCells(domain.get(), particles.get());
    //   boundaryHandler->markFaces(domain.get());
    //   pic.runExtrapolateToDomain();
    //   pic.runComputePressure(dt);
    //   domainModel.updatePressure();
    //   pic.runCorrectVelocities(dt);
    //   domainModel.updateDivergence();
    //   pic.runExtrapolateToDomain();
    //   pic.runApplyBoundaryCondition();
    //   pic.runTransferFromGridToParticles();
    //   boundaryHandler->markCells(domain.get(), particles.get());
    //   boundaryHandler->markFaces(domain.get());
    //   pic.runComputeAdvection(dt);
    // }

    char filename[100];
    sprintf(filename, "picqt/pg%06d.png", curFrame);
    saveFramebuffer(filename, WIDTH, HEIGHT);
    curFrame++;
  };
  app.run();
  return 0;
}
