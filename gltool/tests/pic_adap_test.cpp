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

void updateQuadTree(PICSolver2 *pic, BBox2D region, size_t quadtreeLevel) {
  pic->setSimulationDomain(new SimQuadTree(region, [&](QuadTree::Node &node) -> bool {
    int count = 0;
    pic->particles()->iterateParticles(node.region(), [&count](unsigned int id, Point2d p) {
      UNUSED_VARIABLE(id);
      UNUSED_VARIABLE(p);
      count++;
    });
    return count > 0 && node.level() < quadtreeLevel;
  }));
  dynamic_cast<SimQuadTree *>(pic->domain())
      ->setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                              new PointZGrid2(16));
  pic->setBoundaryHandler(new DomainBoundaryHandler2(
      DomainBoundaryHandler2::BoundaryType::NEUMANN));
  pic->domain()->setInterpolationMethod(
      SimulationDomain2::InterpolationMethod::RBF);
}

int main(int argc, char **argv) {
  // Configuration
  if (argc < 5) {
    std::cout
        << "Usage: <quadtree level> <particle spacing> <# of frames> <dt>\n";
    return 0;
  }
  size_t quadtreeLevel = 0;
  double particleSpacing = 0.;
  size_t nframes = 0;
  double dt = 0.;
  sscanf(argv[1], "%zu", &quadtreeLevel);
  sscanf(argv[2], "%lf", &particleSpacing);
  sscanf(argv[3], "%zu", &nframes);
  sscanf(argv[4], "%lf", &dt);
  // Initialization
  BBox2D region = BBox2D::make_unit_bbox();
  std::shared_ptr<ParticleSystem2> particles(new ParticleSystem2(
      new PointZGrid2(32u, 32u, region), new ShepardInterpolant2()));
  // Inject Particles
  Injector injector(particles.get(), particleSpacing);
  // DROP
  injector.setupCircleShape(Point2d(0.5, 0.8), 0.08);
  injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.999, 0.2));
  // LEFT DAM
  // injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.2, 0.999));
  // RIGHT DAM
  // injector.setupBoxShape(Point2d(0.9, 0.001), Point2d(.999, 0.999));
  // PIC
  PICSolverDebuggerClass pic;
  pic.setParticleSystem(particles.get());
  updateQuadTree(&pic, region, quadtreeLevel);
  // VISUALIZATION
  aergia::SceneApp<> app(WIDTH, HEIGHT, "PIC ADAP 2D", false);
  app.init();
  app.addViewport2D(0, 0, WIDTH, HEIGHT);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(region), 1.1f);
  // HELPERS
  ParticleSystem2Model particlesModel(particles.get(), 0.003);
  SimulationQuadTreeDomainModel domainModel(dynamic_cast<SimQuadTree *>(pic.domain()),
                                            pic.boundaryHandler());
  app.scene.add(&particlesModel);
  app.scene.add(&domainModel);
  // for (size_t i = 0; i < 5; i++) {
  //   pic.runAdvanceTimeStep(dt);
  // }
  particlesModel.update();
  app.renderCallback = [&]() {
    static int curFrame = 0;
    if (curFrame > 799)
      return;
    // return;
    // pic.runAdvanceTimeStep(dt);
    unsigned int subdt = pic.getNumberOfSubTimeSteps(dt);
    double DT = dt / subdt;
    for (size_t i = 0; i < subdt; i++) {
      try {
        updateQuadTree(&pic, region, quadtreeLevel);
        pic.boundaryHandler()->markCells(pic.domain(), particles.get());
        pic.boundaryHandler()->markFaces(pic.domain());
        pic.runTransferFromParticlesToGrid();
        pic.runComputeExternalForces(DT);
        pic.runApplyBoundaryCondition();
        pic.boundaryHandler()->markCells(pic.domain(), particles.get());
        pic.boundaryHandler()->markFaces(pic.domain());
        pic.runExtrapolateToDomain();
        pic.runComputePressure(DT);
        //domainModel.updatePressure();
        pic.runCorrectVelocities(DT);
        //domainModel.updateDivergence();
        pic.runExtrapolateToDomain();
        pic.runApplyBoundaryCondition();
        pic.runTransferFromGridToParticles();
        pic.boundaryHandler()->markCells(pic.domain(), particles.get());
        pic.boundaryHandler()->markFaces(pic.domain());
        pic.runComputeAdvection(DT);
      } catch (const char *e) {
        std::cerr << "Thrown from " << e << '\n';
        exit(-1);
      }
    }

    domainModel.set(dynamic_cast<SimQuadTree *>(pic.domain()), pic.boundaryHandler());
    particlesModel.update();

    char filename[100];
    sprintf(filename, "picqt/pg%06d.png", curFrame);
    saveFramebuffer(filename, WIDTH, HEIGHT);
    curFrame++;
  };
  app.run();
  return 0;
}
