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
  // DROP
  // injector.setupCircleShape(Point2d(0.5, 0.8), 0.08);
  // injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.999, 0.2));
  // LEFT DAM
  injector.setupBoxShape(Point2d(0.001, 0.001), Point2d(.2, 0.999));
  // RIGHT DAM
  // injector.setupBoxShape(Point2d(0.9, 0.001), Point2d(.999, 0.999));
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
  app.keyCallback = [&](int k, int a) {
    if (a != GLFW_RELEASE)
      return;
    if (k == GLFW_KEY_Q)
      app.exit();
    if (k == GLFW_KEY_P) {
      domainModel.drawPressure();
    } else if (k == GLFW_KEY_D) {
      domainModel.updateDivergence();
      domainModel.drawDivergence();
    } else if (k == GLFW_KEY_M) {
      domainModel.drawMaterial();
    } else if (k == GLFW_KEY_1) {
      std::cerr << "mark cells and faces\n";
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
    } else if (k == GLFW_KEY_2) {
      std::cerr << "compute advection\n";
      pic.runComputeAdvection(dt);
    } else if (k == GLFW_KEY_3) {
      std::cerr << "mark cells and faces\n";
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
    } else if (k == GLFW_KEY_4) {
      std::cerr << "particles -> grid\n";
      pic.runTransferFromParticlesToGrid();
    } else if (k == GLFW_KEY_5) {
      std::cerr << "external forces\n";
      pic.runComputeExternalForces(dt);
    } else if (k == GLFW_KEY_6) {
      std::cerr << "boundary condition\n";
      pic.runApplyBoundaryCondition();
    } else if (k == GLFW_KEY_7) {
      std::cerr << "mark cells and faces\n";
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
    } else if (k == GLFW_KEY_8) {
      std::cerr << "extrapolate\n";
      pic.runExtrapolateToDomain();
    } else if (k == GLFW_KEY_9) {
      std::cerr << "computePressure .. ";
      pic.runComputePressure(dt);
      domainModel.updatePressure();
      std::cerr << " finished.\n";
    } else if (k == GLFW_KEY_0) {
      std::cerr << " correct velocities.\n";
      pic.runCorrectVelocities(dt);
      domainModel.updateDivergence();
    } else if (k == GLFW_KEY_Z) {
      std::cerr << "extrapolate\n";
      pic.runExtrapolateToDomain();
    } else if (k == GLFW_KEY_X) {
      std::cerr << "boundary conditions\n";
      pic.runApplyBoundaryCondition();
    } else if (k == GLFW_KEY_C) {
      std::cerr << "grid -> particles\n";
      pic.runTransferFromGridToParticles();
    } else if (k == GLFW_KEY_SPACE) {
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runTransferFromParticlesToGrid();
      pic.runComputeExternalForces(dt);
      pic.runApplyBoundaryCondition();
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runExtrapolateToDomain();
      pic.runComputePressure(dt);
      domainModel.updatePressure();
      pic.runCorrectVelocities(dt);
      domainModel.updateDivergence();
      pic.runExtrapolateToDomain();
      pic.runApplyBoundaryCondition();
      pic.runTransferFromGridToParticles();
      boundaryHandler->markCells(domain.get(), particles.get());
      boundaryHandler->markFaces(domain.get());
      pic.runComputeAdvection(dt);
      static int curFrame = 0;
      char filename[100];
      sprintf(filename, "picqt/pg%06d.png", curFrame);
      saveFramebuffer(filename, WIDTH, HEIGHT);
      curFrame++;
    }
  };
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
        boundaryHandler->markCells(domain.get(), particles.get());
        boundaryHandler->markFaces(domain.get());
        pic.runTransferFromParticlesToGrid();
        pic.runComputeExternalForces(DT);
        pic.runApplyBoundaryCondition();
        boundaryHandler->markCells(domain.get(), particles.get());
        boundaryHandler->markFaces(domain.get());
        pic.runExtrapolateToDomain();
        pic.runComputePressure(DT);
        domainModel.updatePressure();
        pic.runCorrectVelocities(DT);
        domainModel.updateDivergence();
        pic.runExtrapolateToDomain();
        pic.runApplyBoundaryCondition();
        pic.runTransferFromGridToParticles();
        particlesModel.update();
        boundaryHandler->markCells(domain.get(), particles.get());
        boundaryHandler->markFaces(domain.get());
        pic.runComputeAdvection(DT);
      } catch (const char *e) {
        std::cerr << "Thrown from " << e << '\n';
        exit(-1);
      }
    }

    char filename[100];
    sprintf(filename, "picqt/pg%06d.png", curFrame);
    saveFramebuffer(filename, WIDTH, HEIGHT);
    curFrame++;
  };
  app.run();
  return 0;
}
