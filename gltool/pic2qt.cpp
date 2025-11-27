#include <furoo.h>
#include "boundary_handler2_model.h"
#include "particle_system2_model.h"
#include "ponos_wrapper.h"
#include "simulation_quad_tree_domain_model.h"
#include "SearchRenderer.h"
#include <aergia/aergia.h>

using namespace furoo;

class PICSolverDebuggerClass : public PICSolver2 {
public:
  PICSolverDebuggerClass() {}
  PICSolverDebuggerClass(DomainBoundaryHandler2 *b, SimulationDomain2 *d,
                         ParticleSystem2 *ps)
      : PICSolver2(b, d, ps) {}

  void runOnInitialize() { this->onInitialize(); }

  void runAdvanceFrame() { this->advanceFrame(); }

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

int main() {
  aergia::SceneApp<> app(800, 800, "PIC 2D", false);
  app.init();
  Timer timer;
  // SIMULATION
  auto h = 1. / 128.;
  auto region = BBox2D::make_unit_bbox();
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region),
                                new RBFInterpolant<Point2d>(new GaussianKernel<Point2d>(2.)));
  for (size_t i = 37; i < 69; i++)
    for (size_t j = 2; j < 15; j++) {
      ps->addParticle(
          Point2d(i * h + h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h - h / 4. + h / 2., j * h + h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h + h / 4. + h / 2., j * h - h / 4. + h / 2.));
      ps->addParticle(
          Point2d(i * h - h / 4. + h / 2., j * h - h / 4. + h / 2.));
    }
  std::cerr << timer.ellapsedTimeInSeconds() << "s --> particle system" << std::endl;
  timer.reset();
  SimQuadTree domain(BBox2D::make_unit_bbox(), [&](QuadTree::Node &node) -> bool {
    int count = 0;
    ps->iterateParticles(node.region(), [&](unsigned int i, Point2d p) {
      UNUSED_VARIABLE(p);
      UNUSED_VARIABLE(i);
      count++;
    });
    return count > 0 && node.level() < 7;
  });
  domain.setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                               new PointZGrid2(16));
  std::cerr << timer.ellapsedTimeInSeconds() << "s --> tree init" << std::endl;
  timer.reset();
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverDebuggerClass pic(bh, &domain, ps);
  double dt = 1. / 60.;
  // GRAPHIC MODELS
  ParticleSystem2Model psm(ps, h / 8.);
  SimulationQuadTreeDomainModel dm(&domain, bh);
  SearchRenderer sr(ps, &domain);
  std::cout << domain._vEdgeIndices.size() << std::endl;
  std::cout << domain._vFaceCenters->size() << std::endl;
  // INIT
  app.addViewport2D(0, 0, 800, 800);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(region), 1.2f);
  aergia::Text *text = new aergia::Text(FONT_PATH);
  app.viewports[0].renderCallback = [&]() {
    std::ostringstream stringStream;
    stringStream << "Q - exit | M - next | N - prev | P - pressure | V - vel | "
        "G - transf. 2 grid | E - extrapolate";
    text->render(stringStream.str(), ponos::Point3(-0.9f, 0.9f, 0), 0.4f,
                 aergia::Color(1.f, 0.f, 0.f, 0.5f));

  };
  app.keyCallback = [&](int k, int a) {
    if (a != GLFW_PRESS)
      return;
    switch (k) {
    case GLFW_KEY_P:
      for (size_t i = 0; i < ps->size(); i++) {
        ps->setScalarProperty(0, i, 0.1);
        ps->setScalarProperty(1, i, 0.2);
      }
      break;
    case GLFW_KEY_Q:app.exit();
      break;
    case GLFW_KEY_1:pic.runBeginAdvanceTimeStep(dt);
      break;
    case GLFW_KEY_2:pic.runComputeExternalForces(dt);
      break;
    case GLFW_KEY_3:pic.runApplyBoundaryCondition();
      break;
    case GLFW_KEY_4:bh->markCells(&domain, ps);
      break;
    case GLFW_KEY_5:bh->markFaces(&domain);
      break;
    case GLFW_KEY_6:pic.runExtrapolateToDomain();
      break;
    case GLFW_KEY_7:pic.runComputePressure(dt);
      break;
    case GLFW_KEY_8:pic.runCorrectVelocities(dt);
      break;
    case GLFW_KEY_9:pic.runExtrapolateToDomain();
      break;
    case GLFW_KEY_0:pic.runApplyBoundaryCondition();
      break;
    case GLFW_KEY_MINUS:pic.runTransferFromGridToParticles();
      break;
    case GLFW_KEY_EQUAL:pic.runComputeAdvection(dt);
      break;
    case GLFW_KEY_D:dm.updateDivergence();
      dm.drawDivergence = !dm.drawDivergence;
      break;
    case GLFW_KEY_O:pic.runEndAdvanceTimeStep(dt);
    default:break;
    }
  };
  app.scene.add(&psm);
  app.scene.add(&dm);
  app.scene.add(&sr);
  app.run();
  return 0;
}
