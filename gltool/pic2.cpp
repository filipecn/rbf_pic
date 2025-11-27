#include <furoo.h>
#include "boundary_handler2_model.h"
#include "particle_system2_model.h"
#include "ponos_wrapper.h"
#include "simulation_regular_domain_model.h"
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

int main(int argc, char **argv) {
  bool loadingFromFile = false;
  aergia::SceneApp<> app(800, 800, "PIC 2D", false);
  app.init();
  char path[20] = "simulation2/";
  if (argc > 1) {
    std::strcpy(path, argv[1]);
    loadingFromFile = true;
  }
  // SIMULATION
  int frameNumber = 0;
  auto h = 1. / 100.;
  auto dt = 1. / 60.;
  auto region = BBox2D::make_unit_bbox();
  auto dm = new SimRegularDomain2(100);
  dm->setInterpolationMethod(SimulationDomain2::InterpolationMethod::BILINEAR);
  auto ps = new ParticleSystem2(new PointZGrid2(16u, 16u, region));
  if (!loadingFromFile) {
    for (size_t i = 0; i < 3; i++)
      for (size_t j = 0; j < 3; j++) {
        ps->addParticle(
            Point2d(i * h + h / 4. + h / 2., j * h + h / 4. + h / 2.));
        ps->addParticle(
            Point2d(i * h - h / 4. + h / 2., j * h + h / 4. + h / 2.));
        ps->addParticle(
            Point2d(i * h + h / 4. + h / 2., j * h - h / 4. + h / 2.));
        ps->addParticle(
            Point2d(i * h - h / 4. + h / 2., j * h - h / 4. + h / 2.));
      }
  }
  auto bh =
      new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN);
  PICSolverDebuggerClass pic(bh, dm, ps);
  bh->markFaces(dm);
  // GRAPHIC MODELS
  SimulationRegularDomain2Model sdm(dm);
  DomainBoundaryHandler2Model bhm(dm, bh);
  ParticleSystem2Model psm(ps, h / 8.);
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
    case GLFW_KEY_Q:app.exit();
      break;
    case GLFW_KEY_1:pic.runBeginAdvanceTimeStep(dt);
      break;
    case GLFW_KEY_2:pic.runComputeExternalForces(dt);
      break;
    case GLFW_KEY_3:pic.runApplyBoundaryCondition();
      break;
    case GLFW_KEY_4:bh->markCells(dm, ps);
      break;
    case GLFW_KEY_5:bh->markFaces(dm);
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
    case GLFW_KEY_O:pic.runEndAdvanceTimeStep(dt);
      break;
    case GLFW_KEY_SPACE:pic.runBeginAdvanceTimeStep(dt);

      pic.runComputeExternalForces(dt);

      pic.runApplyBoundaryCondition();

      bh->markCells(dm, ps);

      bh->markFaces(dm);

      pic.runExtrapolateToDomain();

      pic.runComputePressure(dt);

      pic.runCorrectVelocities(dt);

      pic.runExtrapolateToDomain();

      pic.runApplyBoundaryCondition();

      pic.runTransferFromGridToParticles();

      pic.runComputeAdvection(dt);

      pic.runEndAdvanceTimeStep(dt);
      break;

    case GLFW_KEY_M: {
      frameNumber++;
      char filename[100];
      sprintf(filename, "%s%06d", path, frameNumber);
      std::cout << "loading " << filename << std::endl;
      IO::loadFromPV(filename, ps);
      bh->markCells(dm, ps);
      bh->markFaces(dm);
      IO::loadFromCF(filename, dm, 0);
      sdm.updatePressure();
    }
      break;
    case GLFW_KEY_N: {
      char filename[100];
      frameNumber = std::max(0, frameNumber - 1);
      sprintf(filename, "%s%06d", path, frameNumber);
      IO::loadFromPV(filename, ps);
      bh->markCells(dm, ps);
      bh->markFaces(dm);
      IO::loadFromCF(filename, dm, 0);
      sdm.updatePressure();
    }
      break;
    case GLFW_KEY_P:sdm.drawPressure = !sdm.drawPressure;
      break;
    case GLFW_KEY_D:sdm.drawDivergence = !sdm.drawDivergence;
      if (sdm.drawDivergence)
        sdm.updateDivergence();
      break;
    case GLFW_KEY_C:bhm.drawCells = !bhm.drawCells;
      break;
    case GLFW_KEY_V:sdm.drawVelocities = !sdm.drawVelocities;
      break;
    case GLFW_KEY_G:pic.runTransferFromParticlesToGrid();
      break;
    case GLFW_KEY_E:pic.runExtrapolateToDomain();
      break;
    default:break;
    }
  };
  app.scene.add(&sdm);
  app.scene.add(&bhm);
  app.scene.add(&psm);
  app.run();
  return 0;
}
