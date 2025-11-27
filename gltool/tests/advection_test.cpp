#include <furoo.h>
#include "boundary_handler2_model.h"
#include "particle_system2_model.h"
#include "ponos_wrapper.h"
#include "simulation_quad_tree_domain_model.h"
#include "SearchRenderer.h"
#include <aergia/aergia.h>
#include "test_utils.h"

using namespace furoo;

#define WIDTH 1000
#define HEIGHT 1000

ParticleSystem2 *numericalParticles;
std::shared_ptr<DomainBoundaryHandler2> bh;
std::shared_ptr<SimQuadTree> domain;
// exact
ParticleSystem2 *analyticalParticles;
// field
Vector2d velocity(Point2d p) {
  return Vector2d(
      2. * SQR(sin(PI * p.x())) * sin(PI * p.y()) * cos(PI * p.y()),
      -2. * SQR(sin(PI * p.y())) * sin(PI * p.x()) * cos(PI * p.x()));
  // return Vector2d(-(PI / 2.) * (0.5 - p.y()), -(PI / 2.) * (p.x() - 0.5));
}

void setupStructures() {

  HaltonSequence rngX(3), rngY(5);
  std::vector<Point2d> samples;
  for (int i = 0; i < 47; i++)
    samples.emplace_back(rngX.randomDouble(), rngY.randomDouble());

  domain.reset(new SimQuadTree(BBox2D::make_unit_bbox(), [&](QuadTree::Node &node) -> bool {
    if (node.level() > 6)
      return false;
    return true;

    for (auto &p : samples)
      if (node.region().inside(p))
        return true;
    return false;
    //if (node.level() < 2)
    //  return true;
    //return static_cast<size_t>(ponos::lerp(node.region().center().x(), 3.f, 6.f)) >= node.level();
    int count = 0;
    numericalParticles->iterateParticles(node.region(), [&](unsigned int i, Point2d p) {
      UNUSED_VARIABLE(p);
      UNUSED_VARIABLE(i);
      count++;
    });
    return count > 0 && node.level() < 7;
  }));
  domain->setPointSetStructures(new PointZGrid2(16), new PointZGrid2(16),
                                new PointZGrid2(16));
  bh.reset(new DomainBoundaryHandler2(DomainBoundaryHandler2::BoundaryType::NEUMANN));
  // fill velocity field
  std::cout << domain->cellCount() << std::endl;
}

int main() {
  aergia::SceneApp<> app(WIDTH, HEIGHT, "PIC 2D", false);
  app.init();
  // SIMULATION
  auto h = 1. / 128.;
  auto region = BBox2D::make_unit_bbox();
  numericalParticles = new ParticleSystem2(new PointZGrid2(16u, 16u, region),
                           new RBFInterpolant<Point2d>(new GaussianKernel<Point2d>(2.)));
  Injector injector(numericalParticles, h);
  injector.setupCircleShape(Point2d(0.5, 0.75), 0.15);
  analyticalParticles = new ParticleSystem2(new PointZGrid2(16u, 16u, region),
                            new RBFInterpolant<Point2d>(new GaussianKernel<Point2d>(2.)));
  setupStructures();
  PICSolverDebuggerClass pic(bh.get(), domain.get(), numericalParticles);
  //injector.loadFromFile("simulation/particles_000360");
  numericalParticles->iterateParticles([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(i);
    analyticalParticles->addParticle(p);
  });
  double dt = 1. / 200.;
  // GRAPHIC MODELS
  ParticleSystem2Model psm(numericalParticles, h / 4.);
  ParticleSystem2Model epsm(analyticalParticles, h / 5.);
  epsm.particleColor = aergia::COLOR_RED;
  epsm.particleColor.a = 0.4;
  SimulationQuadTreeDomainModel dm(domain.get(), bh.get());
  SearchRenderer sr(numericalParticles, domain.get());
  domain->iterateDxScalar(0, [](double &v, Point2d p) {
    v = velocity(p).x();
  });
  domain->iterateDyScalar(0, [](double &v, Point2d p) {
    v = velocity(p).y();
  });
  // INIT
  app.addViewport2D(0, 0, WIDTH, HEIGHT);
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
    case GLFW_KEY_4:bh->markCells(domain.get(), numericalParticles);
      break;
    case GLFW_KEY_5:bh->markFaces(domain.get());
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
      break;
    default:
      break;
    }
  };
  int curFrame = 0;
  app.renderCallback = [&]() {
    if (curFrame == 400) {
      domain->iterateDxScalar(0, [](double &v, Point2d p) {
        v = -velocity(p).x();
      });
      domain->iterateDyScalar(0, [](double &v, Point2d p) {
        v = -velocity(p).y();
      });
    }
    if (curFrame > 799)
      return;
    // initiate structures
    /*setupStructures();
    std::cout << domain->region();
    pic.setSimulationDomain(domain.get());
    pic.setBoundaryHandler(bh.get());
    domain->iterateDxScalar(0, [](double &v, Point2d p) {
      v = -(PI / 2.) * (0.5 - p.y());
    });
    domain->iterateDyScalar(0, [](double &v, Point2d p) {
      v = -(PI / 2.) * (p.x() - 0.5);
    });
    // update graphic models
    dm.set(domain.get(), bh.get());
    sr.set(ps, domain.get());
     */
    bh->markCells(domain.get(), numericalParticles);
    bh->markFaces(domain.get());
    analyticalParticles->iterateParticles([&](int i, Point2d p) {
      Point2d np = advectExact(p, dt, velocity);
      if(curFrame >= 400)
        np = advectExact(p, -dt, velocity);
      analyticalParticles->setPosition(i, np);
    });
    pic.runComputeAdvection(dt);
    char filename[100];
    sprintf(filename, "simulation/%06d.png", curFrame);
    saveFramebuffer(filename, WIDTH, HEIGHT);
    sprintf(filename, "simulation/particles_%06d", curFrame);
    IO::saveToPV(filename, numericalParticles, true);
    curFrame++;
  };
  app.scene.add(&psm);
  app.scene.add(&epsm);
  app.scene.add(&dm);
  app.scene.add(&sr);
  app.run();
  delete numericalParticles;
  return 0;
}