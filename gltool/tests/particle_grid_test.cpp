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
double maxError = 0., minError = INFINITY;
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
    if (node.level() > 4)
      return false;
    return true;
    //for (auto &p : samples)
    //  if (node.region().inside(p))
    //    return true;
    //return false;
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

void updateParticleVelocity(PICSolverDebuggerClass &pic) {
  for (size_t i = 0; i < numericalParticles->size(); i++) {
    auto p = (*numericalParticles)[i];
    numericalParticles->setScalarProperty(0, i, velocity(p).x());
    numericalParticles->setScalarProperty(1, i, velocity(p).y());
  }
  bh->markCells(domain.get(), numericalParticles);
  bh->markFaces(domain.get());
  pic.runTransferFromParticlesToGrid();
  domain->iterateDxScalar(0, [&](double &v, Point2d p) {
    double error = std::abs(velocity(p).x() - v);
    minError = std::min(minError, error);
    maxError = std::max(maxError, error);
  });
  domain->iterateDyScalar(0, [&](double &v, Point2d p) {
    double error = std::abs(velocity(p).y() - v);
    minError = std::min(minError, error);
    maxError = std::max(maxError, error);
  });
}

void step(PICSolverDebuggerClass &pic, double dt, uint curFrame) {
  // mark grid with numerical particles
  bh->markCells(domain.get(), numericalParticles);
  bh->markFaces(domain.get());
  // update analytical particles velocities
  for (size_t i = 0; i < numericalParticles->size(); i++) {
    auto p = (*numericalParticles)[i];
    numericalParticles->setScalarProperty(0, i, velocity(p).x() * ((curFrame >= 400) ? -1. : 1.));
    numericalParticles->setScalarProperty(1, i, velocity(p).y() * ((curFrame >= 400) ? -1. : 1.));
  }
  // numerical particles -> grid
  pic.runTransferFromParticlesToGrid();
  // advect numerical particles
  pic.runComputeAdvection(dt);
  // advect analytical particles
  analyticalParticles->iterateParticles([&](int i, Point2d p) {
    double tdt = dt;
    if (curFrame >= 400)
      tdt = -dt;
    auto np = advectExact(p, tdt, velocity);
    analyticalParticles->setPosition(i, np);
  });
}

int main() {
  aergia::SceneApp<> app(WIDTH, HEIGHT, "PIC 2D", false);
  std::cerr << "init aergia...";
  app.init();
  std::cerr << "finished.\n";
  // SIMULATION
  auto h = 1. / 128.;
  auto region = BBox2D::make_unit_bbox();
  std::cerr << "create particle system and injector.\n";
  numericalParticles = new ParticleSystem2(new PointZGrid2(32u, 32u, region),
                                           new RBFInterpolant<Point2d>(new QuinticKernel<Point2d>()));
  std::cerr << "inject.\n";
  Injector injector(numericalParticles, h);
  injector.setupCircleShape(Point2d(0.5, 0.75), 0.15);
  analyticalParticles = new ParticleSystem2(new PointZGrid2(32u, 32u, region),
                                            new RBFInterpolant<Point2d>(new QuinticKernel<Point2d>()));
  std::cerr << "add property\n";
  analyticalParticles->addScalarProperty();
  analyticalParticles->addScalarProperty();
  numericalParticles->iterateParticles([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(i);
    analyticalParticles->addParticle(p);
  });
  std::cerr << "build quadtree\n";
  setupStructures();
  std::cerr << "create pic\n";
  PICSolverDebuggerClass pic(bh.get(), domain.get(), numericalParticles);
  double dt = 1. / 200.;
  // GRAPHIC MODELS
  ParticleSystem2Model psm(numericalParticles, h / 3.);
  ParticleSystem2Model epsm(analyticalParticles, h / 5.);
  epsm.particleColor = aergia::COLOR_RED;
  epsm.particleColor.a = 0.3;
  SimulationQuadTreeDomainModel dm(domain.get(), bh.get());
  SearchRenderer sr(numericalParticles, domain.get());
  // INIT
  app.addViewport2D(0, 0, WIDTH, HEIGHT);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(region), 0.9f);
  std::cerr << "draw!\n";
  app.keyCallback = [&](int k, int a) {
    if (a != GLFW_PRESS)
      return;
    switch (k) {
    case GLFW_KEY_Q:app.exit();
      break;
    default:step(pic, dt, 0);
      updateParticleVelocity(pic);
      ponos::BBox2D cameraRegion((*numericalParticles)[0].x(), (*numericalParticles)[1].y());
      numericalParticles->iterateParticles([&](uint i, Point2d p) {
        UNUSED_VARIABLE(i);
        ponos::make_union(cameraRegion, ponos::Point2(p.x(), p.y()));
      });
      app.getCamera<aergia::Camera2D>(0)->fit(cameraRegion, 1.1f);
      break;
    }
  };
  // updateParticleVelocity(pic);
  int curFrame = 0;
 // ponos::BBox2D cameraRegion((*numericalParticles)[0].x(), (*numericalParticles)[1].y());
  //numericalParticles->iterateParticles([&](uint i, Point2d p) {
  //  UNUSED_VARIABLE(i);
  //  ponos::make_union(cameraRegion, ponos::Point2(p.x(), p.y()));
  //});
  //app.getCamera<aergia::Camera2D>(0)->fit(cameraRegion, 1.1f);

  app.renderCallback = [&]() {
    /*domain->iterateDxScalar(0, [&](size_t id, double &v) {
      if (bh->faceType(id) != DomainBoundaryHandler2::MaterialType::FLUID)
        return;
      Point2d p = domain->faceCenterPosition(id);
      double error = std::abs(velocity(p).x() - v);
      double d = h / 3;
      ponos::BBox2D bb(ponos::Point2(p.x() - d, p.y() - d), ponos::Point2(p.x() + d, p.y() + d));
      aergia::draw_bbox(bb, aergia::COLOR_BLACK, aergia::Color(0, 0, 0, (error - minError) / (maxError - minError)));
    });
    domain->iterateDyScalar(0, [&](size_t id, double &v) {
      if (bh->faceType(id) != DomainBoundaryHandler2::MaterialType::FLUID)
        return;
      Point2d p = domain->faceCenterPosition(id);
      double error = std::abs(velocity(p).y() - v);
      double d = h / 3;
      ponos::BBox2D bb(ponos::Point2(p.x() - d, p.y() - d), ponos::Point2(p.x() + d, p.y() + d));
      aergia::draw_bbox(bb, aergia::COLOR_BLACK, aergia::Color(0, 0, 0, (error - minError) / (maxError - minError)));
    });*/
    if (curFrame > 799)
      return;
    step(pic, dt, curFrame);
    char filename[100];
    sprintf(filename, "simulation/pg%06d.png", curFrame);
    saveFramebuffer(filename, WIDTH, HEIGHT);
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