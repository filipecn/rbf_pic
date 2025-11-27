#include "particle_system_model.h"
#include "cell_graph2_model.h"
#include "test_utils.h"

using namespace furoo;

int main() {
  auto domainRegion = BBox2d::squareBox();
  AdaptivePic2 simulation(domainRegion, 6, 7);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  {
    Injector2 injector(particles, 0.007);
//    injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.2, 0.5)));
    injector.setupBoxShape(BBox2d(Point2d(0.45, 0.6), Point2d(0.55, 0.7)));
    injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.999, 0.2)));
//    injector.setupBoxShape(BBox2d(Point2d(0.4, 0.001), Point2d(0.6, 0.2)));
  }
  simulation.buildTree();
  aergia::SceneApp<> app(800, 800, "", false);
  app.init();
  app.addViewport2D(0, 0, 800, 800);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(domainRegion));
  ParticleSystemModel2 particleSystemModel(particles, 0.003);
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  CellGraph2Model cellGraph2Model(dynamic_cast<CellGraph2 *>(structure));
//  cellGraph2Model.debugMode = true;
  app.scene.add(&cellGraph2Model);
  app.scene.add(&particleSystemModel);
  app.keyCallback = [&](int k, int a) {
    UNUSED_VARIABLE(k);
    if (a == GLFW_RELEASE) {
      static int step = 0;
      TimeLogger timer;
      double dt = 0.001;
      switch (step) {
        case 0: simulation.updateTree();
          timer.log("update");
          break;
        case 1: simulation.markCells();
          timer.log("markCells");
          break;
        case 2: simulation.markFaces();
          timer.log("markFaces");
          break;
        case 3: simulation.transferFromParticlesToGrid();
          timer.log("part -> grid");
          break;
        case 4: simulation.computeExternalForces(dt);
          timer.log("gravity");
          break;
        case 5: simulation.applyBoundaryCondition();
          timer.log("boundary");
          break;
        case 6: simulation.solvePressure(dt);
          timer.log("pressure");
          cellGraph2Model.updatePressure();
          break;
        case 7: simulation.transferFromGridToParticles();
          timer.log("grid -> part");
          break;
        case 8: simulation.advectParticles(dt);
          timer.log("advect");
          break;
        default: step = 0;
      }
      timer.report();
      step++;
    }
  };
  app.renderCallback = [&]() {
//    static int first = 0;
//    if (first++)
//      return;
    double dt = 1 / 60.;
    // simulation loop
    size_t cflSteps = std::ceil(std::max(dt / simulation.cfl(), 1.));
    double cflDT = dt / cflSteps;
    std::cerr << cflSteps << " substeps, " << cflDT << std::endl;
    TimeLogger timer;
    for (size_t i = 0; i < cflSteps; i++) {
      simulation.updateTree();
      timer.log("update");
      simulation.markCells();
      timer.log("markCells");
      simulation.markFaces();
      timer.log("markFaces");
      simulation.transferFromParticlesToGrid();
      timer.log("part -> grid");
      simulation.computeExternalForces(cflDT);
      timer.log("gravity");
      simulation.applyBoundaryCondition();
      timer.log("boundary");
      simulation.solvePressure(cflDT);
      timer.log("pressure");
      cellGraph2Model.updatePressure();
      simulation.transferFromGridToParticles();
      timer.log("grid -> part");
      simulation.advectParticles(cflDT);
      timer.log("advect");
    }
    timer.report();
  };
  app.run();
  return 0;
}