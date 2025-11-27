#include "cell_centered_graph2_model.h"
#include "particle_system_model.h"
#include "test_utils.h"

using namespace furoo;

int main(int argc, char **argv) {
  UNUSED_VARIABLE(argc);
  UNUSED_VARIABLE(argv);
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 4, 6);
  auto particles = simulation.particles();
  Injector2 injector2(particles, 0.00390625);
//  injector2.setupCircleShape(Point2d(0.5), 0.3);
  //injector2.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.2, 0.3)));
  injector2.setupBoxShape(BBox2d(Point2d(0.8, 0.001), Point2d(0.999, 0.5)));
  simulation.buildTree();
  simulation.updateTree();
  aergia::SceneApp<> app(1000, 1000, "", false);
  app.init();
  app.addViewport2D(0, 0, 1000, 1000);
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(domainRegion));
  ParticleSystemModel2 particleSystemModel(particles, 0.001);
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  CellCenteredGraph2Model cellGraph2Model(
      dynamic_cast<CellGraph2 *>(structure));
  app.scene.add(&particleSystemModel);
  app.scene.add(&cellGraph2Model);
  app.renderCallback = [&]() {
  };
  app.run();
  return 0;
}
