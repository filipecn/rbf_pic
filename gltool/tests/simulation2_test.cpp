#include "particle_system_model.h"
#include "cell_centered_graph2_model.h"
#include "test_utils.h"

using namespace furoo;

int main() {
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 5, 6);
  simulation.setRegion(domainRegion);
  auto particles = simulation.particles();
  {
    Injector2 injector(particles, 0.007);
    injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  }
  simulation.buildTree();
  aergia::SceneApp<> app(800, 800, "", false);
  app.init();
  app.addViewport2D(0, 0, 800, 800);
  ParticleSystemModel2 particleSystemModel(particles, 0.003);
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &velocityXField = *structure->field<double>(8);
  auto &velocityYField = *structure->field<double>(9);
  CellCenteredGraph2Model
      cellGraph2Model(dynamic_cast<CellGraph2 *>(structure));
  app.scene.add(&cellGraph2Model);
  app.scene.add(&particleSystemModel);
  app.renderCallback = [&]() {
    simulation.updateTree();
    structure->iterateCells([&](size_t i) {
      auto v = furoo::Functions::enright(structure->cellCenterPosition(i));
      velocityXField[i] = v.x();
      velocityYField[i] = v.y();
    });
    simulation.transferFromGridToParticles();
    simulation.advectParticles(0.001);
  };

  app.run();
  return 0;
}