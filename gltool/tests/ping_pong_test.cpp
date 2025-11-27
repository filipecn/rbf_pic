#include "particle_system_model.h"
#include "cell_centered_graph2_model.h"
#include "test_utils.h"

using namespace furoo;

int main() {
  ParticleSystem2d exactParticleSystem(new PointZGrid2d(16));
  {
    exactParticleSystem.addScalarProperty(0.);
    exactParticleSystem.addScalarProperty(0.);
    Injector2 injector(&exactParticleSystem, 0.007);
    injector.setupBoxShape(BBox2d(Point2d(0.2, 0.7), Point2d(0.6, 0.9)));
  }
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 3, 5);
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
  app.getCamera<aergia::Camera2D>(0)->fit(Ponos::bbox2D(domainRegion));
  ParticleSystemModel2 particleSystemModel(particles, 0.003);
  ParticleSystemModel2 exactParticleSystemModel(&exactParticleSystem, 0.003);
  exactParticleSystemModel.particleColor = aergia::COLOR_GREEN;
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  auto &faceVelocityField = *structure->field<double>(10);
  CellCenteredGraph2Model
      cellGraph2Model(dynamic_cast<CellGraph2 *>(structure));
  app.scene.add(&cellGraph2Model);
  app.scene.add(&exactParticleSystemModel);
  app.scene.add(&particleSystemModel);
  bool paused = false;
  app.keyCallback = [&](int k, int a) {
    if (a == GLFW_RELEASE && k == GLFW_KEY_SPACE)
      paused = !paused;
  };
  app.renderCallback = [&]() {
    if (paused)
      return;
    simulation.updateTree();
    exactParticleSystem.iterateParticles([&](size_t id, Point2d p) {
      exactParticleSystem.setPosition(id,
                                      advectExact(p,
                                                  0.001,
                                                  Functions::enright));
    });
    // set particles velocities from analytical function
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setScalarProperty(0, id, Functions::enright(p).x());
      particles->setScalarProperty(1, id, Functions::enright(p).y());
    });
//    simulation.transferFromParticlesToGrid();
    structure->iterateFaces([&](size_t i) {
      auto v = Functions::enright(structure->faceCenterPosition(i));
      if (structure->faceOrientation(i) == Definitions::Orientation::HORIZONTAL)
        faceVelocityField[i] = v.y();
      else
        faceVelocityField[i] = v.x();
    });
    simulation.transferFromGridToParticles();
    particles->iterateParticles([&](size_t id, Point2d p) {
      particles->setPosition(id, advectExact(p,
                                             0.001,
                                             [&](Point2d point) -> Vector2d {
                                               UNUSED_VARIABLE(point);
                                               return Vector2d(particles->getScalarProperty(
                                                   0, id),
                                                               particles->getScalarProperty(
                                                                   1, id));
                                             }));
    });
    // simulation loop
    //simulation.markFaces();
    //simulation.computeExternalForces(0.001);
    //simulation.applyBoundaryCondition();
  };

  app.run();
}