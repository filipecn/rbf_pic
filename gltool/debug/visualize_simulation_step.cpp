#include "cell_centered_graph2_model.h"
#include "particle_system_model.h"
#include "test_utils.h"

using namespace furoo;

int main(int argc, char **argv) {
  UNUSED_VARIABLE(argc);
  furoo::ArgsParser parser;
  parser.addStringArgument("particlesFile", "");
  parser.addStringArgument("cellGraphFile", "");
  parser.parse(argv[1]);

  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 4, 6);
  auto particles = simulation.particles();
  IO::loadFromPV(parser.getString("particlesFile").c_str(), particles);
  simulation.loadFromFile(parser.getString("cellGraphFile").c_str());
  aergia::SceneApp<> app(1000, 1000, "", false);

  app.init();
  app.addViewport2D(0, 0, 1000, 1000);
  app.getCamera<aergia::UserCamera2D>(0)->fit(Ponos::bbox2D(domainRegion));
  ParticleSystemModel2 particleSystemModel(particles, 0.001);
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  CellCenteredGraph2Model cellGraph2Model(
      dynamic_cast<CellGraph2 *>(structure));
  app.scene.add(&particleSystemModel);
  app.scene.add(&cellGraph2Model);
  //  static int step = 9;

  app.renderCallback = [&]() {
    glPointSize(5.f);
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    glVertex2f(0.779779, 0.00117042);
    glEnd();
    glPointSize(1.f);
    return;
  };
  app.run();
  return 0;
}
