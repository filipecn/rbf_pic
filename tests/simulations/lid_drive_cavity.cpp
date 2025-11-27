#include <furoo.h>

using namespace furoo;

int main(int argc, char const *argv[]) {

  // Simulation structures
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, 4 /*minlevel*/, 6 /*maxlevel*/);
  auto particles = simulation.particles();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());

  // Initializing setup
  Injector2 injector(particles, 0.00390625 /*particle spacing*/);
  injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.999, 0.99)));
  simulation.buildTree();
  simulation.updateTree();

  // Simulation steps

  return 0;
}
