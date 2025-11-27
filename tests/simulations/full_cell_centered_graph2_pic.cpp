#include <furoo.h>
#include <iomanip>

using namespace furoo;

int main(int argc, char **argv) {
  UNUSED_VARIABLE(argc);
  furoo::ArgsParser parser;
  parser.addIntArgument("minLevel", 5);
  parser.addIntArgument("maxLevel", 7);
  parser.addStringArgument("config", std::string());
  parser.addStringArgument("particlesOutput", std::string("."));
  parser.addStringArgument("cellGraphOutput", std::string("."));
  parser.addStringArgument("subSteps", std::string("."));
  parser.addDoubleArgument("particleSpacing", 0.00390625);
  parser.addIntArgument("frames", 60);
  parser.parse(argv[1]);

  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation(domainRegion, parser.getInt("minLevel"),
                                  parser.getInt("maxLevel"));
  auto particles = simulation.particles();
  static int first = 0;
  {
    // Injector2 injector(particles, 0.0078125);
    Injector2 injector(particles, parser.getDouble("particleSpacing"));

    if (parser.getString("config") == "drop") {
      // TANK
      injector.setupBoxShape(
          BBox2d(Point2d(0.001, 0.001), Point2d(0.999, 0.2)));
      // DROP
      injector.setupBoxShape(BBox2d(Point2d(0.4, 0.4), Point2d(0.6, 0.6)));
    } else if (parser.getString("config") == "rightDam") {
      // DAM BREAK DIREITA
      injector.setupBoxShape(BBox2d(Point2d(0.8, 0.001), Point2d(0.999, 0.5)));
    }
  }
  simulation.buildTree();
  simulation.updateTree();
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  //  simulation.reseedParticles();
  for (int frames = 0; frames < parser.getInt("frames"); ++frames) {
    CodeProfiler::instance().reset();
    //    if (first > 3)
    //      return;
    std::cerr << "frame " << first++ << "particles " << particles->size()
              << std::endl;
    double dt = 0.01;
    // simulation loop
    size_t cflSteps = static_cast<size_t>(
        std::ceil(std::max(dt / (simulation.cfl() * .75), 1.)));
    double cflDT = dt / cflSteps;
    std::cerr << cflSteps << " substeps, " << cflDT << std::endl;
    TimeLogger timer;
    for (size_t i = 0; i < cflSteps; i++) {
      simulation.advectParticles(cflDT);
      timer.log("advect");
      simulation.updateTree();
      timer.log("update");
      //      simulation.markCells();
      //      timer.log("markCells");
      simulation.markFaces();
      timer.log("markFaces");
      try {
        simulation.transferFromParticlesToGrid();
        timer.log("part -> grid");
      } catch (std::string &e) {
        std::cerr << e << std::endl;
      } catch (const char *e) {
        std::cerr << e << std::endl;
      }
      simulation.computeExternalForces(cflDT);
      timer.log("gravity");
      simulation.applyBoundaryCondition();
      timer.log("boundary");
      try {
        simulation.solvePressure(cflDT);
        timer.log("pressure");
      } catch (std::string &e) {
        std::cerr << "PPPPP " << e << std::endl;
      } catch (const char *e) {
        std::cerr << "PPPPP " << e << std::endl;
      }
      simulation.reseedParticles();
      timer.log("reseed");
      simulation.transferFromGridToParticles();
      timer.log("grid -> part");

      std::ostringstream velFile;
      velFile << parser.getString("substepsOutput") << "/" << std::setfill('0')
              << std::setw(6) << first << i;
      IO::saveToPV(velFile.str().c_str(), particles);
      IO::saveCellGraph(velFile.str().c_str(), structure, 2);
      std::cerr << "===========================" << std::endl;
    }
    timer.report();
    CodeProfiler::instance().createSnapshot("profile");

    std::ostringstream particlesFile;
    particlesFile << parser.getString("particlesOutput") << "/"
                  << std::setfill('0') << std::setw(6) << first;
    IO::saveToPV(particlesFile.str().c_str(), particles);
    std::ostringstream cellGraphFile;
    particlesFile << parser.getString("cellGraphOutput") << "/"
                  << std::setfill('0') << std::setw(6) << first;
    IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
    std::cerr << "===========================" << std::endl;
  }
  return 0;
}
