#include "cell_centered_graph2_model.h"
#include "particle_system_model.h"
#include "test_utils.h"
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
  parser.addStringArgument("substepsOutput", std::string("."));
  parser.addDoubleArgument("particleSpacing", 0.00390625);
  parser.addIntArgument("frames", 60);

  if (argc > 1)
    parser.parse(argv[1]);
  auto domainRegion = BBox2d::squareBox();
  AdaptiveCenteredPic2 simulation
      (domainRegion, parser.getInt("minLevel"), parser.getInt("maxLevel"));
  auto particles = simulation.particles();
  //  static int first = 0;
  {
    // Injector2 injector(particles, 0.0078125);
    Injector2 injector(particles, 0.00390625);
    //    injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.2,
    //                                                                 0.5)));
    //    injector.setupBoxShape(BBox2d(Point2d(0.45, 0.001), Point2d(0.55,
    //    0.2)));
    // TANK
//    injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.999,
//    0.2)));
    // DROP
    //    injector.setupBoxShape(BBox2d(Point2d(0.4, 0.4),
    //                                  Point2d(0.6, 0.6)));

    // DAM BREAK DIREITA
    injector.setupBoxShape(BBox2d(Point2d(0.001, 0.001), Point2d(0.2, 0.5)));
//    injector.setupCircleShape(Point2d(0.5, 0.4921875), 0.1);
    //    injector.setupBoxShape(BBox2d(Point2d(0.3, 0.52), Point2d(0.7,
    //    0.65))); particles->iterateParticles([&](size_t id, Point2d pos) {
    //      UNUSED_VARIABLE(pos);
    //      particles->setScalarProperty(1, id, -2 * 0.4);
    //    });
    //    injector.setupBoxShape(BBox2d(Point2d(0.3, 0.361), Point2d(0.7,
    //    0.481))); particles->iterateParticles([&](size_t id, Point2d pos) {
    //      UNUSED_VARIABLE(pos);
    //      double vel = particles->getScalarProperty(1, id);
    //      particles->setScalarProperty(1, id, vel + 0.4);
    //    });

    //    injector.loadFromFile("particles/p000021");
    //    first = 21;
  }
  simulation.buildTree();
//  simulation.resetTree(5);
  simulation.updateTree();
  //  simulation.reseedParticles();
//  simulation.loadFromFile("cellGraph/000053.cg");
  aergia::SceneApp<> app(800, 800, "", false);
  app.init();
  app.addViewport2D(0, 0, 800, 800);
  app.getCamera<aergia::UserCamera2D>(0)->fit(Ponos::bbox2D(domainRegion));

  ParticleSystemModel2 particleSystemModel(particles, 0.0015);
  auto structure = dynamic_cast<CellGraph2 *>(simulation.structure());
  IO::saveCellGraph("cellGraph.cg", structure, 2);

  CellCenteredGraph2Model cellGraph2Model(
      dynamic_cast<CellGraph2 *>(structure));

  app.scene.add(&particleSystemModel);
  app.scene.add(&cellGraph2Model);
  static int step = 9;
  app.keyCallback = [&](int k, int a, int t1, int t2) {
    UNUSED_VARIABLE(t1);
    UNUSED_VARIABLE(t2);
    if (a == GLFW_RELEASE) {
      if (k == GLFW_KEY_D) {
        particleSystemModel.debugMode = !particleSystemModel.debugMode;
        cellGraph2Model.debugMode = !cellGraph2Model.debugMode;
        return;
      } else if (k == GLFW_KEY_SPACE) {
        TimeLogger timer;
        double dt = 1 / 60.;
        size_t cflSteps = static_cast<size_t>(
            std::ceil(std::max(dt / (simulation.cfl() * .75), 1.)));
        double cflDT = dt / cflSteps;
        dt = cflDT;
        switch (step) {
          case 0:simulation.updateTree();
            timer.log("update");
            break;
          case 1:
            //          simulation.markCells();
            timer.log("markCells");
            break;
          case 2:simulation.markFaces();
            timer.log("markFaces");
            break;
          case 3:std::cerr << "PArticles to grid\n";
            try {
              simulation.transferFromParticlesToGrid();
              timer.log("part -> grid");
            } catch (std::string e) {

              std::cerr << e << std::endl;
              exit(-1);
            }
            break;
          case 4:simulation.computeExternalForces(dt);
            timer.log("gravity");
            break;
          case 5:simulation.applyBoundaryCondition();
            timer.log("boundary");
            break;
          case 6:
            try {
              simulation.solvePressure(dt);
              timer.log("pressure");
              cellGraph2Model.updatePressure();
            } catch (const char *e) {
              std::cerr << e << std::endl;
              THROW(false, "Solve pressure error!");
            }
            break;
          case 7:simulation.reseedParticles();
            timer.log("reseed");
            break;
          case 8:simulation.transferFromGridToParticles();
            timer.log("grid -> part");
            break;
          case 9:simulation.advectParticles(dt);
            timer.log("advect");
            break;
          default:step = -1;
        }
        timer.report();
        step++;
      }
    }
  };
  size_t frame = 0;
  app.renderCallback = [&]() {
//    if (frame >= 13)
//      return;
    CodeProfiler::instance().reset();
    std::cerr << "frame " << frame << "particles " << particles->size()
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
      {
        std::ostringstream particlesFile;
        particlesFile << parser.getString("particlesOutput") << "/afterAdvect";
        IO::saveToPV(particlesFile.str().c_str(), particles);
        std::ostringstream cellGraphFile;
        cellGraphFile << parser.getString("cellGraphOutput") << "/afterAdvect";
        IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
      }
      timer.log("advect");
      try {
        simulation.updateTree();
        timer.log("update");
      } catch (std::string e) {
        std::cerr << e;
        exit(12);
      }
      {
        std::ostringstream particlesFile;
        particlesFile << parser.getString("particlesOutput") << "/afterUpdate";
        IO::saveToPV(particlesFile.str().c_str(), particles);
        std::ostringstream cellGraphFile;
        cellGraphFile << parser.getString("cellGraphOutput") << "/afterUpdate";
        IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
      }
      //      simulation.markCells();
      //      timer.log("markCells");
      simulation.markFaces();
      timer.log("markFaces");
      try {
        simulation.transferFromParticlesToGrid();
        timer.log("part -> grid");
      } catch (std::string &e) {
        std::ostringstream particlesFile;
        particlesFile << parser.getString("particlesOutput") << "/error"
                      << frame;
        IO::saveToPV(particlesFile.str().c_str(), particles);
        std::ostringstream cellGraphFile;
        cellGraphFile << parser.getString("cellGraphOutput") << "/error"
                      << frame;
        IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
        std::cerr << e << std::endl;
        exit(13);
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
        exit(14);
      }
      simulation.reseedParticles();
      timer.log("reseed");
      try {
        simulation.transferFromGridToParticles();
        timer.log("grid -> part");
      } catch (std::string e) {
        std::ostringstream particlesFile;
        particlesFile << parser.getString("particlesOutput") << "/errorgp"
                      << frame;
        IO::saveToPV(particlesFile.str().c_str(), particles);
        std::ostringstream cellGraphFile;
        cellGraphFile << parser.getString("cellGraphOutput") << "/errorgp"
                      << frame;
        IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
        std::cerr << e;
        exit(15);
      }

      std::ostringstream velFile;
      velFile << parser.getString("substepsOutput") << "/" << std::setfill('0')
              << std::setw(6) << frame << "." << i;
      IO::saveToPV(velFile.str().c_str(), particles);
      IO::saveCellGraph(velFile.str().c_str(), structure, 2);
      std::cerr << "===========================" << std::endl;
    }
    timer.report();
    CodeProfiler::instance().createSnapshot("profile");

    std::ostringstream particlesFile;
    particlesFile << parser.getString("particlesOutput") << "/"
                  << std::setfill('0') << std::setw(6) << frame;
    IO::saveToPV(particlesFile.str().c_str(), particles);
    std::ostringstream cellGraphFile;
    cellGraphFile << parser.getString("cellGraphOutput") << "/"
                  << std::setfill('0') << std::setw(6) << frame;
    IO::saveCellGraph(cellGraphFile.str().c_str(), structure, 2);
    std::cerr << "===========================" << std::endl;
    saveFramebuffer(cellGraphFile.str().c_str(), 800, 800);
    frame++;
  };
  app.run();
  return 0;
}
