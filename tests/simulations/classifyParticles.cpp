#include <furoo.h>
//#include <iomanip>

int main(int argc, char **argv) {
  if (argc < 3)
    std::cout << "custom_pic <Dimensions> <simFolder> <n>\n";
  int D = atoi(argv[1]);
  int n = atoi(argv[3]);
  if (D == 2) {
    furoo::CustomSimulation<2> sim;
    sim.loadFromFile(furoo::concat(argv[2], "/configFile"));
    for (int frame = 0; frame < n; frame++) {
      auto filename = furoo::concat(argv[2], '/', std::setfill('0'),
                                    std::setw(6), frame, ".pvsl");
      std::ofstream fp(filename, std::ofstream::out);
      std::cerr << "FRAME" << frame << std::endl;
      sim.loadFrame(argv[2], static_cast<size_t>(atoi(argv[3])));
      auto structure = sim.cellGraph();
      auto particles = sim.particles();
      fp << particles->size() << std::endl;
      particles->iterateParticles([&](size_t id, furoo::Point2d p) {
        int cellId = structure->cellId(p);
        fp << id << " " << p.x() << " " << p.y() << " "
           << particles->getScalarProperty(0, id) << " "
           << particles->getScalarProperty(1, id) << " "
           << static_cast<int>((*structure->template field<unsigned char>(
                  sim.fieldId("cellSurfaceMaskField")))[cellId])
           << " "
           << dynamic_cast<furoo::CellGraph2 *>(structure)->node(cellId).level()
           << std::endl;
      });
      fp.close();
    }
  } else {
    furoo::CustomSimulation<3> sim;
    sim.loadFromFile(furoo::concat(argv[2], "/configFile"));
    for (int frame = 0; frame < n; frame++) {
      auto filename = furoo::concat(argv[2], '/', std::setfill('0'),
                                    std::setw(6), frame, ".pvsl");
      std::ofstream fp(filename, std::ofstream::out);
      std::cerr << "FRAME" << frame << std::endl;
      sim.loadFrame(argv[2], static_cast<size_t>(atoi(argv[3])));
      auto structure = sim.cellGraph();
      auto particles = sim.particles();
      fp << particles->size() << std::endl;
      particles->iterateParticles([&](size_t id, furoo::Point3d p) {
        int cellId = structure->cellId(p);
        fp << id << " " << p.x() << " " << p.y() << " " << p.z() << " "
           << particles->getScalarProperty(0, id) << " "
           << particles->getScalarProperty(1, id) << " "
           << particles->getScalarProperty(2, id) << " "
           << static_cast<int>((*structure->template field<unsigned char>(
                  sim.fieldId("cellSurfaceMaskField")))[cellId])
           << " "
           << dynamic_cast<furoo::CellGraph2 *>(structure)->node(cellId).level()
           << std::endl;
      });
      fp.close();
    }
  }
  return 0;
}
