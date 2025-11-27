#include <furoo.h>
//#include <iomanip>

bool intersects(int res, int s, double a, double b) {
  if (s < 0)
    return false;
  double step = 1. / (1 << res);
  return a <= s * step + 0.5 * step && b <= s * step + 0.5 * step;
}

int main(int argc, char **argv) {
  if (argc < 7)
    std::cout << "section <simFolder> <n> <res> <x> <y> <z>\n";
  int n = atoi(argv[3]);
  int res = atoi(argv[4]);
  int x = atoi(argv[5]);
  int y = atoi(argv[6]);
  int z = atoi(argv[7]);
  furoo::CustomSimulation<3> sim;
  sim.loadFromFile(furoo::concat(argv[2], "/configFile"));
  for (int frame = 0; frame < n; frame++) {
    auto p_filename = furoo::concat(argv[2], '/', std::setfill('0'),
                                    std::setw(6), frame, ".psec");
    auto c_filename = furoo::concat(argv[2], '/', std::setfill('0'),
                                    std::setw(6), frame, ".csec");
    std::ofstream fp(p_filename, std::ofstream::out);
    std::ofstream fpc(c_filename, std::ofstream::out);
    std::cerr << "FRAME" << frame << std::endl;
    sim.loadFrame(argv[2], static_cast<size_t>(atoi(argv[3])));
    auto structure = dynamic_cast<furoo::CellGraph3 *>(sim.cellGraph());
    auto particles = sim.particles();
    auto &materialField =
        structure->template field<furoo::Definitions::Material>(
            sim.fielId("cellMaterialField"));
    structure->iterateCells([&](size_t cellId) {
      if (materialField[cellId] != furoo::Definitions::Material::FLUID)
        return;
      fpc << structure->node(cellId) auto region =
          structure->cellRegion(cellId);
      auto level = structure->node(cellId).level();
      if (intersects(res, x, region.lower().x(), region.upper().x()) ||
          intersects(res, y, region.lower().y(), region.upper().y()) ||
          intersects(res, z, region.lower().z(), region.upper().z()))
        particles->iterateParticles(
            region, level, [&](size_t id, furoo::Point3d p) {
              int cellId = structure->cellId(p);
              fp << id << " " << p.x() << " " << p.y() << " " << p.z() << " "
                 << particles->getScalarProperty(0, id) << " "
                 << particles->getScalarProperty(1, id) << " "
                 << particles->getScalarProperty(2, id) << " " << level
                 << std::endl;
            });
    });
    fp.close();
    fpc.close();
  }
  return 0;
}
