#include <furoo.h>
//#include <iomanip>

int main(int argc, char **argv) {
  if (argc < 3)
    std::cout << "section <simFolder> <n>\n";
  int n = atoi(argv[2]);
  furoo::CustomSimulation<3> sim;
  sim.loadFromFile(furoo::concat(argv[1], "/configFile"));
  for (int frame = 0; frame < n; frame++) {
    auto c_filename = furoo::concat(argv[1], '/', std::setfill('0'),
                                    std::setw(6), frame, ".cquad");
    std::ofstream fp(c_filename, std::ofstream::out);
    std::cerr << "FRAME" << frame << std::endl;
    sim.loadFrame(argv[1], static_cast<size_t>(atoi(argv[2])));
    auto structure = dynamic_cast<furoo::CellGraph3 *>(sim.cellGraph());
    auto &materialField =
        *structure->template field<furoo::Definitions::Material>(
            sim.fieldId("cellMaterialField"));
    fp << structure->vertexCount() << std::endl;
    structure->iterateVertices([&](size_t vertexId) {
      /// vertexId positionX_i positionY_j
      auto &vertex = structure->vertex(vertexId);
      fp << vertex.id() << " " << vertex.position().x() << " "
         << vertex.position().y() << " " << vertex.position().z() << std::endl;
    });
    structure->iterateFaces([&](size_t faceId) {
      auto cells = structure->faceCells(faceId);
      if (cells.size() != 2)
        return;
      bool hasFluid = false;
      for (auto cell : cells)
        if (materialField[cell.id] == furoo::Definitions::Material::FLUID)
          hasFluid = true;
      if (hasFluid) {
        size_t minCell = 0;
        if (structure->node(cells[0].id).level() <
            structure->node(cells[1].id).level())
          minCell = 1;
        auto node = structure->node(cells[minCell].id);
        switch (cells[minCell].side) {
        case furoo::Definitions::Side::LEFT:
          fp << node.vertexAt(1) << " ";
          fp << node.vertexAt(5) << " ";
          fp << node.vertexAt(7) << " ";
          fp << node.vertexAt(3) << " ";
          break;
        case furoo::Definitions::Side::RIGHT:
          fp << node.vertexAt(0) << " ";
          fp << node.vertexAt(4) << " ";
          fp << node.vertexAt(6) << " ";
          fp << node.vertexAt(2) << " ";
          break;
        case furoo::Definitions::Side::BOTTOM:
          fp << node.vertexAt(2) << " ";
          fp << node.vertexAt(3) << " ";
          fp << node.vertexAt(7) << " ";
          fp << node.vertexAt(6) << " ";
          break;
        case furoo::Definitions::Side::TOP:
          fp << node.vertexAt(0) << " ";
          fp << node.vertexAt(1) << " ";
          fp << node.vertexAt(5) << " ";
          fp << node.vertexAt(4) << " ";
          break;
        case furoo::Definitions::Side::BACK:
          fp << node.vertexAt(4) << " ";
          fp << node.vertexAt(5) << " ";
          fp << node.vertexAt(7) << " ";
          fp << node.vertexAt(6) << " ";
          break;
        case furoo::Definitions::Side::FRONT:
          fp << node.vertexAt(0) << " ";
          fp << node.vertexAt(1) << " ";
          fp << node.vertexAt(3) << " ";
          fp << node.vertexAt(2) << " ";
          break;
        default:
          break;
        }
        fp << std::endl;
      }
    });
    fp.close();
  }
  return 0;
}
