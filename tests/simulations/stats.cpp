#include <furoo.h>
#include <ctime>
//#include <iomanip>

int main(int argc, char **argv) {
  if (argc < 2)
    std::cout << "section <simFolder> <n>\n";
  int n = atoi(argv[2]);
  furoo::CustomSimulation<3> sim;
  sim.loadFromFile(furoo::concat(argv[1], "/configFile"));
  auto p_filename = furoo::concat(argv[1], "/stats");
  std::ofstream fp(p_filename, std::ofstream::out);
  for (int frame = n; frame <= n; frame++) {
    std::cerr << "FRAME" << frame << std::endl;
    sim.loadFrame(argv[1], frame);
    sim.setOutputPath("stats_output");
    sim.advanceFrame();
    // return 0;
    // sim.run("updateTree");
    // sim.run("reseedParticles");
    auto structure = dynamic_cast<furoo::CellGraph3 *>(sim.cellGraph());
    auto particles = sim.particles();
    auto &materialField =
        *structure->template field<furoo::Definitions::Material>(
            sim.fieldId("cellMaterialField"));
    auto &faceMaterialField =
        *structure->template field<furoo::Definitions::Material>(
            sim.fieldId("faceMaterialField"));
    size_t fluidCells = 0, rFluidCells = 0, nParticles = 0, rNParticles = 0;
    size_t divStencils = 0, gradStencils = 0, laplacianStencils = 0;
    structure->iterateCells([&](size_t cellId) {
      if (materialField[cellId] != furoo::Definitions::Material::FLUID)
        return;
      fluidCells++;
      auto level = structure->node(cellId).level();
      rFluidCells +=
          (1 << (6 - level)) * (1 << (6 - level)) * (1 << (6 - level));
      auto faces = structure->cellFaces(cellId);
      divStencils += faces.size();
      std::vector<size_t> neighbors;
      auto directNeighbors = structure->cellNeighbors(cellId);
      for (auto n : directNeighbors)
        if (n.id >= 0 &&
            materialField[n.id] == furoo::Definitions::Material::SOLID)
          neighbors.emplace_back(n.id);
      structure->iterateCellRing(cellId, [&](size_t nid) {
        if (materialField[nid] != furoo::Definitions::Material::SOLID)
          neighbors.emplace_back(nid);
      });
      laplacianStencils += neighbors.size();
    });
    size_t fluidFaces = 0;
    structure->iterateFaces([&](size_t faceId) {
      if (faceMaterialField[faceId] != furoo::Definitions::Material::FLUID)
        return;
      fluidFaces++;
      auto faceCells = structure->faceCells(faceId);
      // FIXME: what about domain faces???
      if (structure->cellRegion(faceCells[0].id).size() ==
          structure->cellRegion(faceCells[1].id).size()) {
        for (auto cell : faceCells) {
          gradStencils++;
        }
      } else {
        std::set<size_t> cellsSet;
        cellsSet.insert(faceCells[0].id);
        cellsSet.insert(faceCells[1].id);
        for (auto cell : faceCells) {
          auto neighborCells = structure->cellNeighbors(cell.id);
          for (auto nCell : neighborCells) {
            if (nCell.id >= 0)
              cellsSet.insert(nCell.id);
          }
        }
        for (auto cell : cellsSet) {
          if (materialField[cell] != furoo::Definitions::Material::SOLID) {
            gradStencils++;
          }
        }
      }
    });
    nParticles = particles->size();
    rNParticles = rFluidCells * 8;
    std::cerr << "FRAME " << frame << std::endl;
    std::cerr << nParticles << " ";
    std::cerr << rNParticles << " ";
    std::cerr << fluidCells << " ";
    std::cerr << rFluidCells << " ";
    std::cerr << (1.f * divStencils) / fluidCells << " ";
    std::cerr << (1.f * laplacianStencils) / fluidCells << " ";
    std::cerr << (1.f * gradStencils) / fluidFaces << " ";
    std::cerr << std::endl;
    fp.close();
  }
  return 0;
}
