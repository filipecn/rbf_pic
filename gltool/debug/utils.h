#ifndef FUROO_GLTOOL_UTILS_H
#define FUROO_GLTOOL_UTILS_H

#include <aergia/aergia.h>
#include <cell_graph2_model.h>
#include <iomanip>
#include <particle_system_model.h>
#include <string>

inline void flipImage(std::vector<unsigned char> &data, size_t w, size_t h) {
  for (uint x = 0; x < w; x++)
    for (uint y = 0; y < h / 2; y++)
      for (uint k = 0; k < 4; k++) {
        unsigned char tmp = data[4 * w * (h - 1 - y) + 4 * x + k];
        data[4 * w * (h - 1 - y) + 4 * x + k] = data[4 * w * y + 4 * x + k];
        data[4 * w * y + 4 * x + k] = tmp;
      }
}

inline size_t countFrameFiles(const std::string &path, std::string prefix = "",
                              std::string suffix = "") {
  size_t frame(5);
  std::string fileName;
  while (true) {
    fileName = furoo::concat(path, '/', prefix, std::setfill('0'), std::setw(6),
                             frame++, suffix);
    std::ifstream fp(fileName, std::istream::in);
    if (!fp.good())
      break;
  }
  return frame - 1;
}

template <int D> class SimulationModel : public aergia::SceneObject {
public:
  SimulationModel(furoo::CustomSimulation<D> *sim) : _sim(sim) {
    particleSystemModel.set(sim->particles(), particleRadius);
    cellGraph2Model.set(sim->cellGraph());
    cellGraphTopologyModel.set(sim->cellGraph());
  }
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    particleSystemModel.draw(camera, t);
    // cellGraph2Model.draw(camera, t);
    cellGraphTopologyModel.draw(camera, t);
  }
  void updateRegions() {
    cellGraph2Model.updateCellRegions();
    cellGraph2Model.updateFaceRegions();
    cellGraph2Model.updateVertexRegions();
  }
  void updateEdges() {
    cellGraph2Model.updateCellsEdges();
    cellGraph2Model.updateFacesEdges();
    cellGraph2Model.updateVerticesEdges();
    cellGraphTopologyModel.updateNodes();
    cellGraphTopologyModel.updateEdges();
  }
  void updateParticles() { particleSystemModel.update(); }
  void loadFrame(std::string inputDirectory, size_t frame) {
    auto end = furoo::concat(std::setfill('0'), std::setw(6), frame);
    furoo::IO::loadFromPV(furoo::concat(inputDirectory, "/p", end),
                          _sim->particles());
    updateParticles();
    furoo::IO::loadCellGraph(furoo::concat(inputDirectory, "/c", end),
                             _sim->cellGraph());
    updateEdges();
    furoo::IO::loadFields(furoo::concat(inputDirectory, "/f", end),
                          _sim->cellGraph());
    updateRegions();
    updateEdges();
    updateParticles();
  }
  ParticleSystemModel<D> particleSystemModel;
  CellGraphModel<D> cellGraph2Model;
  CellGraphTopologyModel<D> cellGraphTopologyModel;
  double particleRadius = 0.0015;
  furoo::CustomSimulation<D> *_sim = nullptr;
};

#endif // FUROO_UTILS_H