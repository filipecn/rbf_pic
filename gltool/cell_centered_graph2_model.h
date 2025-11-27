#ifndef FUROO_CELL_CENTERED_GRAPH2_MODEL_H
#define FUROO_CELL_CENTERED_GRAPH2_MODEL_H

#include <furoo.h>
#include <aergia/aergia.h>

class CellCenteredGraph2Model : public aergia::SceneObject {
 public:
  explicit CellCenteredGraph2Model(furoo::CellGraph2 *graph);
  void draw() override;
  void updatePressure();

  double dt;
  bool drawGraph;
  bool debugMode;
  std::vector<size_t> selectedCells;
  aergia::Color nodeColor;
  aergia::Color edgeColor;
  aergia::Color regionColor;
 private:
  double _minPressure, _maxPressure;
  furoo::CellGraph2 *_graph;
};

#endif //FUROO_CELL_CENTERED_GRAPH2_MODEL_H
