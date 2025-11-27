#include "../cell_graph2_model.h"

int main() {
  furoo::CellGraph2 graph(furoo::BBox2d::squareBox(1.),
                          [](const furoo::CellGraph2::Node &node) {
                            return node.level() < 3;
                          });
  size_t cid = graph.addCellField<double>(0.);
  graph.iterateCells([&](size_t i) {
    auto &field = *graph.field<double>(cid);
    field[i] = graph.cellCenterPosition(i).x();
  });
  size_t fid = graph.addFaceField<double>(0.);
  graph.iterateFaces([&](size_t i) {
    auto &field = *graph.field<double>(fid);
    field[i] = graph.faceCenterPosition(i).x();
  });
  size_t vid = graph.addVertexField<double>(0.);
  graph.iterateVertices([&](size_t i) {
    auto &field = *graph.field<double>(vid);
    field[i] = graph.vertexPosition(i).x();
  });
  aergia::SceneApp<> app(800, 800, "", false);
  app.addViewport2D(0, 0, 800, 800);
  CellGraph2Model model(&graph);
  auto cfield = new FieldModel<double>(graph.field<double>(cid));
  model.addCellField(cfield);
  auto ffield = new FieldModel<double>(graph.field<double>(fid));
  model.addFaceField(ffield);
  auto vfield = new FieldModel<double>(graph.field<double>(vid));
  model.addVertexField(vfield);
  model.selectCellField(0);
  model.selectFaceField(0);
  model.selectVertexField(0);
  model.showText = true;
  model.updateCellRegions();
  model.updateCellsEdges();
  model.updateFaceRegions();
  model.updateFacesEdges();
  model.updateVertexRegions();
  model.updateVerticesEdges();
  app.scene.add(&model);
  app.run();
  return 0;
}