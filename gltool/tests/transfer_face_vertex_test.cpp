#include <cell_graph2_model.h>

int main() {
  aergia::SceneApp<> app(800, 800, "", false);
  app.addViewport2D(0, 0, 800, 800);
  furoo::CellGraph2 graph(furoo::BBox2d::squareBox(1.),
                          [](const furoo::CellGraph2::Node &node) {
                            return node.id() % 3 && node.level() <= 5;
                          });
  size_t vertexFieldId = graph.addVertexField<double>(0.);
  size_t faceFieldId = graph.addFaceField<double>(0.);
  graph.iterateFaces([&](size_t id) {
    auto &faceField = *graph.field<double>(faceFieldId);
    faceField[id] =
        graph.faceCenterPosition(id).x() + graph.faceCenterPosition(id).y();
  });
  graph.transferFaceFieldToVertexField(faceFieldId,
                                       vertexFieldId,
                                       furoo::Definitions::Orientation::HORIZONTAL,
                                       nullptr);
  CellGraph2Model cm(&graph);
  auto vertexField = new FieldModel<double>(graph.field<double>(vertexFieldId));
  cm.selectVertexField(static_cast<int>(cm.addVertexField(vertexField)));
  auto faceField = new FieldModel<double>(graph.field<double>(faceFieldId));
  cm.selectFaceField(static_cast<int>(cm.addFaceField(faceField)));
  cm.updateVertexRegions();
  cm.updateFaceRegions();
  app.scene.add(&cm);
  app.run();
  return 0;
}