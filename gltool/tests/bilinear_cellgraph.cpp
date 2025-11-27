#include <cell_graph2_model.h>

int main() {
  aergia::SceneApp<> app(800, 800, "", false);
  app.addViewport2D(0, 0, 800, 800);
  furoo::CellGraph2 graph(furoo::BBox2d::squareBox(1.),
                          [](const furoo::CellGraph2::Node &node) {
                            return node.id() % 4 && node.level() <= 6;
                          });
  size_t fieldId = graph.addVertexField<double>(0.);
  auto &field = *graph.field<double>(fieldId);
  graph.iterateVertices([&](size_t id) {
    field[id] = graph.vertexPosition(id).x() + graph.vertexPosition(id).y();
  });
  CellGraph2Model cm(&graph);
  auto vertexField = new FieldModel<double>();
  vertexField->field = &field;
  size_t vertexFieldId = cm.addVertexField(vertexField);
  cm.selectVertexField(vertexFieldId);
//  cm.updateVertexRegions();
  // sample points
  GraphicElementSet samplesSet(ponos::create_icosphere_mesh(
      ponos::Point3(0, 0, 0), 1.f, 3, false, false));
  size_t w = 100, h = 100;
  samplesSet.resize((w + 1) * (h + 1));
  aergia::ColorPalette palette = aergia::HEAT_MATLAB_PALETTE;
  size_t id = 0;
  for (int i = 0; i <= w; i++)
    for (int j = 0; j <= h; j++) {
      furoo::Point2d pos(i * (1.f / w), j * (1.f / h));
      auto p = samplesSet.transformBuffer(id);
      p[0] = pos[0]; // x
      p[1] = pos[1]; // y
      p[2] = 0.005f;    // z
      p[3] = 0.004;
      aergia::Color color = aergia::COLOR_BLACK;
      color = palette(1.f - ponos::linearStep(
          graph.sampleVertexField(fieldId, pos), 0, 1),
                      1.f);
      auto c = samplesSet.colorBuffer(id);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = samplesSet.activeBuffer(id);
      a[0] = 1;
      id++;
    }
  app.scene.add(&cm);
  app.scene.add(&samplesSet);
  app.run();
  return 0;
}