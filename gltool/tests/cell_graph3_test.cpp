#include "cell_graph2_model.h"

int main() {
  furoo::CellGraph3 graph(furoo::BBox3d::squareBox(1.f));
  aergia::SceneApp<> app(1500, 1500, "");
  CellGraphModel<3> model;
  model.set(&graph);
  graph.refine(1);
  graph.coarse(2);
  model.showText = true;
  model.showCellText = true;
  model.showFaceText = true;
  model.showVertexText = true;
  model.updateCellsEdges();
  model.updateCellRegions();
  model.updateFacesEdges();
  model.updateFaceRegions();
  model.updateVerticesEdges();
  model.updateVertexRegions();
  app.scene.add(&model);
  app.buttonCallback = [&](int button, int action, int modes) {
    if (action != GLFW_RELEASE)
      return;
    auto ray = app.getCamera(0)->pickRay(app.viewports[0].getMouseNPos());
    int cellId = -1;
    graph.iterateCells([&](size_t i) {
      auto center = graph.cellCenterPosition(i);
      if (ponos::sphere_ray_intersection(ponos::Sphere(ponos::Point3(center.x(),
                                                                     center.y(),
                                                                     center.z()),
                                                       0.01), ray))
        cellId = i;
    });
    if (cellId <= 0)
      return;
    model.selectedCells.clear();
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
      std::cerr << "refine(";
      graph.refine(cellId);
    } else {
      std::cerr << "coarse(";
      graph.coarse(cellId);
    }
    std::cerr << cellId << ")" << std::endl;
    model.updateCellsEdges();
    model.updateCellRegions();
    model.updateFacesEdges();
    model.updateFaceRegions();
    model.updateVerticesEdges();
    model.updateVertexRegions();
  };
  app.keyCallback = [&](int key, int scancode, int action, int mods) {
    if (action == GLFW_RELEASE) {
      if (key == GLFW_KEY_W)
        model.textSize += 0.0005;
      if (key == GLFW_KEY_Q)
        model.textSize -= 0.0005;
      if(key == GLFW_KEY_V)
        model.showVertexText = !model.showVertexText;
      if(key == GLFW_KEY_F)
        model.showFaceText = !model.showFaceText;
      if(key == GLFW_KEY_C)
        model.showCellText = !model.showCellText;
      if (key == GLFW_KEY_ESCAPE)
        app.exit();
    }
  };
  app.run();
  return 0;
}