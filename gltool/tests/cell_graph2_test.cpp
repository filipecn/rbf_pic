#include "cell_graph2_model.h"

int main() {
  furoo::CellGraph2 graph(furoo::BBox2d::squareBox(1.f));
  aergia::SceneApp<> app(800, 800, "", false);
  app.addViewport2D(0, 0, 800, 800);
  CellGraphModel<2> model;
  model.set(&graph);
  model.showText = true;
  model.showFaceText = true;
  model.showCellText = true;
  model.updateCellsEdges();
  model.updateCellRegions();
  model.updateFacesEdges();
  model.updateFaceRegions();
  model.updateVerticesEdges();
  model.updateVertexRegions();
  app.scene.add(&model);
  enum class mode { EDIT = 0, NEIGHBORS = 1 };
  mode curMode = mode::EDIT;
  graph.addFaceField(0.,
                     "",
                     furoo::Definitions::MeshLocation::HORIZONTAL_FACE_CENTER);
  app.buttonCallback = [&](int button, int action, int modes) {
    if (action != GLFW_RELEASE)
      return;
    auto clickPoint =
        app.getCamera(0)->pickRay(app.viewports[0].getMouseNPos()).o.xy();
    auto cellId = graph.cellId(furoo::Point2d(clickPoint.x, clickPoint.y));
    if (cellId <= 0)
      return;
    model.selectedCells.clear();
    if (curMode == mode::EDIT) {
      if (button == GLFW_MOUSE_BUTTON_LEFT) {
        std::cerr << "refine(";
        graph.refine(cellId);
      } else {
        std::cerr << "coarse(";
        graph.coarse(cellId);
      }
      std::cerr << cellId << ")" << std::endl;
    }
    if (curMode == mode::NEIGHBORS) {
      graph.sampleFaceField(0,
                            furoo::Point2d(clickPoint.x, clickPoint.y),
                            furoo::Definitions::Orientation::HORIZONTAL);
    }
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
      if (key == GLFW_KEY_E)
        curMode = mode::EDIT;
      if (key == GLFW_KEY_N)
        curMode = mode::NEIGHBORS;
      if (key == GLFW_KEY_ESCAPE)
        app.exit();
      if (key == GLFW_KEY_S)
        furoo::IO::saveCellGraph("cellGraphFile", &graph);
      if (key == GLFW_KEY_L) {
        furoo::IO::loadCellGraph("cellGraphFile", &graph);
        graph.removeFields();
        graph.addFaceField(0.,
                           "",
                           furoo::Definitions::MeshLocation::HORIZONTAL_FACE_CENTER);

      }
      if (key == GLFW_KEY_P)
        graph.print();
      if(key == GLFW_KEY_T) {
        model.showFaceText = !model.showFaceText;
        model.showCellText = !model.showCellText;
      }
      model.updateCellsEdges();
    }
  };
  app.run();
  return 0;
}