#include "external/lodepng.h"
#include "ui.h"
#include "utils.h"

using namespace furoo;

int main(int argc, char **argv) {
  aergia::SceneApp<> app(2000, 2000, "CellGraph PIC Debugger", false);
  app.addViewport2D(0, 0, 2000, 2000);
  UI<2> ui;
  // UI CELL FIELDS
  aergia::Quad q;
  ponos::Point2 point(0.148438, 0.59375);
  q.set(point - ponos::vec2(0.001), point + ponos::vec2(0.001));
  app.scene.add(&q);
  ui.saveFrameImageCallback = [&](std::string path) {
    size_t w = 0, h = 0;
    std::vector<unsigned char> data;
    app.viewports[0].renderer->currentPixels(data, w, h);
    flipImage(data, w, h);
    unsigned error =
        lodepng::encode(path, &data[0], static_cast<unsigned int>(w),
                        static_cast<unsigned int>(h));
    if (error)
      std::cout << "encoder error " << error << ": "
                << lodepng_error_text(error) << std::endl;
  };
  try {
    ui.init();
    if (argc == 4)
      ui.loadConfigFile(argv[1], argv[2], atoi(argv[3]));
  } catch (std::string &e) {
    std::cerr << e << std::endl;
    return -1;
  }
  app.scene.add(&ui);
  app.scene.add(new aergia::CartesianGrid(5));
  app.dropCallback = [&](size_t count, const char **filenames) {
    ui.dropPath(filenames[0]);
  };
  app.keyCallback = [&](int key, int scancode, int action, int modifiers) {
    if (action == GLFW_PRESS) {
      if (key == GLFW_KEY_X && modifiers & GLFW_MOD_SHIFT)
        ui.movePlane(0, 1);
      else if (key == GLFW_KEY_X)
        ui.movePlane(0, -1);
      else if (key == GLFW_KEY_Y && modifiers & GLFW_MOD_SHIFT)
        ui.movePlane(1, 1);
      else if (key == GLFW_KEY_Y)
        ui.movePlane(1, -1);
      else if (key == GLFW_KEY_Z && modifiers & GLFW_MOD_SHIFT)
        ui.movePlane(2, 1);
      else if (key == GLFW_KEY_Z)
        ui.movePlane(2, -1);
      else if (key == GLFW_KEY_ESCAPE)
        app.exit();
    }
  };
  auto region = BBox2d::squareBox();
  app.getCamera<aergia::UserCamera2D>(0)->fit(Ponos::bbox2D(region));
  try {
    app.run();
  } catch (std::string &e) {
    std::cerr << e << std::endl;
    return -1;
  }
  return 0;
}
