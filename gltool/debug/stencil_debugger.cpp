#include <graphic_element_set.h>
#include <nanogui_utils.h>

using namespace furoo;

class StencilDebugger : public aergia::SceneObject {
 public:
  StencilDebugger() {
    _screen.reset(new aergia::NanoGUIScreen());
    auto *_gui = new nanogui::FormHelper(_screen.get());
    nanogui::ref<nanogui::Window>
        win = _gui->addWindow(Eigen::Vector2i(10, 10), "Settings");
    win->setLayout(new nanogui::GroupLayout());
    auto *b = new nanogui::Button(win, "reset");
    b->setCallback([&]() {
      points.clear();
    });
    b = new nanogui::Button(win, "run");
    b->setCallback([&]() {
      results.clear();
      std::vector<Point2d> stencilPoints = points;
      stencilPoints.emplace_back(center);
      std::vector<double> stencilValues = values;
      stencilValues.emplace_back(centerValue);
      for (auto rbf : rbfs)
        results.emplace_back(rbf->gradientAt(center,
                                             0,
                                             stencilPoints,
                                             stencilValues));
    });
    functionsBox.set(win, [&](std::string id) {
      values.clear();
      for (auto p : points)
        values.emplace_back(functions[id](p));
      centerValue = functions[id](center);
    });
    _screen->setVisible(true);
    _screen->performLayout();
    text.reset(new aergia::TextRenderer(FONT_PATH));
  }
  void update() {
    for (auto it : functions)
      functionsBox.add(it.first, it.first);
    _screen->performLayout();
  }
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (!points.empty())
      pointSet->draw(camera, t);
    text->setCamera(camera);
    text->textSize = 0.0008;
    if (values.size() == points.size())
      for (size_t i = 0; i < points.size(); i++)
        text->at(ponos::Point3(points[i].x(), points[i].y(), 0)) << values[i];
    text->at(ponos::Point3(center.x(), center.y(), 0)) << centerValue;
    q.draw(camera, t);
    for (size_t i = 0; i < results.size(); i++) {
      std::ostringstream s;
      s << results[i];
      text->render(s.str(), 0.1, 0.1);
    }
  }
  void click(ponos::Point2 clickPoint, int button) {
    if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
      q.set(clickPoint - ponos::vec2(0.005), clickPoint + ponos::vec2(0.005));
      center = Point2d(clickPoint.x, clickPoint.y);
      return;
    }
    if (points.empty()) {
      auto cubeMesh = ponos::RawMeshes::cube();
      pointSet = std::make_shared<GraphicElementSet>(cubeMesh.get());
    }
    points.emplace_back(clickPoint[0], clickPoint[1]);
    pointSet->resize(points.size());
    for (size_t i = 0; i < points.size(); i++) {
      auto p = pointSet->transformBuffer(i + 1);
      p[0] = points[i][0];
      p[1] = points[i][1];
      p[2] = 0.f;
      p[3] = 0.005;
      aergia::Color color = aergia::COLOR_BLACK;
      auto c = pointSet->colorBuffer(i + 1);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = pointSet->activeBuffer(i + 1);
      a[0] = 1;
    }
  }

  std::shared_ptr<aergia::TextRenderer> text;
  aergia::Quad q;
  std::shared_ptr<aergia::NanoGUIScreen> _screen;
  std::shared_ptr<GraphicElementSet> pointSet;
  std::vector<Point2d> points;
  Point2d center;
  ComboBox<std::string> functionsBox;
  std::map<std::string, std::function<double(Point2d)>> functions;
  std::vector<double> values;
  double centerValue;
  std::vector<double> results;
  std::vector<std::shared_ptr<DifferentialRBF<2>>> rbfs;
};

int main() {
  aergia::SceneApp<> app(2000, 2000, "Stencil Debugger", false);
  app.addViewport2D(0, 0, 2000, 2000);
  StencilDebugger debugger;
  debugger.rbfs.emplace_back(new DifferentialRBF<2>(new GaussianKernel<Point2d>(
      0.1)));
  debugger.functions["linearX"] = [](Point2d p) { return p.x(); };
  debugger.update();
  app.scene.add(&debugger);
  app.scene.add(new aergia::CartesianGrid(5));
  app.getCamera<aergia::UserCamera2D>(0)->fit(ponos::BBox2D::unitBox());
  app.buttonCallback = [&](int button, int action, int modes) {
    if (action != GLFW_RELEASE)
      return;
    if (button == GLFW_MOUSE_BUTTON_LEFT)
      return;
    auto clickPoint =
        app.getCamera(0)->pickRay(app.viewports[0].getMouseNPos()).o.xy();
    debugger.click(clickPoint, button);
  };
  app.run();
  return 0;
}