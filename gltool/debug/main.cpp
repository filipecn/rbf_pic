#include "external/lodepng.h"
#include "ui.h"
#include "utils.h"
#include <furoo.h>

using namespace furoo;

template <int D> class SimVis : public aergia::SceneObject {
public:
  explicit SimVis() {
    _particles.reset(new ParticleSystem<double, D>(
        new PointZGrid<double, D>(16),
        new MLS<D>(Definitions::PolynomialType::QUADRATIC)));
    _particles->addScalarProperty(0.);
    _particles->addScalarProperty(0.);
    particleSystemModel.set(_particles.get(), 0.001);
    if (D == 2)
      _structure.reset(dynamic_cast<StructureInterface<D> *>(
          new CellGraph2(BBox2d::squareBox(1.))));
    else {
      _particles->addScalarProperty(0.);
      _structure.reset(dynamic_cast<StructureInterface<D> *>(
          new CellGraph3(BBox3d::squareBox(1.))));
    }
    cellGraph2Model.set(_structure.get());

    _screen.reset(new aergia::NanoGUIScreen());
    {
      auto *_gui = new nanogui::FormHelper(_screen.get());
      nanogui::ref<nanogui::Window> win =
          _gui->addWindow(Eigen::Vector2i(10, 10), "Settings");
      win->setLayout(new nanogui::GroupLayout());
      auto *b =
          new nanogui::Button(win, "Images output Dir", ENTYPO_ICON_ARCHIVE);
      b->setCallback([&]() { outputDirectory = nanogui::folder_dialog()[0]; });
      b = new nanogui::Button(win, "Simulation Directory", ENTYPO_ICON_ARCHIVE);
      b->setCallback([&]() {
        inputDirectory = nanogui::folder_dialog()[0];
        std::string fileName;
        // count particle files
        _slider.maxFrame = countFrameFiles(inputDirectory, "p", ".pv");
        _pCount->setCaption(concat("Particle Files Count: ", _slider.maxFrame));
        // count cell files
        size_t count = countFrameFiles(inputDirectory, "c");
        //        _slider.maxFrame = std::min(_slider.maxFrame, count);
        _cCount->setCaption(concat("Cell Files Count: ", count));
        // count particle files
        count = countFrameFiles(inputDirectory, "f");
        _fCount->setCaption(concat("Field Files Count: ", count));
        //        _slider.maxFrame = std::min(_slider.maxFrame, count);
        _slider.curFrame = 1;
        loadFrame(inputDirectory, 1);
        _screen->performLayout();
      });
      _pCount = new nanogui::Label(win, "Particle Files Count: 0");
      _cCount = new nanogui::Label(win, "CellGraph Files Count: 0");
      _fCount = new nanogui::Label(win, "Field Files Count: 0");
      _save = new nanogui::CheckBox(win, "save images");
      _playButton = new nanogui::Button(win, "", ENTYPO_ICON_CONTROLLER_PLAY);
      _playButton->setCallback([&]() {
        if (_playButton->icon() == ENTYPO_ICON_CONTROLLER_PLAY)
          _playButton->setIcon(ENTYPO_ICON_CONTROLLER_PAUS);
        else
          _playButton->setIcon(ENTYPO_ICON_CONTROLLER_PLAY);
      });
      _curFrameLabel = new nanogui::Label(win, "Frame: X");
      new nanogui::Label(win, "Cell Field");
      cellFields.set(win, [&](size_t id) {
        cellGraph2Model.selectCellField(id);
        cellGraph2Model.updateCellRegions();
      });
      new nanogui::Label(win, "Face Field");
      faceFields.set(win, [&](size_t id) {
        cellGraph2Model.selectFaceField(id);
        cellGraph2Model.updateFaceRegions();
      });
      new nanogui::Label(win, "Vertex Field");
      vertexFields.set(win, [&](size_t id) {
        cellGraph2Model.selectVertexField(id);
        cellGraph2Model.updateVertexRegions();
      });
    }
    {
      auto *_gui = new nanogui::FormHelper(_screen.get());
      nanogui::ref<nanogui::Window> win =
          _gui->addWindow(Eigen::Vector2i(270, 10), "Frames");
      win->setLayout(new nanogui::GroupLayout());
      auto *wid = new nanogui::Widget(win);
      wid->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
                                            nanogui::Alignment::Middle, 0, 2));
      _slider.set(wid);
      _slider.goTo(1);
    }
    _screen->setVisible(true);
    _screen->performLayout();
  }
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    static size_t curFrame = _slider.curFrame;
    if (_playButton->icon() == ENTYPO_ICON_CONTROLLER_PAUS &&
        _slider.curFrame < _slider.maxFrame) {
      _slider.goTo(_slider.curFrame + 1);
    }
    if (curFrame != _slider.curFrame) {
      curFrame = _slider.curFrame;
      _curFrameLabel->setCaption(concat("Frame ", _slider.curFrame));
      loadFrame(inputDirectory, _slider.curFrame);
      _screen->performLayout();
    }
    if (_particles->size()) {
      particleSystemModel.draw(camera, t);
      cellGraph2Model.draw(camera, t);
      if (_playButton->icon() == ENTYPO_ICON_CONTROLLER_PAUS &&
          _save->checked()) {
        frameCallback(concat(outputDirectory, "/", std::setfill('0'),
                             std::setw(6), _slider.curFrame, ".png"));
      }
    }
  }
  void nextFrame() { _slider.goTo(_slider.curFrame + 1); }
  void previousFrame() {
    if (_slider.curFrame)
      _slider.goTo(_slider.curFrame - 1);
  }
  void loadFrame(std::string inputDirectory, size_t frame) {
    auto end = concat(std::setfill('0'), std::setw(6), frame);
    furoo::IO::loadFromPV(concat(inputDirectory, "/p", end), _particles.get());
    particleSystemModel.update();
    return;
    furoo::IO::loadCellGraph(concat(inputDirectory, "/c", end),
                             _structure.get());
    furoo::IO::loadFields(concat(inputDirectory, "/f", end), _structure.get());
    cellGraph2Model.clearFields();
    connectStructureFieldsToModelAndUI<D>(_structure.get(), cellGraph2Model,
                                          cellFields, faceFields, vertexFields,
                                          cellPointFields);
    cellGraph2Model.updateCellRegions();
    cellGraph2Model.updateFaceRegions();
    cellGraph2Model.updateVertexRegions();
    cellGraph2Model.updateCellsEdges();
    cellGraph2Model.updateFacesEdges();
    cellGraph2Model.updateVerticesEdges();
  }
  std::function<void(const std::string &)> frameCallback;

private:
  ComboBox<size_t> cellFields, faceFields, vertexFields, cellPointFields;
  nanogui::CheckBox *_save;
  FrameSlider _slider;
  std::shared_ptr<aergia::NanoGUIScreen> _screen;
  nanogui::Label *_pCount, *_cCount, *_fCount, *_curFrameLabel;
  nanogui::Button *_playButton;
  std::string inputDirectory, outputDirectory;
  std::shared_ptr<ParticleSystem<double, D>> _particles;
  std::shared_ptr<StructureInterface<D>> _structure;
  ParticleSystemModel<D> particleSystemModel;
  CellGraphModel<D> cellGraph2Model;
};

int main() {
  aergia::SceneApp<> app(2000, 2000, "Furoo Debugger", false);
  app.addViewport(0, 0, 2000, 2000);
  // app.addViewport2D(0, 0, 2000, 2000);
  aergia::Quad q;
  ponos::Point2 point(0.999694, 0.311804);
  q.set(point - ponos::vec2(0.001), point + ponos::vec2(0.001));
  app.scene.add(&q);
  app.scene.add(new aergia::CartesianGrid(5));
  SimVis<3> simVis;
  app.scene.add(&simVis);
  app.scene.add(&q);
  simVis.frameCallback = [&](const std::string &path) {
    std::cerr << "saving frame image to " << path << std::endl;
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
  // auto region = BBox2d::squareBox();
  // app.getCamera<aergia::UserCamera2D>(0)->fit(Ponos::bbox2D(region));
  app.keyCallback = [&](int key, int coded, int action, int mods) {
    if (action != GLFW_RELEASE)
      return;
    if (key == GLFW_KEY_ESCAPE)
      app.exit();
    if (key == GLFW_KEY_LEFT)
      simVis.previousFrame();
    if (key == GLFW_KEY_RIGHT)
      simVis.nextFrame();
  };
  app.run();
  return 0;
}
