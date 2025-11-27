#ifndef FUROO_NANOGUI_UTILS_H
#define FUROO_NANOGUI_UTILS_H

#include "utils.h"
#include <aergia/aergia.h>
#include <cell_graph2_model.h>
#include <particle_system_model.h>

template <typename T> struct ComboBox {
  void set(nanogui::Widget *widget,
           const std::function<void(T)> &callback = nullptr) {
    _data = new nanogui::ComboBox(widget, _items);
    if (callback) {
      _valueCallback = callback;
      _callback = [&](int id) { _valueCallback(_values[id]); };
      _data->setCallback(_callback);
    }
  }
  void clear() {
    _items.clear();
    _values.clear();
  }
  void add(std::string item, T value) {
    _items.emplace_back(item);
    _values.emplace_back(value);
    if (_data)
      _data->setItems(_items);
  }
  void copy(const ComboBox<T> &other) {
    clear();
    for (size_t i = 0; i < other._items.size(); i++)
      add(other._items[i], other._values[i]);
  }
  T selectedItem() { return _values[std::max(0, _data->selectedIndex())]; }

private:
  std::function<void(T)> _valueCallback;
  std::function<void(int)> _callback;
  std::vector<std::string> _items;
  std::vector<T> _values;
  nanogui::ComboBox *_data = nullptr;
};

struct FrameSlider {
  void set(nanogui::Widget *widget) {
    _slider = new nanogui::Slider(widget);
    _slider->setCallback([&](float f) {
      curFrame = static_cast<size_t>(f * maxFrame);
      if (callback)
        callback(curFrame);
    });
    _slider->setFixedWidth(400);
  }
  void goTo(size_t frame) {
    curFrame = frame;
    _slider->setValue(static_cast<float>(curFrame) /
                      static_cast<float>(maxFrame));
  }
  size_t curFrame = 0;
  size_t maxFrame = 0;
  std::function<void(size_t)> callback = nullptr;

private:
  nanogui::Slider *_slider = nullptr;
};

template <int D> struct TestFunctionPanel {
  TestFunctionPanel(
      nanogui::FormHelper *gui, nanogui::Widget *win,
      std::function<void(std::function<double(furoo::Point<double, D>)>)>
          cellCallback,
      std::function<void(std::function<double(furoo::Point<double, D>)>)>
          particleCallback,
      std::function<void(std::function<double(furoo::Point<double, D>)>)>
          functionCallback)
      : _cellCallback(cellCallback), _particleCallback(particleCallback),
        _functionCallback(functionCallback) {
    // UI
    functions.set(win);
    _testFunction = [](furoo::Point<double, D> p) -> double {
      return std::cos(10 * p.x()) * std::sin(20 * p.y());
      return (1. / 200.) * std::sin(10. * p.x()) * std::exp(10. * p.x());
    };
    functions.add("gradX", [](furoo::Point<double, D> p) -> double {
      return -10. * std::sin(10. * p.x()) * std::sin(20. * p.y());
    });
    functions.add("gradY", [](furoo::Point<double, D> p) -> double {
      return 20. * std::cos(10. * p.x()) * std::cos(20. * p.y());
    });
    functions.add("laplacian", [](furoo::Point<double, D> p) -> double {
      return std::exp(10. * p.x()) * std::cos(10. * p.x());
    });
    auto b = new nanogui::Button(win, "Apply To Cells");
    b->setCallback([&]() { _cellCallback(functions.selectedItem()); });
    b = new nanogui::Button(win, "Apply To Particles");
    b->setCallback([&]() { _particleCallback(functions.selectedItem()); });
    b = new nanogui::Button(win, "Run GradX");
    b->setCallback([&]() { runGradXCallback(); });
    b = new nanogui::Button(win, "Run GradY");
    b->setCallback([&]() { runGradYCallback(); });
    b = new nanogui::Button(win, "Run Laplacian");
    b->setCallback([&]() { runLaplacianCallback(); });
    b = new nanogui::Button(win, "load test function");
    b->setCallback([&]() { _functionCallback(_testFunction); });
  }
  ComboBox<std::function<double(furoo::Point<double, D>)>> functions;
  std::function<void()> runLaplacianCallback;
  std::function<void()> runGradXCallback;
  std::function<void()> runGradYCallback;

private:
  std::function<double(furoo::Point<double, D>)> _testFunction;
  std::function<void(std::function<double(furoo::Point<double, D>)>)>
      _cellCallback;
  std::function<void(std::function<double(furoo::Point<double, D>)>)>
      _particleCallback;
  std::function<void(std::function<double(furoo::Point<double, D>)>)>
      _functionCallback;
};

template <int D> struct VisPanel {
  VisPanel(SimulationModel<D> *model) : _model(model) {}
  void set(nanogui::FormHelper *gui, nanogui::Widget *win) {
    textSize = gui->addVariable("Text size", _textSize);
    useCutPlanes = new nanogui::CheckBox(win, "use cut planes");
    for (int d = 0; d < D; d++)
      cutPlanes[d] = new nanogui::CheckBox(win, furoo::concat("Plane", d));
    showTopology = new nanogui::CheckBox(win, "show topology");
    showGraph = new nanogui::CheckBox(win, "show cellgraph");
    showParticles = new nanogui::CheckBox(win, "show particles");
    showSolids = new nanogui::CheckBox(win, "show solids");
    showGraph->setChecked(true);
    showParticles->setChecked(true);
    particleVelocityScale =
        gui->addVariable("Velocity scale", _particleVelocityScale);
    particleFields.set(win, [&](size_t id) {
      _model->particleSystemModel.selectPropertyGroup(id);
    });
    showParticleVelocities = new nanogui::CheckBox(win, "show vectors");
    new nanogui::Label(win, "Cell Field");
    showCellText = new nanogui::CheckBox(win, "values in text");
    showCellFiels = new nanogui::CheckBox(win, "draw field");
    cellPointFields.set(win, [&](size_t id) {
      _model->cellGraph2Model.selectCellPointField(id);
      _model->cellGraph2Model.updateCellRegions();
    });
    cellFields.set(win, [&](size_t id) {
      _model->cellGraph2Model.selectCellField(id);
      _model->cellGraph2Model.updateCellRegions();
    });
    showCellVectors = new nanogui::CheckBox(win, "show vector field");
    cellVectorScale = gui->addVariable("Vector scale", _cellVectorScale);
    cellVectorFields[0].set(win, [&](size_t id) {
      _model->cellGraph2Model.selectCellVectorField(id, 0);
      _model->cellGraph2Model.updateCellsVector();
    });
    cellVectorFields[1].set(win, [&](size_t id) {
      _model->cellGraph2Model.selectCellVectorField(id, 1);
      _model->cellGraph2Model.updateCellsVector();
    });
    if (D == 3)
      cellVectorFields[2].set(win, [&](size_t id) {
        _model->cellGraph2Model.selectCellVectorField(id, 2);
        _model->cellGraph2Model.updateCellsVector();
      });
    new nanogui::Label(win, "Face Field");
    showFaceText = new nanogui::CheckBox(win, "values in text");
    showFaceFields = new nanogui::CheckBox(win, "draw field");
    faceFields.set(win, [&](size_t id) {
      _model->cellGraph2Model.selectFaceField(id);
      _model->cellGraph2Model.updateFaceRegions();
    });
    new nanogui::Label(win, "Vertex Field");
    showVertexText = new nanogui::CheckBox(win, "values in text");
    showVertexFields = new nanogui::CheckBox(win, "draw field");
    vertexFields.set(win, [&](size_t id) {
      _model->cellGraph2Model.selectVertexField(id);
      _model->cellGraph2Model.updateVertexRegions();
    });
  }
  ComboBox<size_t> cellFields, faceFields, vertexFields, particleFields;
  ComboBox<size_t> cellVectorFields[D];
  ComboBox<size_t> cellPointFields, facePointFields, vertexPointFields;
  nanogui::CheckBox *useCutPlanes, *cutPlanes[D],
      *showVertexText = nullptr, *showFaceText = nullptr,
      *showCellText = nullptr, *showParticles = nullptr,
      *showCellFiels = nullptr, *showFaceFields = nullptr,
      *showVertexFields = nullptr, *showParticleVelocities = nullptr,
      *showCellVectors = nullptr, *showGraph = nullptr, *showTopology = nullptr,
      *showSolids = nullptr;
  nanogui::detail::FormWidget<float> *textSize = nullptr,
                                     *particleVelocityScale = nullptr,
                                     *cellVectorScale = nullptr;

private:
  float _textSize = 0.00005;
  float _particleVelocityScale = 1.f;
  float _cellVectorScale = 1.f;
  SimulationModel<D> *_model = nullptr;
};

template <int D> struct IOMenu {
  explicit IOMenu(furoo::CustomSimulation<D> *sim,
                  std::function<void()> loadCFCallback,
                  std::function<void(std::string)> loadSimCallback)
      : _sim(sim), _loadCFCallback(std::move(loadCFCallback)),
        _loadSimCallback(std::move(loadSimCallback)) {}
  void set(nanogui::Widget *widget) {
    auto *popupButton =
        new nanogui::PopupButton(widget, "File", ENTYPO_ICON_DRIVE);
    auto *popup = popupButton->popup();
    popup->setLayout(new nanogui::GroupLayout());
    auto b = new nanogui::Button(popup, "Load Sim", ENTYPO_ICON_FOLDER);
    b->setCallback([&] {
      auto folder = nanogui::folder_dialog();
      if (!folder.empty()) {
        _sim->loadFromFile(folder[0] + "/configFile");
        _sim->setOutputPath(folder[0]);
        outputPathLabel->setValue(_sim->outputPath());
        _loadSimCallback(folder[0]);
      }
    });
    b = new nanogui::Button(popup, "Load CF", ENTYPO_ICON_PAPER_PLANE);
    b->setCallback([&] {
      _sim->loadFromFile(nanogui::file_dialog(
          {{"png", "Portable Network Graphics"}, {"*", "Text file"}}, false));
      _loadCFCallback();
      outputPathLabel->setValue(_sim->outputPath());
    });
    b = new nanogui::Button(popup, "Save CF", ENTYPO_ICON_SAVE);
    b->setCallback([&] {
      _sim->saveConfigFile(nanogui::file_dialog(
          {{"png", "Portable Network Graphics"}, {"txt", "Text file"}}, true));
    });
    b = new nanogui::Button(popup, "Output Files Path", ENTYPO_ICON_ARCHIVE);
    b->setCallback([&]() {
      auto f = nanogui::folder_dialog();
      if (!f.empty()) {
        outputPathLabel->setTooltip(f[0]);
        outputPathLabel->setValue(
            f[0].substr(std::max(static_cast<int>(f[0].size()) - 15, 0), 15));
        outputPath = f[0];
        _sim->setOutputPath(f[0]);
      }
    });
    outputPathLabel = new nanogui::TextBox(popup, "");
    _saveFramesImages = new nanogui::CheckBox(widget, "Store frames images");
    _saveFramesImages->setChecked(true);
    saveSubFrames = new nanogui::CheckBox(widget, "Store sub frames");
    saveParticles = new nanogui::CheckBox(widget, "Store particles");
    saveStructure = new nanogui::CheckBox(widget, "Store structure");
  }
  bool saveFrameImages() const { return _saveFramesImages->checked(); }
  std::string outputPath;

private:
  nanogui::CheckBox *saveParticles = nullptr, *saveStructure = nullptr,
                    *saveSubFrames = nullptr, *_saveFramesImages = nullptr;
  furoo::CustomSimulation<D> *_sim = nullptr;
  nanogui::TextBox *outputPathLabel = nullptr;
  std::function<void()> _loadCFCallback = nullptr;
  std::function<void(std::string)> _loadSimCallback = nullptr;
};

template <int D> struct SceneMenu {
  explicit SceneMenu(furoo::CustomSimulation<D> *simulation,
                     std::function<void()> loadScene)
      : sim(simulation), loadSceneCallback(std::move(loadScene)) {
    scenes.add("DROP", "DROP");
    scenes.add("LDB", "LDB");
    scenes.add("RDB", "RDB");
    scenes.add("DDB", "DDB");
  }
  void set(nanogui::FormHelper *gui, nanogui::Widget *widget) {
    gui->addGroup("Scene");
    particleSpacingLabel =
        gui->addVariable("Particle Spacing", particleSpacing);
    scenes.set(widget, [&](std::string value) {
      sim->loadScene(value, particleSpacing);
      loadSceneCallback();
    });
  }
  std::function<void()> loadSceneCallback = nullptr;

private:
  nanogui::detail::FormWidget<double> *particleSpacingLabel = nullptr;
  furoo::CustomSimulation<D> *sim;
  ComboBox<std::string> scenes;
  double particleSpacing = 0.03;
};

template <int D> struct SimulationPanel {
  SimulationPanel(furoo::CustomSimulation<D> *sim) : _sim(sim) {}
  void set(nanogui::FormHelper *gui) {
    gui->addGroup("Simulation Parameters");
    gui->addVariable("Time Step (dt)", dt);
    minLevelLabel = gui->addVariable("MinLevel", minLevel);
    minLevelLabel->setSpinnable(true);
    minLevelLabel->setMinMaxValues(1, 10);
    maxLevelLabel = gui->addVariable("MaxLevel", maxLevel);
    maxLevelLabel->setSpinnable(true);
    maxLevelLabel->setMinMaxValues(1, 10);
  }
  void updateLabels() {
    dt = _sim->timeStep();
    minLevel = _sim->minLevel();
    minLevelLabel->setValue(minLevel);
    maxLevel = _sim->maxLevel();
    maxLevelLabel->setValue(maxLevel);
  }

private:
  double dt = 0.001;
  int minLevel = 5, maxLevel = 6;
  nanogui::detail::FormWidget<int> *minLevelLabel = nullptr,
                                   *maxLevelLabel = nullptr;
  furoo::CustomSimulation<D> *_sim = nullptr;
};

template <int D> class SimulationStatistics {
public:
  SimulationStatistics(furoo::CustomSimulation<D> *sim) : _sim(sim) {}
  void set(nanogui::Widget *widget) {
    _particlesCount = new nanogui::Label(widget, "0 particles");
    _cellsCount = new nanogui::Label(widget, "0 cells");
    onScreen = new nanogui::CheckBox(widget, "show on screen");
  }
  void updateLabels() {
    _particlesCount->setCaption(
        furoo::concat(_sim->particles()->size(), " particles"));
    _cellsCount->setCaption(
        furoo::concat(_sim->cellGraph()->cellCount(), " cells"));
  }
  nanogui::CheckBox *onScreen = nullptr;

private:
  furoo::CustomSimulation<D> *_sim;
  nanogui::Label *_particlesCount, *_cellsCount, *_fluidCellsCount,
      *_minCellValueLabel, *_maxCellValueLabel, *_minVertexValueLabel,
      *_maxVertexValueLabel, *_minFaceValueLabel, *_maxFaceValueLabel;
};

template <int D> struct StepControl {
  explicit StepControl(furoo::CustomSimulation<D> *simulation,
                       std::function<void()> updateCallback)
      : _sim(simulation), _updateCallback(std::move(updateCallback)) {}
  void set(aergia::NanoGUIScreen *screen, int x, int y) {
    auto *_gui = new nanogui::FormHelper(screen);
    nanogui::ref<nanogui::Window> win =
        _gui->addWindow(Eigen::Vector2i(x, y), "Control");
    win->setLayout(new nanogui::GroupLayout());
    auto *wid = new nanogui::Widget(win);
    wid->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
                                          nanogui::Alignment::Middle, 0, 2));
    auto b = new nanogui::Button(wid, "", ENTYPO_ICON_CAMERA);
    b->setTooltip("Take screenshot.");
    b->setCallback([&]() {
      if (saveFrameImageCallback)
        saveFrameImageCallback(_sim->curFrame());
    });
    new nanogui::Label(wid, "Frame / Subframe");
    b = new nanogui::Button(wid, "", ENTYPO_ICON_CONTROLLER_JUMP_TO_START);
    b->setCallback([&]() {
      if (_sim->isRunning())
        return;
      _sim->reset();
      updateLabels();
      _updateCallback();
    });
    b = new nanogui::Button(wid, "", ENTYPO_ICON_CHEVRON_THIN_LEFT);
    b->setCallback([&]() {
      if (_sim->isRunning())
        return;
      updateLabels();
      _updateCallback();
    });
    curFrameLabel = new nanogui::TextBox(wid);
    curFrameLabel->setFixedWidth(80);
    curFrameLabel->setValue("");
    b = new nanogui::Button(wid, "", ENTYPO_ICON_CHEVRON_THIN_RIGHT);
    b->setCallback([&]() {
      if (_sim->isRunning())
        return;
      _sim->reset();
      updateLabels();
      _updateCallback();
    });
    b = new nanogui::Button(wid, "", ENTYPO_ICON_CONTROLLER_NEXT);
    b->setCallback([&]() {
      if (_sim->isRunning())
        return;
      _sim->advanceFrame();
      _updateCallback();
    });
    new nanogui::Label(wid, "Simulation");
    curStepLabel = new nanogui::TextBox(wid, std::string(""));
    curStepLabel->setFixedWidth(200);
    nextButton = new nanogui::Button(
        wid, (_sim->stepsCount() > 1) ? _sim->stepName(1) : std::string(""));
    nextButton->setFixedWidth(200);
    nextButton->setCallback([&]() {
      if (_sim->isRunning())
        return;
      _sim->advanceStep();
      _updateCallback();
      updateLabels();
      _sim->saveConfigFile("/configFile.txt");
    });
    playButton = new nanogui::Button(wid, "", ENTYPO_ICON_CONTROLLER_PLAY);
    playButton->setCallback([&]() {
      _sim->toggle();
      if (_sim->isRunning())
        playButton->setIcon(ENTYPO_ICON_CONTROLLER_PAUS);
      else
        playButton->setIcon(ENTYPO_ICON_CONTROLLER_PLAY);
      _sim->saveConfigFile("/configFile.txt");
    });
    slider.set(wid);
  }
  void updateLabels() {
    curFrameLabel->setValue(
        furoo::concat(_sim->curFrame(), '.', _sim->curSubFrame()));
    curStepLabel->setValue(_sim->stepName(
        _sim->curStep() > 0 ? _sim->curStep() - 1 : _sim->stepsCount() - 1));
    nextButton->setCaption(_sim->stepName(_sim->curStep()));
  }
  FrameSlider slider;

private:
  std::function<void(size_t)> saveFrameImageCallback;
  furoo::CustomSimulation<D> *_sim = nullptr;
  std::function<void()> _updateCallback;
  nanogui::Button *nextButton{}, *playButton{};
  nanogui::TextBox *curFrameLabel{};
  nanogui::TextBox *curStepLabel{};
};

template <int D>
void connectParticleFieldsToModelAndUI(
    furoo::ParticleSystem<double, D> *particles, ParticleSystemModel<D> &model,
    ComboBox<size_t> &fields) {
  for (size_t fieldId = 0; fieldId < particles->propertyCount(); fieldId++) {
  }
}

template <int D>
void connectStructureFieldsToModelAndUI(furoo::StructureInterface<D> *structure,
                                        CellGraphModel<D> &model,
                                        ComboBox<size_t> &cellFields,
                                        ComboBox<size_t> &faceFields,
                                        ComboBox<size_t> &vertexFields,
                                        ComboBox<size_t> &cellPointFields) {
  cellFields.clear();
  faceFields.clear();
  vertexFields.clear();
  cellPointFields.clear();
  auto graph = structure;
  for (size_t fieldId = 0; fieldId < graph->fieldCount(); fieldId++) {
    std::string fieldType = graph->fieldDataTypeName(fieldId);
    auto location = graph->fieldLocation(fieldId);
    size_t id = 0;
    bool isPointField = false;
    if (fieldType == "double") {
      auto field =
          new FieldModel<double, D>(graph->template field<double>(fieldId));
      field->structureFieldId = fieldId;
      id = model.addField(field, location);
    } else if (fieldType == "int") {
      auto field = new FieldModel<int, D>(graph->template field<int>(fieldId));
      field->structureFieldId = fieldId;
      id = model.addField(field, location);
    } else if (fieldType == "Definitions::Material") {
      auto field = new FieldModel<furoo::Definitions::Material, D>(
          graph->template field<furoo::Definitions::Material>(fieldId));
      field->structureFieldId = fieldId;
      field->isCategorical = true;
      field->colorMap[furoo::Definitions::Material::AIR] =
          aergia::COLOR_TRANSPARENT;
      field->colorMap[furoo::Definitions::Material::FLUID] = aergia::COLOR_BLUE;
      field->colorMap[furoo::Definitions::Material::SOLID] = aergia::COLOR_RED;
      id = model.addField(field, location);
    } else if (fieldType == "unsigned_char") {
      auto field = new FieldModel<unsigned char, D>(
          graph->template field<unsigned char>(fieldId));
      field->structureFieldId = fieldId;
      field->isCategorical = true;
      field->colorMap[0] = aergia::COLOR_TRANSPARENT;
      field->colorMap[1] = aergia::COLOR_GREEN;
      id = model.addField(field, location);
    } else if (fieldType == "Definitions::Boundary") {
      auto field = new FieldModel<furoo::Definitions::Boundary, D>(
          graph->template field<furoo::Definitions::Boundary>(fieldId));
      field->structureFieldId = fieldId;
      field->isCategorical = true;
      field->colorMap[furoo::Definitions::Boundary::NONE] =
          aergia::COLOR_TRANSPARENT;
      field->colorMap[furoo::Definitions::Boundary::NEUMANN] =
          aergia::COLOR_BLUE;
      field->colorMap[furoo::Definitions::Boundary::DIRICHLET] =
          aergia::COLOR_RED;
      id = model.addField(field, location);
    } else if (fieldType == "Point<double,D>") {
      auto field = new PointFieldModel<D>(
          graph->template field<furoo::Point<double, D>>(fieldId));
      id = model.addField(field, location, true);
      isPointField = true;
    }
    switch (location) {
    case furoo::Definitions::MeshLocation::CELL_CENTER:
      if (isPointField)
        cellPointFields.add(graph->fieldName(fieldId), id);
      else
        cellFields.add(graph->fieldName(fieldId), id);
      break;
    case furoo::Definitions::MeshLocation::FACE_CENTER:
    case furoo::Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    case furoo::Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    case furoo::Definitions::MeshLocation::DEPTH_FACE_CENTER:
      faceFields.add(graph->fieldName(fieldId), id);
      break;
    case furoo::Definitions::MeshLocation::VERTEX_CENTER:
      vertexFields.add(graph->fieldName(fieldId), id);
      break;
    default:
      exit(-1);
      break;
    }
  }
}

#endif // FUROO_NANOGUI_UTILS_H
