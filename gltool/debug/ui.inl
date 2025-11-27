template <int D> UI<D>::UI() = default;

template <int D>
void UI<D>::loadConfigFile(std::string filename, std::string frameDir,
                           int frameId) {
  _sim.loadFromFile(filename);
  if (frameId >= 0)
    _sim.loadFrame(frameDir, frameId);
  updateAll();
  // register solids
  auto solids = _sim.solids();
  for (size_t i = 0; i < solids.size(); i++) {
    _solids.emplace_back(std::make_shared<SolidModel<D>>(solids[i].get()));
  }
}

template <int D> void UI<D>::init() {
  _model.reset(new SimulationModel<D>(&_sim));
  _screen.reset(new aergia::NanoGUIScreen());
  _sceneMenu.reset(new SceneMenu<D>(&_sim, [&]() { updateAll(); }));
  _visPanel.reset(new VisPanel<D>(_model.get()));
  _stepControl.reset(new StepControl<D>(&_sim, [&]() { updateAll(); }));
  _stepControl->slider.callback = [&](size_t frame) {
    _sim.setCurFrameNumber(frame);
    _model->loadFrame(_inputDirectory, frame);
    _model->cellGraph2Model.clearFields();
    connectStructureFieldsToModelAndUI<D>(
        _sim.cellGraph(), _model->cellGraph2Model, _visPanel->cellFields,
        _visPanel->faceFields, _visPanel->vertexFields,
        _visPanel->cellPointFields);
    for (int d = 0; d < D; d++)
      _visPanel->cellVectorFields[d].copy(_visPanel->cellFields);
    updateAll();
    _screen->performLayout();
  };
  _simPanel.reset(new SimulationPanel<D>(&_sim));
  _simStats.reset(new SimulationStatistics<D>(&_sim));
  _ioMenu.reset(new IOMenu<D>(
      &_sim,
      [&]() {
        updateAll();
        // register solids
        auto solids = _sim.solids();
        for (size_t i = 0; i < solids.size(); i++) {
          _solids.emplace_back(
              std::make_shared<SolidModel<D>>(solids[i].get()));
        }
      },
      [&](std::string inputDirectory) {
        _inputDirectory = inputDirectory;
        // count particle files
        _stepControl->slider.maxFrame =
            countFrameFiles(inputDirectory, "p", ".pv");
        // count cell files
        size_t count = countFrameFiles(inputDirectory, "c");
        _stepControl->slider.maxFrame =
            std::min(_stepControl->slider.maxFrame, count);
        // count particle files
        count = countFrameFiles(inputDirectory, "f");
        _stepControl->slider.maxFrame =
            std::min(_stepControl->slider.maxFrame, count);
        _stepControl->slider.curFrame = 1;
        _model->loadFrame(inputDirectory, 1);
        _model->cellGraph2Model.clearFields();
        connectStructureFieldsToModelAndUI<D>(
            _sim.cellGraph(), _model->cellGraph2Model, _visPanel->cellFields,
            _visPanel->faceFields, _visPanel->vertexFields,
            _visPanel->cellPointFields);
        for (int d = 0; d < D; d++)
          _visPanel->cellVectorFields[d].copy(_visPanel->cellFields);
        updateAll();
        _screen->performLayout();
      }));
  connectStructureFieldsToModelAndUI<D>(
      _sim.cellGraph(), _model->cellGraph2Model, _visPanel->cellFields,
      _visPanel->faceFields, _visPanel->vertexFields,
      _visPanel->cellPointFields);
  for (int d = 0; d < D; d++)
    _visPanel->cellVectorFields[d].copy(_visPanel->cellFields);
  setupOptions();
  _screen->setVisible(true);
  _screen->performLayout();
}

template <int D> void UI<D>::setupOptions() {
  auto *_gui = new nanogui::FormHelper(_screen.get());
  nanogui::ref<nanogui::Window> win =
      _gui->addWindow(Eigen::Vector2i(10, 10), "Settings");
  win->setLayout(new nanogui::GroupLayout());
  // FILE POPUP MENU
  _ioMenu->set(win);
  _simPanel->set(_gui);
  _sceneMenu->set(_gui, win);
  _gui->addGroup("Tools & Visualization");
  {
    furoo::Vector<size_t, D> g;
    g[0] = _sim.fieldId("particleXVelocityProperty");
    g[1] = _sim.fieldId("particleYVelocityProperty");
    if (D == 3)
      g[2] = _sim.fieldId("particleZVelocityProperty");
    _model->particleSystemModel.addVectorGroup(g);
    _visPanel->particleFields.add("velocity", 0);
  }
  {
    furoo::Vector<size_t, D> g;
    g[0] = _sim.fieldId("particleXPressureGradientProperty");
    g[1] = _sim.fieldId("particleYPressureGradientProperty");
    if (D == 3)
      g[2] = _sim.fieldId("particleZPressureGradientProperty");
    _model->particleSystemModel.addVectorGroup(g);
    _visPanel->particleFields.add("PG", 1);
  }
  _visPanel->set(_gui, win);
  new nanogui::Button(win, "", ENTYPO_ICON_MOUSE_POINTER);
  //  _gui->addVariable("radius", this->_particleRadius);
  _testFunctionPanel.reset(new TestFunctionPanel<D>(
      _gui, win,
      [&](std::function<double(furoo::Point<double, D>)> f) {
        furoo::SimulationSteps::applyFunctionToCells(
            _sim.cellGraph(), _sim.fieldId("cellTestSolution"),
            _sim.fieldId("cellParticleCentroidFIeld"), f);
        updateAll();
        // apply to cells
      },
      [&](std::function<double(furoo::Point<double, D>)> f) {
        // apply to particles
        furoo::SimulationSteps::applyFunctionToParticles(
            _sim.particles(), _sim.fieldId("particleTestSolution"), f);
        updateAll();
      },
      [&](std::function<double(furoo::Point<double, D>)> f) {
        furoo::SimulationSteps::applyFunctionToCells(
            _sim.cellGraph(), _sim.fieldId("cellTestFunction"),
            _sim.fieldId("cellParticleCentroidFIeld"), f);
        updateAll();
      }));
  _testFunctionPanel->runGradXCallback = [&]() { _sim.run("testGradX"); };
  _testFunctionPanel->runGradYCallback = [&]() { _sim.run("testGradY"); };
  _testFunctionPanel->runLaplacianCallback = [&]() {
    _sim.run("testLaplacian");
  };
  _gui->addGroup("Statistics");
  _simStats->set(win);
  _gui->addGroup("Log");
  _logLines.resize(5, nullptr);
  for (size_t i = 0; i < 5; i++)
    _logLines[i] = new nanogui::Label(win, " ");
  _stepControl->set(_screen.get(), 201, 10);
}

template <int D> void UI<D>::log(std::string line) {
  for (size_t i = _logLines.size() - 1; i > 0; i--) {
    _logLines[i]->setCaption(_logLines[i - 1]->caption());
  }
  _logLines[0]->setCaption(line);
}

template <int D>
void UI<D>::draw(const aergia::CameraInterface *camera, ponos::Transform t) {
  if (_sim.isRunning()) {
    _sim.advanceStep();
    _model->cellGraph2Model.updateCellRegions();
    _model->cellGraph2Model.updateCellsVector();
    _model->cellGraph2Model.updateCellsEdges();
    _model->cellGraph2Model.updateFaceRegions();
    _model->cellGraphTopologyModel.updateEdges();
    _model->cellGraphTopologyModel.updateNodes();
    _model->particleSystemModel.update();
    _stepControl->updateLabels();
    _simStats->updateLabels();
  }
  if (_model->cellGraph2Model.cellVectorScale !=
      _visPanel->cellVectorScale->value()) {
    _model->cellGraph2Model.cellVectorScale =
        _visPanel->cellVectorScale->value();
    _model->cellGraph2Model.updateCellsVector();
  }
  if (_model->particleSystemModel.velocityScale !=
      _visPanel->particleVelocityScale->value()) {
    _model->particleSystemModel.velocityScale =
        _visPanel->particleVelocityScale->value();
    _model->particleSystemModel.updateVelocities();
  }
  if (_visPanel->showSolids->checked())
    for (size_t i = 0; i < _solids.size(); i++)
      _solids[i]->draw(camera, t);
  // solid->draw(camera, t);
  if (_model->cellGraph2Model.useCutPlanes !=
      _visPanel->useCutPlanes->checked()) {
    _model->cellGraph2Model.useCutPlanes = _visPanel->useCutPlanes->checked();
    _model->cellGraph2Model.updateCellsEdges();
    _model->cellGraph2Model.updateCellRegions();
  }
  for (int d = 0; d < D; d++)
    if (_model->cellGraph2Model.cutPlanes[d] !=
        _visPanel->cutPlanes[d]->checked()) {
      _model->cellGraph2Model.useCut[d] = _visPanel->cutPlanes[d]->checked();
      _model->cellGraph2Model.updateCellsEdges();
      _model->cellGraph2Model.updateCellRegions();
    }
  _model->cellGraph2Model.textSize = _visPanel->textSize->value();
  _model->cellGraph2Model.showText = _simStats->onScreen->checked();
  _model->cellGraph2Model.showCellText = _visPanel->showCellText->checked();
  _model->cellGraph2Model.showFaceText = _visPanel->showFaceText->checked();
  _model->cellGraph2Model.showVertexText = _visPanel->showVertexText->checked();
  _model->cellGraph2Model.showCellRegions = _visPanel->showCellFiels->checked();
  _model->cellGraph2Model.showCellVectors =
      _visPanel->showCellVectors->checked();
  _model->cellGraph2Model.showFaceRegions =
      _visPanel->showFaceFields->checked();
  _model->cellGraph2Model.showVertexRegions =
      _visPanel->showVertexFields->checked();
  if (_visPanel->showTopology->checked())
    _model->cellGraphTopologyModel.draw(camera, t);
  if (_visPanel->showGraph->checked())
    _model->cellGraph2Model.draw(camera, t);
  _model->particleSystemModel.showVelocities =
      _visPanel->showParticleVelocities->checked();
  if (_visPanel->showParticles->checked())
    _model->particleSystemModel.draw(camera, t);
  if (_sim.isRunning() && _ioMenu->saveFrameImages() && !_sim.curStep())
    saveFrameImage();
}

template <int D> bool UI<D>::intersect(const ponos::Ray3 &r, float *t) {
  return SceneObject::intersect(r, t);
}

template <int D> void UI<D>::saveFrameImage() {
  saveFrameImageCallback(
      furoo::concat(_sim.outputPath(), "/img/frame_", _sim.curFrame(), ".png"));
}

template <int D> void UI<D>::updateAll() {
  _simStats->updateLabels();
  _stepControl->updateLabels();
  _simPanel->updateLabels();
  _simStats->updateLabels();
  _model->updateParticles();
  _model->updateRegions();
  _model->updateEdges();
}

template <int D> void UI<D>::dropPath(std::string filename) {
  _sim.loadFromFile(filename);
  updateAll();
}

template <int D> void UI<D>::movePlane(int axis, int direction) {
  double step = 0.01;
  _model->cellGraph2Model.cutPlanes[axis] += direction * step;
  std::cerr << _model->cellGraph2Model.cutPlanes[axis] << std::endl;
  _model->cellGraph2Model.updateCellsEdges();
  _model->cellGraph2Model.updateCellRegions();
}