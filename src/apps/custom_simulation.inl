template <int D>
void CustomSimulation<D>::setCFLCallback(
    const std::function<double()> &cflCallback) {
  _cflCallback = cflCallback;
}

template <int D>
void CustomSimulation<D>::loadFromFile(const std::string &filename) {
  ArgsParser parser;
  parser.addIntArgument("minLevel", 5);
  parser.addIntArgument("maxLevel", 7);
  parser.addStringArgument("config", std::string());
  parser.addStringArgument("outputDir", std::string("output"));
  parser.addDoubleArgument("particleSpacing", 0.00390625);
  parser.addIntArgument("frames", 60);
  parser.addDoubleArgument("dt", 0.001);
  parser.addStringArgument("scene", "LDB");
  parser.addDoubleListArgument("injectors");
  parser.addDoubleListArgument("injectorTime");
  parser.addStringListArgument("algorithm");
  parser.addDoubleListArgument("fluidBoxes");
  parser.addIntListArgument("iFluidBoxes");
  parser.addIntListArgument("iFluidSpheres");
  parser.addDoubleListArgument("solidBoxes");
  parser.addIntListArgument("iSolidBoxes");
  parser.addBoolArgument("graded", true);
  parser.addStringArgument("rbfKernel", "GAUSS");
  parser.addStringArgument("base", "LINEAR");
  parser.addDoubleArgument("rbfRadius", 0.1);
  parser.addBoolArgument("adaptiveRBFRadius", false);
  parser.addStringArgument("fluidFile", std::string());
  parser.parse(filename);
  setTimeStep(parser.getDouble("dt"));
  _minLevel = parser.getInt("minLevel");
  _maxLevel = parser.getInt("maxLevel");
  _graded = parser.getBool("graded");
  _adaptiveRBFRadius = parser.getBool("adaptiveRBFRadius");
  _rbfRadius = parser.getDouble("rbfRadius");
  _frames = parser.getInt("frames");
  _outputPath = parser.getString("outputDir");
  if (parser.getString("rbfKernel") == "GAUSSIAN")
    _rbfKernel.reset(new GaussianKernel<Point<double, D>>(_rbfRadius));
  else
    _rbfKernel.reset(new CubicKernel<Point<double, D>>());
  _rbf.reset(new DifferentialRBF<D>(*_rbfKernel.get()));
  if (parser.getString("base") == "QUADRATIC")
    _rbf->setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
  loadInjectors(parser.getDoubleList("injectors"));
  loadInjectorTime(parser.getDoubleList("injectorTime"));
  std::string scene = parser.getString("scene");
  double particleSpacing = parser.getDouble("particleSpacing");
  // add domain solid walls
  loadSceneSolids(parser.getDoubleList("solidBoxes"));
  loadSceneSolids(parser.getIntList("iSolidBoxes"));
  if (D == 2) {
    int maxCoord = 1 << _maxLevel;
    int maxLevel = _maxLevel;
    std::vector<int> solidWallsDescription = {
        // resolution
        maxLevel, maxLevel,
        // bottom
        0, 0, maxCoord, 1,
        // left
        0, 0, 1, maxCoord,
        // right
        maxCoord - 1, 0, maxCoord, maxCoord,
        // top
        0, maxCoord - 1, maxCoord, maxCoord};
    loadSceneSolids(solidWallsDescription);
  } else {
    int maxCoord = 1 << _maxLevel;
    int maxLevel = _maxLevel;
    std::vector<int> solidWallsDescription = {
        // resolution
        maxLevel, maxLevel, maxLevel,
        // bottom
        0, 0, 0, maxCoord, 1, maxCoord,
        // left
        0, 0, 0, 1, maxCoord, maxCoord,
        // right
        maxCoord - 1, 0, 0, maxCoord, maxCoord, maxCoord,
        // top
        0, maxCoord - 1, 0, maxCoord, maxCoord, maxCoord,
        // back
        0, 0, 0, maxCoord, maxCoord, 1,
        // front
        0, 0, maxCoord - 1, maxCoord, maxCoord, maxCoord};
    // loadSceneSolids(solidWallsDescription);
    // // TODO: change to load implicit scene later
    // loadScene(parser.getDoubleList("fluidBoxes"), particleSpacing);
    // loadScene(scene, particleSpacing);
  }
  // loadImplicitScene(parser.getIntList("iFluidBoxes"), particleSpacing);
  loadScene(parser.getIntList("iFluidBoxes"),
            parser.getIntList("iFluidSpheres"), particleSpacing);
  if (parser.getString("fluidFile").size()) {
    ParticleInjector<D> injector(_particles.get(), particleSpacing);
    injector.fromFile(parser.getString("fluidFile"));
    SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                               _particles.get());
    SimulationSteps::refineFluid(_particles.get(), _graph.get(),
                                 _fieldIds["cellMaterialField"], _maxLevel);
    SimulationSteps::markFluidSurface(_graph.get(),
                                      _fieldIds["cellSurfaceMaskField"],
                                      _fieldIds["cellMaterialField"]);
    SimulationSteps::markSolids<D>(_graph.get(), _fieldIds["cellMaterialField"],
                                   _solids);
    SimulationSteps::updateTree(
        _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
        _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
        _minLevel, _maxLevel);
    SimulationSteps::reseedParticles(
        _graph.get(), _fieldIds["cellMaterialField"],
        _fieldIds["cellSurfaceMaskField"], _particles.get());
  }
  //  addStep("begin", [](double dt) {});
  for (const auto &s : parser.getStringList("algorithm"))
    addStep(s);
  // save config file
  std::ifstream infp(filename, std::ifstream::in);
  std::stringstream fileContent;
  while (infp.good()) {
    std::string line;
    std::getline(infp, line);
    fileContent << line << "\n";
  }
  infp.close();
  std::ofstream outfp(_outputPath + "/configFile", std::ofstream::out);
  outfp << fileContent.str();
  outfp.close();
  // auto &surfaceField =
  //     *_graph->template field<unsigned
  //     char>(_fieldIds["cellSurfaceMaskField"]);

  // auto &materialField = *_graph->template field<Definitions::Material>(
  //     _fieldIds["cellMaterialField"]);
  // _graph->iterateCells([&](size_t cellId) {
  //   if (materialField[cellId] == Definitions::Material::FLUID)
  //     surfaceField[cellId] = 1;
  // });
}

template <int D>
void CustomSimulation<D>::saveConfigFile(const std::string &filename) {
  NOT_IMPLEMENTED();
}

template <>
inline void CustomSimulation<2>::loadScene(const std::string &scene,
                                           double particleSpacing) {
  if (scene == "NONE")
    return;
  _particles->clear();
  int total = 0;
  Injector2 injector(_particles.get(), particleSpacing);
  if (scene == "DROP") {
    // TANK
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.001, 0.001), Point2d(0.999, 0.2)));
    // DROP
    total +=
        injector.setupBoxShape(BBox2d(Point2d(0.4, 0.4), Point2d(0.6, 0.6)));
  } else if (scene == "LDB")
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.001, 0.001), Point2d(0.2, 0.5)));
  else if (scene == "RDB")
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.8, 0.001), Point2d(0.999, 0.5)));
  else if (scene == "DDB") {
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.001, 0.001), Point2d(0.2, 0.5)));
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.8, 0.001), Point2d(0.999, 0.5)));
  } else if (scene == "TANK") {
    total += injector.setupBoxShape(
        BBox2d(Point2d(0.001, 0.001), Point2d(0.999, 0.3)));
  }
  SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                             _particles.get());
  SimulationSteps::refineFluid(_particles.get(), _graph.get(),
                               _fieldIds["cellMaterialField"], _maxLevel);
  SimulationSteps::markFluidSurface(_graph.get(),
                                    _fieldIds["cellSurfaceMaskField"],
                                    _fieldIds["cellMaterialField"]);
  SimulationSteps::updateTree(
      _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
      _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
      _minLevel, _maxLevel);
  SimulationSteps::reseedParticles(_graph.get(), _fieldIds["cellMaterialField"],
                                   _fieldIds["cellSurfaceMaskField"],
                                   _particles.get());
}

template <>
inline void CustomSimulation<3>::loadScene(const std::string &scene,
                                           double particleSpacing) {
  if (scene == "NONE")
    return;
  _particles->clear();
  int total = 0;
  Injector3 injector(_particles.get(), particleSpacing);
  if (scene == "DROP") {
    // TANK
    //    total +=
    //        injector.setupBoxShape(BBox3d(Point3d(0.001, 0.001, 0.001),
    //                                      Point3d(0.999, 0.2, 0.999)));
    // DROP
    total += injector.setupBoxShape(
        BBox3d(Point3d(0.4, 0.4, 0.4), Point3d(0.6, 0.6, 0.6)));
  } else if (scene == "LDB")
    total += injector.setupBoxShape(
        BBox3d(Point3d(0.001, 0.001, 0.001), Point3d(0.2, 0.5, 0.2)));
  else if (scene == "RDB")
    total += injector.setupBoxShape(
        BBox3d(Point3d(0.8, 0.001, 0.8), Point3d(0.999, 0.5, 0.999)));
  else if (scene == "DDB") {
    total += injector.setupBoxShape(
        BBox3d(Point3d(0.001, 0.001, 0.001), Point3d(0.2, 0.5, 0.2)));
    total += injector.setupBoxShape(
        BBox3d(Point3d(0.8, 0.001, 0.8), Point3d(0.999, 0.5, 0.999)));
  }
  std::cerr << total << " particles injected\n";
  SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                             _particles.get());
  SimulationSteps::refineFluid(_particles.get(), _graph.get(),
                               _fieldIds["cellMaterialField"], _maxLevel);
  SimulationSteps::markFluidSurface(_graph.get(),
                                    _fieldIds["cellSurfaceMaskField"],
                                    _fieldIds["cellMaterialField"]);
  SimulationSteps::updateTree(
      _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
      _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
      _minLevel, _maxLevel);

  SimulationSteps::reseedParticles(_graph.get(), _fieldIds["cellMaterialField"],
                                   _fieldIds["cellSurfaceMaskField"],
                                   _particles.get());
}

template <int D>
void CustomSimulation<D>::loadScene(std::vector<double> boxes,
                                    double particleSpacing) {
  if (boxes.empty())
    return;
  _particles->clear();
  int total = 0;
  ParticleInjector<D> injector(_particles.get(), particleSpacing);
  THROW(boxes.size() % (D * 2) == 0,
        "CustomSimulation<>::loadScene Invalid box description");
  size_t n = boxes.size() / (D * 2);
  for (size_t i = 0; i < n; i++)
    total += injector.setupBoxShape(
        BBox<double, D>(Point<double, D>(&boxes[i * D * 2]),
                        Point<double, D>(&boxes[i * D * 2 + D])));
  std::cerr << total << " particles injected\n";
  SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                             _particles.get());
  SimulationSteps::refineFluid(_particles.get(), _graph.get(),
                               _fieldIds["cellMaterialField"], _maxLevel);
  SimulationSteps::markFluidSurface(_graph.get(),
                                    _fieldIds["cellSurfaceMaskField"],
                                    _fieldIds["cellMaterialField"]);
  SimulationSteps::updateTree(
      _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
      _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
      _minLevel, _maxLevel);

  SimulationSteps::reseedParticles(_graph.get(), _fieldIds["cellMaterialField"],
                                   _fieldIds["cellSurfaceMaskField"],
                                   _particles.get());
}

template <int D>
void CustomSimulation<D>::loadScene(std::vector<int> boxes,
                                    std::vector<int> spheres,
                                    double particleSpacing) {
  if (boxes.size() < D)
    return;
  _particles->clear();
  int total = 0;
  ParticleInjector<D> injector(_particles.get(), particleSpacing);
  THROW((boxes.size() - D) % (D * 2 + 1) == 0,
        "CustomSimulation<>::loadScene Invalid box description");
  size_t n = (boxes.size() - D) / (D * 2 + 1);
  Vector<double, D> step;
  for (int d = 0; d < D; d++)
    step[d] = 1. / (1 << boxes[d]);
  for (size_t i = 0; i < n; i++) {
    Point<double, D> pmin, pmax;
    for (int d = 0; d < D; d++) {
      pmin[d] = boxes[D + i * (D * 2 + 1) + d] * step[d];
      pmax[d] = boxes[D + i * (D * 2 + 1) + D + d] * step[d];
    }
    double group = boxes[D + i * (D * 2 + 1) + D + D];
    total += injector.setupBoxShape(BBox<double, D>(pmin, pmax), group);
    std::cerr << "setup fluid box " << BBox<double, D>(pmin, pmax) << std::endl;
  }
  if (spheres.size()) {
    THROW((spheres.size() - D) % (D + 1) == 0,
          "CustomSimulation<>::loadScene Invalid sphere description");
    n = (spheres.size() - D) / (D + 1);
    for (int d = 0; d < D; d++)
      step[d] = 1. / (1 << spheres[d]);
    for (size_t i = 0; i < n; i++) {
      Point<double, D> center;
      for (int d = 0; d < D; d++)
        center[d] = spheres[D + i * (D + 1) + d] * step[d];
      double radius = spheres[D + i * (D + 1) + D] * step[0];
      total += injector.setupCircleShape(center, radius);
      std::cerr << "setup fluid sphere with radius " << radius << " center "
                << center << std::endl;
    }
  }
  std::cerr << total << " particles injected\n";
  SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                             _particles.get());
  SimulationSteps::refineFluid(_particles.get(), _graph.get(),
                               _fieldIds["cellMaterialField"], _maxLevel);
  SimulationSteps::markFluidSurface(_graph.get(),
                                    _fieldIds["cellSurfaceMaskField"],
                                    _fieldIds["cellMaterialField"]);
  SimulationSteps::updateTree(
      _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
      _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
      _minLevel, _maxLevel);
  SimulationSteps::reseedParticles(_graph.get(), _fieldIds["cellMaterialField"],
                                   _fieldIds["cellSurfaceMaskField"],
                                   _particles.get());
}

template <int D>
void CustomSimulation<D>::loadImplicitScene(std::vector<int> boxes,
                                            double particleSpacing) {
  if (boxes.size() < D)
    return;
  _particles->clear();
  try {
    SimulationSteps::initImplicitScene(
        boxes, _solids, _maxLevel, _minLevel,
        _fieldIds["cellSurfaceTDistanceField"], _graph.get(), _particles.get(),
        _fieldIds["cellMaterialField"], _fieldIds["cellSurfaceMaskField"],
        _graded);
  } catch (std::string &e) {
    std::cerr << e << std::endl;
    THROW(false, "CustomSimulation::loadImplicScene not able to init scene");
  }
  try {
    SimulationSteps::updateTreeCoarseBottom(
        _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
        _fieldIds["cellSurfaceTDistanceField"], _fieldIds["cellMaterialField"],
        _minLevel, _maxLevel);
  } catch (std::string &e) {
    std::cerr << e << std::endl;
    THROW(false, "CustomSimulation;:loadImplicitScene not able to update tree")
  }
}

template <int D>
void CustomSimulation<D>::loadSceneSolids(std::vector<double> boxes) {
  if (boxes.empty())
    return;
  THROW(boxes.size() % (D * 2) == 0,
        "CustomSimulation<>::loadSceneSolids Invalid box description");
  size_t n = boxes.size() / (D * 2);
  for (size_t i = 0; i < n; i++)
    _solids.push_back(std::make_shared<SolidBox<D>>(
        BBox<double, D>(Point<double, D>(&boxes[i * D * 2]),
                        Point<double, D>(&boxes[i * D * 2 + D]))));
}

template <int D>
void CustomSimulation<D>::loadSceneSolids(std::vector<int> boxes) {
  if (boxes.size() < D)
    return;
  THROW((boxes.size() - D) % (D * 2) == 0,
        "CustomSimulation<>::loadSceneSolids Invalid box description");
  size_t n = (boxes.size() - D) / (D * 2);
  Vector<double, D> step;
  for (int d = 0; d < D; d++)
    step[d] = 1. / (1 << boxes[d]);
  for (size_t i = 0; i < n; i++) {
    Point<double, D> pmin, pmax;
    for (int d = 0; d < D; d++) {
      pmin[d] = boxes[D + i * D * 2 + d] * step[d];
      pmax[d] = boxes[D + i * D * 2 + D + d] * step[d];
    }
    _solids.push_back(
        std::make_shared<SolidBox<D>>(BBox<double, D>(pmin, pmax)));
    std::cerr << "// setup solid box " << BBox<double, D>(pmin, pmax)
              << std::endl;
    if (D == 3)
      std::cerr << " auto b" << i << " = Box3::builder()\n"
                << ".withLowerCorner({" << pmin[0] << ", " << pmin[1] << ", "
                << pmin[2] << "})"
                << ".withUpperCorner({" << pmax[0] << ", " << pmax[1] << ", "
                << pmax[2] << "})"
                << ".makeShared();\n";
  }
  if (D == 3) {
    std::cerr << "auto wallSet = "
                 "ImplicitSurfaceSet3::builder().withExplicitSurfaces({";
    for (size_t i = 0; i < n; i++)
      std::cerr << "b" << i << ", ";
    std::cerr << "}).makeShared();";
  }
}

template <int D>
void CustomSimulation<D>::loadFrame(std::string inputDirectory, size_t frame) {
  std::cerr << "loading frame\n";
  auto f = concat(std::setfill('0'), std::setw(6), frame);
  furoo::IO::loadFromPV(concat(inputDirectory, "/pFrame", f), _particles.get());
  furoo::IO::loadCellGraph(concat(inputDirectory, "/cFrame", f), _graph.get());
  _particles->update();
  // furoo::IO::loadFields(concat(inputDirectory, "/f", f), _graph.get());
  _curFrame = frame;
  _curStep = 0;
  _curSubFrame = 0;
  _cflSteps = 0;
  SimulationSteps::markCells(_graph.get(), _fieldIds["cellMaterialField"],
                             _particles.get());
  auto &surface =
      *_graph->template field<unsigned char>(_fieldIds["cellSurfaceMaskField"]);
  _graph->iterateCells([&](size_t cellId) { surface[cellId] = 0; });
  SimulationSteps::markFluidSurface(_graph.get(),
                                    _fieldIds["cellSurfaceMaskField"],
                                    _fieldIds["cellMaterialField"]);
  SimulationSteps::markSolids<D>(_graph.get(), _fieldIds["cellMaterialField"],
                                 _solids);
  SimulationSteps::markFaces<D>(_graph.get(), _fieldIds["cellMaterialField"],
                                _fieldIds["faceMaterialField"],
                                _fieldIds["faceBoundaryTypeField"],
                                _domainBoundaries.data());
}

template <int D>
void CustomSimulation<D>::loadInjectors(std::vector<double> injectors) {
  if (injectors.empty())
    return;
  size_t dataSize = D * 3 + 1;
  THROW(injectors.size() % dataSize == 0,
        "CustomSimulation<>::loadInjectors Invalid injector description");
  size_t n = injectors.size() / dataSize;
  for (size_t i = 0; i < n; i++) {
    Vector<double, D> size, velocity, pos;
    double inflowRate = injectors[(i + 1) * dataSize - 1];
    for (int d = 0; d < D; d++) {
      size[d] = injectors[i * dataSize + d];
      pos[d] = injectors[i * dataSize + D + d];
      velocity[d] = injectors[i * dataSize + 2 * D + d];
    }
    auto m = Transform<D>::translate(pos).matrix() *
             Transform<D>::scale(size).matrix() *
             Transform<D>::translate(Vector<double, D>(-0.5)).matrix();
    Transform<D> transform(m, furoo::inverse(m));
    _injectors.emplace_back(BBox<double, D>::squareBox(1.), transform, velocity,
                            inflowRate, _particles.get());
    std::cerr << "setup injector " << transform(BBox<double, D>::squareBox(1.))
              << " vel " << velocity << " flowRate " << inflowRate << std::endl;
  }
}

template <int D>
void CustomSimulation<D>::loadInjectorTime(std::vector<double> injectorTimes) {
  if (injectorTimes.empty())
    return;
  THROW(injectorTimes.size() % 3 == 0,
        "CustomSimulation<>::loadInjectorTime invalid description");
  for (size_t j = 0; j < injectorTimes.size() / 3; j++) {
    int i;
    double starting, ending;
    i = static_cast<int>(injectorTimes[j * 3 + 0]);
    starting = injectorTimes[j * 3 + 1];
    ending = injectorTimes[j * 3 + 2];
    if (ending < 0) {
      ending = INFINITY;
    }
    THROW(i <= _injectors.size(),
          "CustomSimulation::loadInjectorTime invalid injector index");
    THROW(starting >= 0,
          "CustomSimulation::loadInjectorTime starting time lower than 0");
    _injectors[i].setStartEndTime(starting, ending);
    std::cerr << "Injector " << i << " interval: " << starting << ", " << ending
              << std::endl;
  }
}

template <int D>
void CustomSimulation<D>::addStep(std::string name,
                                  std::function<void(double)> callBack) {
  _steps.emplace_back(name, callBack);
}

template <int D> void CustomSimulation<D>::addStep(std::string name) {
  if (_stepMap.find(name) == _stepMap.end()) {
    std::cerr << "step name: " << name << std::endl;
    THROW(false, "CustomSimulation<D>::addStep step name not registered!\n");
  }
  _steps.emplace_back(name, _stepMap[name]);
}

template <int D> size_t CustomSimulation<D>::stepsCount() const {
  return _steps.size();
}

template <int D> std::string CustomSimulation<D>::stepName(size_t i) const {
  return _steps[i % _steps.size()].first;
}

template <int D> size_t CustomSimulation<D>::curStep() const {
  return _curStep;
}

template <int D> size_t CustomSimulation<D>::curFrame() const {
  return _curFrame;
}

template <int D> size_t CustomSimulation<D>::curSubFrame() const {
  return _curSubFrame;
}

template <int D> void CustomSimulation<D>::toggle() {
  _isSimulating = !_isSimulating;
}

template <int D> bool CustomSimulation<D>::isRunning() const {
  return _isSimulating;
}

template <int D> void CustomSimulation<D>::setTimeStep(double dt) { _dt = dt; }

template <int D> void CustomSimulation<D>::run() {
  omp_set_num_threads(10);
  Eigen::initParallel();
  while (curFrame() < _frames)
    advanceFrame();
}

template <int D>
void CustomSimulation<D>::stepNames(std::vector<std::string> &names) const {
  for (auto &name : _stepMap)
    names.emplace_back(name.first);
}

template <int D>
const std::vector<std::string> &CustomSimulation<D>::cellFieldNames() const {
  return _cFieldsNames;
}

template <int D>
const std::vector<std::string> &CustomSimulation<D>::faceFieldNames() const {
  return _fFieldsNames;
}

template <int D>
const std::vector<std::string> &CustomSimulation<D>::vertexFieldNames() const {
  return _vFieldsNames;
}

template <int D>
const std::vector<std::string> &
CustomSimulation<D>::particleFieldNames() const {
  return _pFieldsNames;
}

template <int D>
const std::vector<std::string> &CustomSimulation<D>::fieldNames() const {
  return _fieldNames;
}

template <int D> ParticleSystem<double, D> *CustomSimulation<D>::particles() {
  return _particles.get();
}

template <int D>
std::vector<std::shared_ptr<SolidInterface<D>>> &CustomSimulation<D>::solids() {
  return _solids;
}

template <int D> StructureInterface<D> *CustomSimulation<D>::cellGraph() {
  return _graph.get();
}

template <int D>
int CustomSimulation<D>::fieldId(const std::string &fieldName) const {
  const auto it = _fieldIds.find(fieldName);
  if (it == _fieldIds.end())
    return -1;
  return static_cast<int>(it->second);
}

template <int D> void CustomSimulation<D>::run(const std::string &stepName) {
  LOG_INFO("run step " << stepName);
  _stepMap[stepName](_dt);
}

template <int D> void CustomSimulation<D>::reset() {
  _curStep = _curSubFrame = _curFrame = 0;
  _cflSteps = static_cast<size_t>(
      std::ceil(std::max(_dt / (_cflCallback() * 0.5), 1.)));
  _cflDT = _dt / _cflSteps;
}

template <int D> void CustomSimulation<D>::setCurFrameNumber(size_t frame) {
  reset();
  _curFrame = frame;
}

template <int D> size_t CustomSimulation<D>::advanceFrame() {
  size_t frame = _curFrame;
  while (frame == _curFrame)
    advanceStep();
  return frame;
}

template <int D> size_t CustomSimulation<D>::advanceStep() {
  LOG_INFO("advace step: curStep " << _curStep << " curSubFrame "
                                   << _curSubFrame << " cflSteps "
                                   << _cflSteps);
  static double curDt = _dt;
  if (!_curStep) {
    {
      CodeProfiler::profile("Cell count", _graph->cellCount());
      CodeProfiler::profile("Particle count", _particles->size());
      CodeProfiler::createSnapshot(concat("", _curFrame, ".", _curSubFrame));
      CodeProfiler::reset();
    }
    _curSubFrame++;
    if (curDt <= 0.) {
      curDt = _dt;
      _curSubFrame = 0;
      _curFrame++;
      { // save frame
        std::cerr << "save final frame " << curFrame() << " " << curSubFrame()
                  << std::endl;
        std::cerr << "====================" << std::endl;
        std::ostringstream particlesFile, cellGraphFile, fieldsFile;
        particlesFile << _outputPath << "/pFrame" << std::setfill('0')
                      << std::setw(6) << curFrame();
        cellGraphFile << _outputPath << "/cFrame" << std::setfill('0')
                      << std::setw(6) << curFrame();
        IO::saveCellGraph<D>(cellGraphFile.str(), _graph.get());
        IO::saveToPV(particlesFile.str().c_str(), _particles.get());
      }
      // #ifdef FUROO_PROFILE
      CodeProfiler::saveHistory(concat(_outputPath, "/profiler.csv"), true);
      // #endif
    }
    LOG_INFO("computing cflSteps, remaining time:" << curDt);
    _cflSteps = static_cast<size_t>(
        std::ceil(std::max(curDt / (_cflCallback() * .5), 1.)));
    _cflDT = std::min(curDt / _cflSteps, _dt);
    LOG_INFO("frame " << _curFrame << " with " << _cflSteps << " sub-frames");
    curDt -= _cflDT;
  }
  LOG_INFO("\tstep " << _steps[_curStep].first << " with dt " << _cflDT);
#ifdef FUROO_DEBUG
  IO::saveFields(concat(_outputPath, "/lf"), _graph.get());
  IO::saveCellGraph<D>(concat(_outputPath, "/lc"), _graph.get());
  IO::saveToPV(concat(_outputPath, "/lp").c_str(), _particles.get());
#endif
  // #ifdef FUROO_PROFILE
  Timer timer;
  // #endif
  timer.reset();
  std::clock_t c_start = std::clock();
  //   double start_time = omp_get_wtime();
  {
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::cout << ctime(&timenow) << std::endl;
  }
  _steps[_curStep].second(_cflDT);
  {
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::cout << ctime(&timenow) << std::endl;
  }
  std::clock_t c_end = std::clock();

  //  double time = omp_get_wtime() - start_time;
  long double time_elapsed_ms =
      1000.0 * (c_end - c_start) / (double)CLOCKS_PER_SEC;
  std::cerr << "CPU time used: " << time_elapsed_ms << " ms\n";
  static long double total = 0;
  total += time_elapsed_ms;
  std::cerr << "total time" << total << std::endl;

  // #ifdef FUROO_PROFILE
  CodeProfiler::profile(concat(_steps[_curStep].first, ""),
                        timer.ellapsedTimeInSeconds());
  // LOG_INFO("\tin " << timer.ellapsedTimeInSeconds() << " seconds");
  // LOG_INFO("\tin " << time << " seconds");
  // #endif
  _curStep = (_curStep + 1) % _steps.size();
  return _curStep;
}

template <int D> void CustomSimulation<D>::setupFields() {
  _fieldIds["particleXVelocityProperty"] = _particles->addScalarProperty(0.);
  _fieldIds["particleYVelocityProperty"] = _particles->addScalarProperty(0.);
  if (D == 3)
    _fieldIds["particleZVelocityProperty"] = _particles->addScalarProperty(0.);
  _fieldIds["particleGroupProperty"] = _particles->addScalarProperty(0.);
  //  _fieldIds["particleDivergenceField"] = _particles->addScalarProperty(0.);
  //  _fieldIds["particleXPressureGradientProperty"] =
  //      _particles->addScalarProperty(0.);
  //  _fieldIds["particleYPressureGradientProperty"] =
  //      _particles->addScalarProperty(0.);
  //  if (D == 3)
  //    _fieldIds["particleZPressureGradientProperty"] =
  //        _particles->addScalarProperty(0.);
  _fieldIds["cellXVelocityField"] =
      _graph->template addCellField<double>(0., "xVel");
  _fieldIds["cellYVelocityField"] =
      _graph->template addCellField<double>(0., "yVel");
  if (D == 3)
    _fieldIds["cellZVelocityField"] =
        _graph->template addCellField<double>(0., "zVel");
  _fieldIds["cellPressureField"] =
      _graph->template addCellField<double>(0., "pressure");
  _fieldIds["cellPressureGradientXField"] =
      _graph->template addCellField<double>(0., "pGradientX");
  _fieldIds["cellPressureGradientYField"] =
      _graph->template addCellField<double>(0., "pGradientY");
  if (D == 3)
    _fieldIds["cellPressureGradientZField"] =
        _graph->template addCellField<double>(0., "pGradientZ");
  _fieldIds["cellDivergenceField"] =
      _graph->template addCellField<double>(0., "div");
  _fieldIds["cellMaterialField"] =
      _graph->template addCellField<Definitions::Material>(
          Definitions::Material::AIR, "material");
  _fieldIds["cellSurfaceMaskField"] =
      _graph->template addCellField<unsigned char>(0, "surface");
  _fieldIds["cellSurfaceTDistanceField"] =
      _graph->template addCellField<int>(INT_INFINITY, "surfDist");
  //  _fieldIds["cellParticleCentroidField"] =
  //      _graph->template addCellField<Point<double, D>>(
  //          Point<double, D>(INFINITY), "cellParticleCentroid");
  //  _fieldIds["cellTestFunction"] =
  //      _graph->template addCellField<double>(0., "test");
  //  _fieldIds["cellTestResult"] =
  //      _graph->template addCellField<double>(0., "test_numeric");
  //  _fieldIds["cellTestSolution"] =
  //      _graph->template addCellField<double>(0., "test_solution");
  _fieldIds["facePressureField"] =
      _graph->template addFaceField<double>(0., "pressure");
  _fieldIds["facePressureGradientXField"] =
      _graph->template addFaceField<double>(
          0., "pGradientX", Definitions::MeshLocation::VERTICAL_FACE_CENTER);
  // 0., "pGradientX", Definitions::MeshLocation::FACE_CENTER);
  _fieldIds["facePressureGradientYField"] =
      _graph->template addFaceField<double>(
          0., "pGradientY", Definitions::MeshLocation::HORIZONTAL_FACE_CENTER);
  // 0., "pGradientY", Definitions::MeshLocation::FACE_CENTER);
  if (D == 3)
    _fieldIds["facePressureGradientZField"] =
        _graph->template addFaceField<double>(
            0., "pGradientZ", Definitions::MeshLocation::DEPTH_FACE_CENTER);
  // 0., "pGradientZ", Definitions::MeshLocation::FACE_CENTER);
  _fieldIds["faceXVelocityField"] = _graph->template addFaceField<double>(
      0., "velX", Definitions::MeshLocation::VERTICAL_FACE_CENTER);
  // 0., "velX", Definitions::MeshLocation::FACE_CENTER);
  _fieldIds["faceYVelocityField"] = _graph->template addFaceField<double>(
      0., "velY", Definitions::MeshLocation::HORIZONTAL_FACE_CENTER);
  // 0., "velY", Definitions::MeshLocation::FACE_CENTER);
  if (D == 3)
    _fieldIds["faceZVelocityField"] = _graph->template addFaceField<double>(
        0., "velZ", Definitions::MeshLocation::DEPTH_FACE_CENTER);
  // 0., "velZ", Definitions::MeshLocation::FACE_CENTER);
  _fieldIds["faceMaterialField"] =
      _graph->template addFaceField<Definitions::Material>(
          Definitions::Material::AIR, "material");
  _fieldIds["faceBoundaryTypeField"] =
      _graph->template addFaceField<Definitions::Boundary>(
          Definitions::Boundary::NONE, "boundary");
  //  _fieldIds["faceFlowField"] =
  //      _graph->template addFaceField<double>(0., "flow");
  //  _fieldIds["faceMidPointField"] =
  //      _graph->template addFaceField<Point<double, D>>(Point<double, D>(0.0),
  //                                                      "faceMidPoint");
  //  _fieldIds["vertexXVelocityField"] =
  //      _graph->template addVertexField<double>(0., "velX");
  //  _fieldIds["vertexYVelocityField"] =
  //      _graph->template addVertexField<double>(0., "velY");
  //  if (D == 3)
  //    _fieldIds["vertexZVelocityField"] =
  //        _graph->template addVertexField<double>(0., "velZ");
}

template <int D> CustomSimulation<D>::CustomSimulation() {
  // TODO break this method into smaller pieces of code
  _pFieldsNames.emplace_back("particleXVelocityProperty");
  _pFieldsNames.emplace_back("particleYVelocityProperty");
  if (D == 3)
    _pFieldsNames.emplace_back("particleZVelocityProperty");
  _pFieldsNames.emplace_back("particleXPressureGradientProperty");
  _pFieldsNames.emplace_back("particleYPressureGradientProperty");
  if (D == 3)
    _pFieldsNames.emplace_back("particleZPressureGradientProperty");
  _pFieldsNames.emplace_back("particleTestFunction");
  _pFieldsNames.emplace_back("particleDivergenceField");
  _cFieldsNames.emplace_back("cellXVelocityField");
  _cFieldsNames.emplace_back("cellYVelocityField");
  if (D == 3)
    _cFieldsNames.emplace_back("cellZVelocityField");
  _cFieldsNames.emplace_back("cellPressureField");
  _cFieldsNames.emplace_back("cellPressureGradientXField");
  _cFieldsNames.emplace_back("cellPressureGradientYField");
  if (D == 3)
    _cFieldsNames.emplace_back("cellPressureGradientZField");
  _cFieldsNames.emplace_back("cellDivergenceField");
  _cFieldsNames.emplace_back("cellMaterialField");
  _cFieldsNames.emplace_back("cellSurfaceMaskField");
  _cFieldsNames.emplace_back("cellSurfaceTDistanceField");
  _cFieldsNames.emplace_back("cellParticleCentroidField");
  _cFieldsNames.emplace_back("cellTestFunction");
  _cFieldsNames.emplace_back("cellTestResult");
  _cFieldsNames.emplace_back("cellTestSolution");
  _fFieldsNames.emplace_back("facePressureField");
  _fFieldsNames.emplace_back("facePressureGradientXField");
  _fFieldsNames.emplace_back("facePressureGradientYField");
  if (D == 3)
    _fFieldsNames.emplace_back("facePressureGradientZField");
  _fFieldsNames.emplace_back("faceXVelocityField");
  _fFieldsNames.emplace_back("faceYVelocityField");
  if (D == 3)
    _fFieldsNames.emplace_back("faceZVelocityField");
  _fFieldsNames.emplace_back("faceMaterialField");
  _fFieldsNames.emplace_back("faceFlowField");
  _fFieldsNames.emplace_back("faceBoundaryTypeField");
  _vFieldsNames.emplace_back("vertexXVelocityField");
  _vFieldsNames.emplace_back("vertexYVelocityField");
  if (D == 3)
    _vFieldsNames.emplace_back("vertexZVelocityField");
  for (auto &name : _cFieldsNames)
    _fieldNames.emplace_back(name);
  for (auto &name : _fFieldsNames)
    _fieldNames.emplace_back(name);
  for (auto &name : _vFieldsNames)
    _fieldNames.emplace_back(name);
  _particles.reset(new ParticleSystem<double, D>(
      new PointZGrid<double, D>(64), new RBFInterpolant<Point<double, D>>{}));
  _particles->setDomainRegion(BBox<double, D>::squareBox());
  // new MLS<D>(Definitions::PolynomialType::QUADRATIC)));
  if (D == 2)
    _graph.reset(dynamic_cast<StructureInterface<D> *>(
        new CellGraph2(BBox2d::squareBox(1.))));
  else
    _graph.reset(dynamic_cast<StructureInterface<D> *>(
        new CellGraph3(BBox3d::squareBox(1.))));
  setupFields();
  for (size_t d = 0; d < 6; d++)
    _domainBoundaries.emplace_back(Definitions::Boundary::NEUMANN);

  _useVerticalFluidFace = [&](size_t faceId) {
    auto &faceMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["faceMaterialField"]);
    if (_graph->faceOrientation(faceId) != Definitions::Orientation::VERTICAL)
      return false;
    return !(faceMaterial[faceId] != Definitions::Material::FLUID);
  };

  _useHorizontalFluidFace = [&](size_t faceId) {
    auto &faceMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["faceMaterialField"]);
    if (_graph->faceOrientation(faceId) != Definitions::Orientation::HORIZONTAL)
      return false;
    return !(faceMaterial[faceId] != Definitions::Material::FLUID);
  };

  _useDepthFluidFace = [&](size_t faceId) {
    auto &faceMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["faceMaterialField"]);
    if (_graph->faceOrientation(faceId) != Definitions::Orientation::DEPTH)
      return false;
    return !(faceMaterial[faceId] != Definitions::Material::FLUID);
  };

  _useFluidFace = [&](size_t faceId) {
    auto &faceMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["faceMaterialField"]);
    return !(faceMaterial[faceId] != Definitions::Material::FLUID);
  };

  _useVelocityFace = [&](size_t faceId) {
    auto &faceMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["faceMaterialField"]);
    return faceMaterial[faceId] != Definitions::Material::SOLID;
  };

  _useVelocityCell = [&](size_t faceId) {
    auto &cellMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["cellMaterialField"]);
    return cellMaterial[faceId] != Definitions::Material::SOLID;
  };

  _useFluidCell = [&](size_t cellId) {
    auto &cellMaterial = *_graph->template field<Definitions::Material>(
        _fieldIds["cellMaterialField"]);
    return !(cellMaterial[cellId] != Definitions::Material::FLUID);
  };

  _cflCallback = [&]() -> double {
    std::vector<size_t> fields;
    fields.emplace_back(_fieldIds["particleXVelocityField"]);
    fields.emplace_back(_fieldIds["particleYVelocityField"]);
    if (D == 3)
      fields.emplace_back(_fieldIds["particleZVelocityField"]);
    return SimulationSteps::cfl<D>(_particles.get(), fields,
                                   1. / (1 << _maxLevel));
  };

  // STEPS ///////////////////////////////////////////////////////////////
  _stepMap["advectParticles"] = [&](double dt) {
    std::vector<size_t> vFields = {_fieldIds["particleXVelocityProperty"],
                                   _fieldIds["particleYVelocityProperty"]};
    if (D == 3)
      vFields.emplace_back(_fieldIds["particleZVelocityProperty"]);
    SimulationSteps::advectParticles<D>(_particles.get(), vFields, _graph.get(),
                                        _solids, &_integrator, dt);
  };
  _stepMap["injectParticles"] = [&](double dt) {
    for (auto &injector : _injectors)
      std::cerr << "injected " << injector.inject(dt) << " particles dt=" << dt
                << "\n";
  };
  _stepMap["handleInjectors"] = [&](double dt) {
    SimulationSteps::handleInjectors(_graph.get(),
                                     _fieldIds["cellMaterialField"],
                                     _fieldIds["cellSurfaceMaskField"],
                                     _maxLevel, _particles.get(), _injectors);
  };
  _stepMap["updateTree"] = [&](double dt) {
    try {
      SimulationSteps::updateTree<D>(
          _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
          _fieldIds["cellSurfaceTDistanceField"],
          _fieldIds["cellMaterialField"], _minLevel, _maxLevel, _graded);
    } catch (const char *e) {
      std::cerr << "##################\n"
                << "save frame " << curFrame() << " " << curSubFrame()
                << std::endl;
      std::ostringstream particlesFile, cellGraphFile, fieldsFile;
      particlesFile << _outputPath << "/_p" << std::setfill('0') << std::setw(6)
                    << curFrame();
      cellGraphFile << _outputPath << "/_c" << std::setfill('0') << std::setw(6)
                    << curFrame();
      fieldsFile << _outputPath << "/_f" << std::setfill('0') << std::setw(6)
                 << curFrame();
      IO::saveFields(fieldsFile.str(), _graph.get());
      IO::saveCellGraph<D>(cellGraphFile.str(), _graph.get());
      IO::saveToPV(particlesFile.str().c_str(), _particles.get());
      throw "failed to update tree!";
    }
  };
  _stepMap["updateTreeCoarseBottom"] = [&](double dt) {
    try {
      SimulationSteps::updateTreeCoarseBottom<D>(
          _particles.get(), _graph.get(), _fieldIds["cellSurfaceMaskField"],
          _fieldIds["cellSurfaceTDistanceField"],
          _fieldIds["cellMaterialField"], _minLevel, _maxLevel, _graded);
    } catch (const char *e) {
      std::cerr << "##################\n"
                << "save frame " << curFrame() << " " << curSubFrame()
                << std::endl;
      std::ostringstream particlesFile, cellGraphFile, fieldsFile;
      particlesFile << _outputPath << "/_p" << std::setfill('0') << std::setw(6)
                    << curFrame();
      cellGraphFile << _outputPath << "/_c" << std::setfill('0') << std::setw(6)
                    << curFrame();
      fieldsFile << _outputPath << "/_f" << std::setfill('0') << std::setw(6)
                 << curFrame();
      IO::saveFields(fieldsFile.str(), _graph.get());
      IO::saveCellGraph<D>(cellGraphFile.str(), _graph.get());
      IO::saveToPV(particlesFile.str().c_str(), _particles.get());
      throw "failed to update tree!";
    }
  };
  _stepMap["markBoundariesOnFaces"] = [&](double dt) {
    SimulationSteps::markFaces<D>(_graph.get(), _fieldIds["cellMaterialField"],
                                  _fieldIds["faceMaterialField"],
                                  _fieldIds["faceBoundaryTypeField"],
                                  _domainBoundaries.data());
  };
  _stepMap["markSolids"] = [&](double dt) {
    SimulationSteps::markSolids<D>(_graph.get(), _fieldIds["cellMaterialField"],
                                   _solids);
  };
  _stepMap["particlesToFaces"] = [&](double dt) { // _USE_OPENMP
    SimulationSteps::transferFromParticlesToGrid<D>(
        _particles.get(), _fieldIds["particleXVelocityProperty"], _graph.get(),
        _fieldIds["faceXVelocityField"], _useVerticalFluidFace, &_interpolator);
    SimulationSteps::transferFromParticlesToGrid<D>(
        _particles.get(), _fieldIds["particleYVelocityProperty"], _graph.get(),
        _fieldIds["faceYVelocityField"], _useHorizontalFluidFace,
        &_interpolator);
    if (D == 3)
      SimulationSteps::transferFromParticlesToGrid<D>(
          _particles.get(), _fieldIds["particleZVelocityProperty"],
          _graph.get(), _fieldIds["faceZVelocityField"], _useDepthFluidFace,
          &_interpolator);
  };
  _stepMap["particlesToFacesHandleLevels"] = [&](double dt) {
    SimulationSteps::transferFromParticlesToGridHandlingCellLevels<D>(
        _particles.get(), _fieldIds["particleXVelocityProperty"], _graph.get(),
        _fieldIds["cellMaterialField"], _fieldIds["faceXVelocityField"],
        _useVerticalFluidFace, &_interpolator);
    SimulationSteps::transferFromParticlesToGridHandlingCellLevels<D>(
        _particles.get(), _fieldIds["particleYVelocityProperty"], _graph.get(),
        _fieldIds["cellMaterialField"], _fieldIds["faceYVelocityField"],
        _useHorizontalFluidFace, &_interpolator);
    if (D == 3)
      SimulationSteps::transferFromParticlesToGridHandlingCellLevels<D>(
          _particles.get(), _fieldIds["particleZVelocityProperty"],
          _graph.get(), _fieldIds["cellMaterialField"],
          _fieldIds["faceZVelocityField"], _useDepthFluidFace, &_interpolator);
  };
  _stepMap["particlesToFacesOneSystem"] = [&](double dt) {

  };
  _stepMap["particlesToMidPointFaces"] = [&](double dt) {
    std::vector<size_t> faceVelocityFieldIds = {
        _fieldIds["faceXVelocityField"], _fieldIds["faceYVelocityField"]};
    std::vector<size_t> particleVelocityFieldIds = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3) {
      faceVelocityFieldIds.emplace_back(_fieldIds["faceZVelocityField"]);
      particleVelocityFieldIds.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    }
    SimulationSteps::transferVelocitiesFromParticlesToMidPointFaces(
        _graph.get(), _particles.get(), particleVelocityFieldIds,
        faceVelocityFieldIds, _fieldIds["faceMaterialField"],
        _fieldIds["cellMaterialField"], _fieldIds["faceMidPointField"]);
  };
  _stepMap["divergenceFromMidPointFaces"] = [&](double dt) {
    std::vector<size_t> faceVelocityFieldIds = {
        _fieldIds["faceXVelocityField"], _fieldIds["faceYVelocityField"]};
    if (D == 3)
      faceVelocityFieldIds.emplace_back(_fieldIds["faceZVelocityField"]);
    SimulationSteps::computeDivergenceFieldFromMidPointFaces(
        _graph.get(), _fieldIds["cellDivergenceField"],
        _fieldIds["faceMidPointField"], faceVelocityFieldIds, _useFluidCell,
        _rbf.get());
  };
  _stepMap["particlesToCells"] = [&](double dt) {
    SimulationSteps::transferFromParticlesToGrid<D>(
        _particles.get(), _fieldIds["particleXVelocityProperty"], 20,
        _graph.get(), _fieldIds["cellXVelocityField"], _useFluidCell);
    SimulationSteps::transferFromParticlesToGrid<D>(
        _particles.get(), _fieldIds["particleYVelocityProperty"], 20,
        _graph.get(), _fieldIds["cellYVelocityField"], _useFluidCell);
    if (D == 3)
      SimulationSteps::transferFromParticlesToGrid<D>(
          _particles.get(), _fieldIds["particleZVelocityProperty"], 20,
          _graph.get(), _fieldIds["cellZVelocityField"], _useFluidCell);
  };
  _stepMap["addGravityToFaces"] = [&](double dt) {
    SimulationSteps::applyExternalForce<D>(_graph.get(),
                                           _fieldIds["faceYVelocityField"],
                                           -9.81, dt, _useHorizontalFluidFace);
  };
  _stepMap["addGravityToMidPointFaces"] = [&](double dt) {
    std::vector<size_t> faceFields = {_fieldIds["faceXVelocityField"],
                                      _fieldIds["faceYVelocityField"]};
    if (D == 3)
      faceFields.emplace_back(_fieldIds["faceZVelocityField"]);
    Vector<double, D> acceleration;
    acceleration[1] = -9.81;
    SimulationSteps::applyExternalForceToMidPointFaces<D>(
        _graph.get(), faceFields, acceleration, dt,
        _fieldIds["faceMaterialField"], _fieldIds["cellMaterialField"]);
  };
  _stepMap["addGravityToParticles"] = [&](double dt) {
    SimulationSteps::applyExternalForce<D>(
        _particles.get(), _fieldIds["particleYVelocityProperty"], -9.81, dt);
  };
  _stepMap["addGravityToCells"] = [&](double dt) {
    SimulationSteps::applyExternalForce<D>(_graph.get(),
                                           _fieldIds["cellYVelocityField"],
                                           -9.81, dt, _useFluidCell);
  };
  _stepMap["divergenceOnCellsFromFaces"] = [&](double dt) {
    std::vector<size_t> faceFields = {_fieldIds["faceXVelocityField"],
                                      _fieldIds["faceYVelocityField"]};
    if (D == 3)
      faceFields.emplace_back(_fieldIds["faceZVelocityField"]);
    SimulationSteps::computeDivergenceField<D>(
        _graph.get(), _fieldIds["cellDivergenceField"], faceFields,
        _useFluidCell, _rbf.get(), _adaptiveRBFRadius);
  };

  _stepMap["divergenceOnCellsFromFacesAndCellCenter"] = [&](double dt) {
    std::vector<size_t> faceFields = {_fieldIds["faceXVelocityField"],
                                      _fieldIds["faceYVelocityField"]};
    std::vector<size_t> cellFields = {_fieldIds["cellXVelocityField"],
                                      _fieldIds["cellYVelocityField"]};
    if (D == 3) {
      faceFields.emplace_back(_fieldIds["faceZVelocityField"]);
      cellFields.emplace_back(_fieldIds["cellZVelocityField"]);
    }
    SimulationSteps::computeDivergenceField<D>(
        _graph.get(), _fieldIds["cellDivergenceField"], faceFields, cellFields,
        _useFluidCell, _rbf.get(), _adaptiveRBFRadius);
  };
  _stepMap["applyFalloffToParticles"] = [&](double dt) {
    std::vector<size_t> particleVelocityFields = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3)
      particleVelocityFields.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    SimulationSteps::applyFalloffToParticles<D>(_graph.get(), _particles.get(),
                                                _fieldIds["cellMaterialField"],
                                                particleVelocityFields);
  };
  _stepMap["divergenceOnParticles"] = [&](double dt) {
    std::vector<size_t> particleVelocityFields = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3)
      particleVelocityFields.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    SimulationSteps::computeDivergenceOnParticles<D>(
        _particles.get(), particleVelocityFields,
        _fieldIds["particleDivergenceField"], _rbf.get());
  };
  _stepMap["divergenceOnCellsFromParticles"] = [&](double dt) {
    std::vector<size_t> particleVelocityFields = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3)
      particleVelocityFields.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    SimulationSteps::computeDivergenceFieldFromParticles<D>(
        _graph.get(), _particles.get(), _fieldIds["cellDivergenceField"],
        _fieldIds["cellMaterialField"], _fieldIds["cellParticleCentroidField"],
        particleVelocityFields, _useFluidCell, _rbf.get());
  };
  _stepMap["pressureOnCellsRing"] = [&](double dt) {
    SimulationSteps::solvePressureRing<D>(
        _graph.get(), _fieldIds["cellPressureField"],
        _fieldIds["cellDivergenceField"], _fieldIds["cellMaterialField"], dt,
        _rbf.get());
  };

  _stepMap["pressureOnCells"] = [&](double dt) {
    SimulationSteps::solvePressure<D>(
        _graph.get(), _fieldIds["cellPressureField"],
        _fieldIds["cellDivergenceField"], _fieldIds["cellMaterialField"], dt,
        _rbf.get(), _adaptiveRBFRadius);
  };
  _stepMap["pressureOnCellParticleCenter"] = [&](double dt) {
    SimulationSteps::solvePressure<D>(
        _graph.get(), _particles.get(), _fieldIds["cellPressureField"],
        _fieldIds["cellDivergenceField"], _fieldIds["cellMaterialField"],
        _fieldIds["cellSurfaceMaskField"],
        _fieldIds["cellParticleCentroidField"], dt, _rbf.get(),
        _adaptiveRBFRadius);
  };
  _stepMap["pressureGradientOnMidPointFaces"] = [&](double dt) {
    std::vector<size_t> facePressureFields = {
        _fieldIds["facePressureGradientXField"],
        _fieldIds["facePressureGradientYField"]};
    if (D == 3)
      facePressureFields.emplace_back(_fieldIds["facePressureGradientZField"]);
    SimulationSteps::computePressureGradientOnMidPointFaces(
        _graph.get(), _fieldIds["cellMaterialField"],
        _fieldIds["cellPressureField"], _fieldIds["faceMaterialField"],
        facePressureFields, _rbf.get());
  };
  _stepMap["pressureGradientOnCellsParticleCenter"] = [&](double dt) {
    std::vector<size_t> cellFields = {_fieldIds["cellPressureGradientXField"],
                                      _fieldIds["cellPressureGradientYField"]};
    if (D == 3)
      cellFields.emplace_back(_fieldIds["cellPressureGradientZField"]);
    SimulationSteps::computePressureGradient(
        _graph.get(), _particles.get(), _fieldIds["cellPressureField"],
        cellFields, _fieldIds["cellMaterialField"],
        _fieldIds["cellParticleCentroidField"], _rbf.get());
  };
  _stepMap["pressureGradientOnCells"] = [&](double dt) {
    std::vector<size_t> cellFields = {_fieldIds["cellPressureGradientXField"],
                                      _fieldIds["cellPressureGradientYField"]};
    if (D == 3)
      cellFields.emplace_back(_fieldIds["cellPressureGradientZField"]);
    SimulationSteps::computePressureGradient(
        _graph.get(), _fieldIds["cellPressureField"], cellFields,
        _fieldIds["cellMaterialField"], _rbf.get());
  };
  _stepMap["cellPressureToFacePressure"] = [&](double dt) {
    SimulationSteps::transferCellFieldToFaceField<D>(
        _graph.get(), _fieldIds["cellPressureField"],
        _fieldIds["facePressureField"], _useFluidCell, _useFluidFace);
  };
  _stepMap["cellPressureGradientToFacePressureGradient"] = [&](double dt) {
    SimulationSteps::transferCellFieldToFaceField<D>(
        _graph.get(), _fieldIds["cellPressureGradientXField"],
        _fieldIds["facePressureGradientXField"], _useFluidCell,
        _useVerticalFluidFace);
    SimulationSteps::transferCellFieldToFaceField<D>(
        _graph.get(), _fieldIds["cellPressureGradientYField"],
        _fieldIds["facePressureGradientYField"], _useFluidCell,
        _useHorizontalFluidFace);
    if (D == 3)
      SimulationSteps::transferCellFieldToFaceField<D>(
          _graph.get(), _fieldIds["cellPressureGradientZField"],
          _fieldIds["facePressureGradientZField"], _useFluidCell,
          _useDepthFluidFace);
  };
  _stepMap["pressureGradientOnFacesFromCellsAndFaces"] = [&](double dt) {
    std::vector<size_t> faceFields = {_fieldIds["facePressureGradientXField"],
                                      _fieldIds["facePressureGradientYField"]};
    if (D == 3)
      faceFields.emplace_back(_fieldIds["facePressureGradientZField"]);
    SimulationSteps::computePressureGradient<D>(
        _graph.get(), _fieldIds["cellPressureField"],
        _fieldIds["facePressureField"], faceFields,
        _fieldIds["faceMaterialField"], _rbf.get());
  };
  _stepMap["pressureGradientOnFacesFromCells"] = [&](double dt) {
    std::vector<size_t> faceFields = {_fieldIds["facePressureGradientXField"],
                                      _fieldIds["facePressureGradientYField"]};
    if (D == 3)
      faceFields.emplace_back(_fieldIds["facePressureGradientZField"]);
    SimulationSteps::computePressureGradient<D>(
        _graph.get(), _fieldIds["faceMaterialField"],
        _fieldIds["cellPressureField"], _fieldIds["cellMaterialField"],
        faceFields, _rbf.get());
  };
  _stepMap["cellVelocitiesToParticles"] = [&](double dt) {
    std::vector<size_t> velocityFields = {_fieldIds["cellXVelocityField"],
                                          _fieldIds["cellYVelocityField"]};
    std::vector<size_t> particleVelocityFields = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3) {
      velocityFields.emplace_back(_fieldIds["cellZVelocityField"]);
      particleVelocityFields.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    }
    SimulationSteps::transferCellVelocitiesToParticles(
        _graph.get(), _particles.get(), _fieldIds["cellMaterialField"],
        velocityFields, particleVelocityFields, &_interpolator);
  };
  _stepMap["correctFaceVelocitiesFromFacePressureGradient"] = [&](double dt) {
    std::vector<size_t> velocityFields = {_fieldIds["faceXVelocityField"],
                                          _fieldIds["faceYVelocityField"]};
    std::vector<size_t> pgFields = {_fieldIds["facePressureGradientXField"],
                                    _fieldIds["facePressureGradientYField"]};
    if (D == 3) {
      velocityFields.emplace_back(_fieldIds["faceZVelocityField"]);
      pgFields.emplace_back(_fieldIds["facePressureGradientZField"]);
    }
    SimulationSteps::correctVelocities<D>(_graph.get(), pgFields,
                                          velocityFields,
                                          _fieldIds["faceMaterialField"], dt);
  };
  _stepMap["correctMidPointFaceVelocitiesFromFacePressureGradient"] =
      [&](double dt) {
        std::vector<size_t> velocityFields = {_fieldIds["faceXVelocityField"],
                                              _fieldIds["faceYVelocityField"]};
        std::vector<size_t> pgFields = {
            _fieldIds["facePressureGradientXField"],
            _fieldIds["facePressureGradientYField"]};
        if (D == 3) {
          velocityFields.emplace_back(_fieldIds["faceZVelocityField"]);
          pgFields.emplace_back(_fieldIds["facePressureGradientZField"]);
        }
        SimulationSteps::correctVelocitiesOnMidPointFaces<D>(
            _graph.get(), pgFields, velocityFields,
            _fieldIds["faceMaterialField"], dt);
      };
  _stepMap["correctVelocityOnCellFromFaceFlow"] = [&](double dt) {
    std::vector<size_t> faceVelocities = {_fieldIds["faceXVelocityField"],
                                          _fieldIds["faceYVelocityField"]};
    std::vector<size_t> cellVelocities = {_fieldIds["cellXVelocityField"],
                                          _fieldIds["cellYVelocityField"]};
    if (D == 3) {
      faceVelocities.emplace_back(_fieldIds["faceZVelocityField"]);
      cellVelocities.emplace_back(_fieldIds["cellZVelocityField"]);
    }

    SimulationSteps::transferMidPointFacesFlowToCellVelocities<D>(
        _graph.get(), _fieldIds["cellPressureField"],
        _fieldIds["cellMaterialField"], _fieldIds["faceMaterialField"],
        _fieldIds["faceMidPointField"], cellVelocities, faceVelocities,
        _rbf.get(), dt);
  };
  _stepMap["correctParticleVelocitiesFromFacePressureGradient"] =
      [&](double dt) {
        std::vector<size_t> particleVelocityFields = {
            _fieldIds["particleXVelocityProperty"],
            _fieldIds["particleYVelocityProperty"]};
        std::vector<size_t> particlePressureGradients = {
            _fieldIds["particleXPressureGradientProperty"],
            _fieldIds["particleYPressureGradientProperty"]};
        SimulationSteps::transferFromGridToParticles<D>(
            _graph.get(), _fieldIds["facePressureGradientXField"],
            _particles.get(), particlePressureGradients[0], _useVelocityFace,
            _useVelocityCell,
            Definitions::BoundaryVelocity::NO_SLIP_HORIZONTAL);
        SimulationSteps::transferFromGridToParticles<D>(
            _graph.get(), _fieldIds["facePressureGradientYField"],
            _particles.get(), particlePressureGradients[1], _useVelocityFace,
            _useVelocityCell, Definitions::BoundaryVelocity::NO_SLIP_VERTICAL);
        if (D == 3)
          SimulationSteps::transferFromGridToParticles<D>(
              _graph.get(), _fieldIds["facePressureGradientZField"],
              _particles.get(), particlePressureGradients[2], _useVelocityFace,
              _useVelocityCell, Definitions::BoundaryVelocity::NO_SLIP_DEPTH);

        SimulationSteps::correctVelocities(_particles.get(),
                                           particlePressureGradients,
                                           particleVelocityFields, dt);
      };
  _stepMap["correctParticleVelocitiesFromCellPressureGradient"] =
      [&](double dt) {
        SimulationSteps::correctVelocities<D>(
            _graph.get(), _fieldIds["cellPressureGradientXField"],
            _particles.get(), _fieldIds["particleXVelocityProperty"],
            _useFluidCell, dt);
        SimulationSteps::correctVelocities<D>(
            _graph.get(), _fieldIds["cellPressureGradientYField"],
            _particles.get(), _fieldIds["particleYVelocityProperty"],
            _useFluidCell, dt);
        if (D == 3)
          SimulationSteps::correctVelocities<D>(
              _graph.get(), _fieldIds["cellPressureGradientZField"],
              _particles.get(), _fieldIds["particleZVelocityProperty"],
              _useFluidCell, dt);
      };
  _stepMap["correctParticleVelocitiesFromCellPressureGradientParticleCenter"] =
      [&](double dt) {
        std::vector<size_t> pressureGradientFields = {
            _fieldIds["cellPressureGradientXField"],
            _fieldIds["cellPressureGradientYField"]};
        std::vector<size_t> particleVelocityFields = {
            _fieldIds["particleXVelocityProperty"],
            _fieldIds["particleYVelocityProperty"]};
        if (D == 3) {
          pressureGradientFields.emplace_back(
              _fieldIds["cellPressureGradientZField"]);
          particleVelocityFields.emplace_back(
              _fieldIds["particleZVelocityProperty"]);
        }
        SimulationSteps::correctVelocities<D>(
            _graph.get(), _particles.get(), pressureGradientFields,
            _fieldIds["cellSurfaceMaskField"], _fieldIds["cellMaterialField"],
            _fieldIds["cellParticleCentroidField"], particleVelocityFields, dt,
            &_interpolator);
        // SimulationSteps::correctVelocitiesOnSurface<D>(
        //     _graph.get(), _particles.get(), pressureGradientFields,
        //     _fieldIds["cellSurfaceMaskField"],
        //     _fieldIds["cellMaterialField"],
        //     _fieldIds["cellParticleCentroidField"], particleVelocityFields,
        //     dt,
        //     &_interpolator);
      };
  _stepMap["correctParticleVelocitiesFromParticlePressureGradient"] =
      [&](double dt) {
        std::vector<size_t> particleVelocityFields = {
            _fieldIds["particleXVelocityProperty"],
            _fieldIds["particleYVelocityProperty"]};
        std::vector<size_t> particlePressureGradients = {
            _fieldIds["particleXPressureGradientProperty"],
            _fieldIds["particleYPressureGradientProperty"]};
        if (D == 3) {
          particlePressureGradients.emplace_back(
              _fieldIds["particleZPressureGradientProperty"]);
          particleVelocityFields.emplace_back(
              _fieldIds["particleZVelocityProperty"]);
        }
        SimulationSteps::computePressureGradient(
            _graph.get(), _particles.get(), _fieldIds["cellMaterialField"],
            _fieldIds["cellPressureField"],
            _fieldIds["cellParticleCentroidField"], particlePressureGradients,
            _rbf.get());
        SimulationSteps::correctVelocities(_particles.get(),
                                           particlePressureGradients,
                                           particleVelocityFields, dt);
      };
  _stepMap["computeCellParticlesCentroid"] = [&](double dt) {
    SimulationSteps::computeCellParticlesCentroid<D>(
        _graph.get(), _particles.get(), _fieldIds["cellSurfaceMaskField"],
        _fieldIds["cellMaterialField"], _fieldIds["cellParticleCentroidField"]);
  };
  _stepMap["computeFaceMidPoint"] = [&](double dt) {
    SimulationSteps::computeFaceMidPoint(_graph.get(),
                                         _fieldIds["faceMidPointField"],
                                         _fieldIds["faceMaterialField"]);
  };
  _stepMap["propagateVelocities"] = [&](double dt) {
    std::vector<size_t> velocityFields = {_fieldIds["faceXVelocityField"],
                                          _fieldIds["faceYVelocityField"]};
    if (D == 3)
      velocityFields.emplace_back(_fieldIds["faceZVelocityField"]);
    SimulationSteps::propagateVelocities<D>(
        _graph.get(), velocityFields, _fieldIds["cellMaterialField"],
        _fieldIds["faceMaterialField"], _fieldIds["cellSurfaceMaskField"]);
  };
  _stepMap["reseedParticles"] = [&](double dt) {
    SimulationSteps::reseedParticles<D>(
        _graph.get(), _fieldIds["cellMaterialField"],
        _fieldIds["cellSurfaceMaskField"], _particles.get());
  };
  _stepMap["reseedParticlesWithVelocities"] = [&](double dt) {
    std::vector<size_t> particleVelocities = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"],
    };
    if (D == 3)
      particleVelocities.emplace_back(_fieldIds["particleZVelocityProperty"]);
    SimulationSteps::reseedParticles<D>(_graph.get(),
                                        _fieldIds["cellMaterialField"],
                                        _fieldIds["cellSurfaceMaskField"],
                                        _particles.get(), particleVelocities);
  };
  _stepMap["facesToParticlesAllComponents"] = [&](double dt) { // _USE_OPENMP
    std::vector<size_t> faceVelocityFieldIds = {
        _fieldIds["faceXVelocityField"], _fieldIds["faceYVelocityField"]};
    std::vector<size_t> particleVelocities = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"],
    };
    if (D == 3) {
      faceVelocityFieldIds.emplace_back(_fieldIds["faceZVelocityField"]);
      particleVelocities.emplace_back(_fieldIds["particleZVelocityProperty"]);
    }
    SimulationSteps::transferFromGridToParticles<D>(
        _graph.get(), _particles.get(), _fieldIds["cellMaterialField"],
        _fieldIds["faceMaterialField"], _fieldIds["cellSurfaceMaskField"],
        faceVelocityFieldIds, particleVelocities, &_interpolator);
  };
  _stepMap["facesToParticles"] = [&](double dt) {
    SimulationSteps::transferFromGridToParticles<D>(
        _graph.get(), _fieldIds["faceXVelocityField"], _particles.get(),
        _fieldIds["particleXVelocityProperty"], _useVelocityFace,
        _useVelocityCell, Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL);
    SimulationSteps::transferFromGridToParticles<D>(
        _graph.get(), _fieldIds["faceYVelocityField"], _particles.get(),
        _fieldIds["particleYVelocityProperty"], _useVelocityFace,
        _useVelocityCell, Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL);
    if (D == 3)
      SimulationSteps::transferFromGridToParticles<D>(
          _graph.get(), _fieldIds["faceZVelocityField"], _particles.get(),
          _fieldIds["particleZVelocityProperty"], _useVelocityFace,
          _useVelocityCell, Definitions::BoundaryVelocity::FREE_SLIP_DEPTH);
  };
  _stepMap["midPointFacesToParticles"] = [&](double dt) {
    std::vector<size_t> faceVelocityFieldIds = {
        _fieldIds["faceXVelocityField"], _fieldIds["faceYVelocityField"]};
    std::vector<size_t> particleVelocityFieldIds = {
        _fieldIds["particleXVelocityProperty"],
        _fieldIds["particleYVelocityProperty"]};
    if (D == 3) {
      faceVelocityFieldIds.emplace_back(_fieldIds["faceZVelocityField"]);
      particleVelocityFieldIds.emplace_back(
          _fieldIds["particleZVelocityProperty"]);
    }
    SimulationSteps::transferMidPointFaceVelocitiesToParticles(
        _graph.get(), _particles.get(), _fieldIds["cellMaterialField"],
        _fieldIds["faceMaterialField"], faceVelocityFieldIds,
        particleVelocityFieldIds, Definitions::BoundaryVelocity::FREE_SLIP,
        &_interpolator);
  };
  _stepMap["saveAllData"] = [&](double dt) {
    std::cerr << "save frame " << curFrame() << " " << curSubFrame()
              << std::endl;
    std::ostringstream particlesFile, cellGraphFile, fieldsFile;
    particlesFile << _outputPath << "/p" << std::setfill('0') << std::setw(6)
                  << curFrame();
    cellGraphFile << _outputPath << "/c" << std::setfill('0') << std::setw(6)
                  << curFrame();
    fieldsFile << _outputPath << "/f" << std::setfill('0') << std::setw(6)
               << curFrame();
    IO::saveFields(fieldsFile.str(), _graph.get());
    IO::saveCellGraph<D>(cellGraphFile.str(), _graph.get());
    IO::saveToPV(particlesFile.str().c_str(), _particles.get());
  };
  _stepMap["testGradX"] = [&](double dt) {
    SimulationSteps::computeGradient<D>(
        0, _graph.get(), _fieldIds["cellTestFunction"],
        _fieldIds["cellTestResult"], _fieldIds["cellMaterialField"],
        _fieldIds["cellParticleCentroidField"], _rbf.get());
    SimulationSteps::computeError<D>(_graph.get(), _fieldIds["cellTestResult"],
                                     _fieldIds["cellTestSolution"]);
  };
  _stepMap["testGradY"] = [&](double dt) {
    SimulationSteps::computeGradient<D>(
        1, _graph.get(), _fieldIds["cellTestFunction"],
        _fieldIds["cellTestResult"], _fieldIds["cellMaterialField"],
        _fieldIds["cellParticleCentroidField"], _rbf.get());
    SimulationSteps::computeError<D>(_graph.get(), _fieldIds["cellTestResult"],
                                     _fieldIds["cellTestSolution"]);
  };
  _stepMap["testLaplacian"] = [&](double dt) {
    SimulationSteps::computeLaplacian<D>(
        _graph.get(), _fieldIds["cellTestFunction"],
        _fieldIds["cellTestResult"], _fieldIds["cellMaterialField"],
        _fieldIds["cellParticleCentroidField"], _rbf.get());
  };
}

template <int D> double CustomSimulation<D>::timeStep() const { return _dt; }
