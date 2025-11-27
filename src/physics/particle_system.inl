template <typename T, int D> ParticleSystem<T, D>::ParticleSystem() = default;

template <typename T, int D>
ParticleSystem<T, D>::ParticleSystem(PointSetInterface<T, D> *pointSet,
                                     Interpolant<Point<T, D>> *interpolant) {
  _pointSet.reset(pointSet);
  _interpolant.reset(interpolant);
}

template <typename T, int D> ParticleSystem<T, D>::~ParticleSystem() = default;

template <typename T, int D>
double ParticleSystem<T, D>::sampleScalarProperty(size_t propertyId,
                                                  const Point<T, D> &p,
                                                  size_t n) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  std::vector<Point<T, D>> points;
  std::vector<double> values;
  //  std::cerr << "sampling on target " << p << std::endl;
  _pointSet->iterateClosestPoints(p, n,
                                  [&](size_t i, Point<T, D> pos) {
                                    //    std::cerr << i << " ";
                                    //    std::cerr << pos << std::endl;
                                    points.emplace_back(pos);
                                    values.emplace_back(
                                        _scalarProperties[propertyId][i]);
                                  },
                                  nullptr);
  return _interpolant->interpolateAt(p, points, values);
}

template <typename T, int D>
Vector<double, D> ParticleSystem<T, D>::sampleVectorFromProperties(
    const std::vector<size_t> &propertyIds, const Point<T, D> &target,
    size_t n) const {
  THROW(propertyIds.size() == D, "ParticleSystem::sampleVectorFromProperties "
                                 "wrong number of property ids.");
  std::vector<Point<T, D>> points;
  std::vector<size_t> pointIds;
  _pointSet->iterateClosestPoints(target, n,
                                  [&](size_t i, Point<T, D> pos) {
                                    points.emplace_back(pos);
                                    pointIds.emplace_back(i);
                                  },
                                  nullptr);
  std::vector<Point<double, D>> values(points.size(), 0.);
  for (int d = 0; d < D; d++) {
    THROW(propertyIds[d] < _scalarProperties.size(),
          "ParticleSystem::sampleVectorFromProperties invalid property id");
    for (size_t i = 0; i < points.size(); i++)
      values[i][d] = _scalarProperties[propertyIds[d]][pointIds[i]];
  }
  auto ans = _interpolant->interpolateAt(target, points, values);
  return Vector<double, D>(D, ans.asConstArray());
}

template <typename T, int D>
size_t ParticleSystem<T, D>::addScalarProperty(double v) {
  _scalarProperties.emplace_back(_pointSet->size(), v);
  return _scalarProperties.size() - 1;
}

template <typename T, int D>
size_t ParticleSystem<T, D>::addParticle(Point<T, D> position) {
  size_t id = _pointSet->add(position);
  for (size_t i = 0; i < _scalarProperties.size(); i++)
    if (_scalarProperties[i].size() < id + 1)
      _scalarProperties[i].resize(id + 1, 0.);
  return id;
}

template <typename T, int D>
void ParticleSystem<T, D>::setDomainRegion(const BBox<T, D> &region) {
  _pointSet->setDomainRegion(region);
}

template <typename T, int D>
void ParticleSystem<T, D>::setScalarProperty(size_t propertyId,
                                             size_t particleId, double value) {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  if (particleId >= _scalarProperties[propertyId].size())
    std::cerr << particleId << " " << _scalarProperties[propertyId].size()
              << std::endl;
  ASSERT_FATAL(particleId < _scalarProperties[propertyId].size());
  _scalarProperties[propertyId][particleId] = value;
}

template <typename T, int D>
double ParticleSystem<T, D>::getScalarProperty(size_t propertyId,
                                               size_t particleId) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  ASSERT_FATAL(particleId < _scalarProperties[propertyId].size());
  return _scalarProperties[propertyId][particleId];
}

template <typename T, int D>
Vector<T, D>
ParticleSystem<T, D>::getVector(const std::vector<size_t> &propertyIds,
                                size_t particleId) const {
  std::vector<T> v;
  for (size_t d = 0; d < D; d++) {
    v.emplace_back(_scalarProperties[propertyIds[d]][particleId]);
  }
  return Vector<double, D>(D, v.data());
}

template <typename T, int D>
Point<T, D> ParticleSystem<T, D>::operator[](size_t id) const {
  return (*_pointSet)[id];
}

template <typename T, int D>
double ParticleSystem<T, D>::maxPropertyValue(size_t propertyId) const {
  // TODO we should call update instead of size!
  _pointSet->size();
  double m = -INFINITY;
  _pointSet->iteratePoints([&](size_t i, Point<T, D>) {
    m = std::max(m, _scalarProperties[propertyId][i]);
  });
  return m;
}

template <typename T, int D> inline void ParticleSystem<T, D>::update() {
  _pointSet->update();
}

template <typename T, int D>
double ParticleSystem<T, D>::minPropertyValue(size_t propertyId) const {
  // TODO we should call update instead of size!
  _pointSet->size();
  double m = INFINITY;
  _pointSet->iteratePoints([&](size_t i, Point<T, D>) {
    m = std::min(m, _scalarProperties[propertyId][i]);
  });
  return m;
}

template <typename T, int D> void ParticleSystem<T, D>::clear() {
  _pointSet->clear();
}

template <typename T, int D> size_t ParticleSystem<T, D>::size() {
  return _pointSet->size();
}

template <typename T, int D>
void ParticleSystem<T, D>::setPosition(size_t particleId,
                                       Point<T, D> position) {
  _pointSet->setPosition(particleId, position);
}

template <typename T, int D>
void ParticleSystem<T, D>::iterateParticles(
    const BBox<T, D> &searchRegion,
    const std::function<void(size_t, Point<T, D>)> &f) const {
#ifdef FUROO_PROFILE
  CodeProfiler::count("boxSearchCount");
#endif
  _pointSet->search(searchRegion, [&](size_t i) { f(i, (*_pointSet)[i]); });
}

template <typename T, int D>
void ParticleSystem<T, D>::iterateParticlesZ(
    const BBox<T, D> &b, size_t level,
    const std::function<void(size_t, Point<T, D>)> &f) const {
  _pointSet->searchZ(b, level, [&](size_t i) { f(i, (*_pointSet)[i]); });
}

template <typename T, int D>
bool ParticleSystem<T, D>::containsParticle(
    const BBox<T, D> &searchRegion) const {
  return _pointSet->containsPoint(searchRegion);
}

template <typename T, int D>
bool ParticleSystem<T, D>::containsParticleZ(
    const BBox<T, D> &searchRegion, size_t level) const {
  return _pointSet->containsPointZ(searchRegion, level);
}

template <typename T, int D>
void ParticleSystem<T, D>::iterateClosestParticles(
    Point<T, D> center, size_t neighborCount,
    const std::function<void(size_t, Point<T, D>)> &f) const {
#ifdef FUROO_PROFILE
  CodeProfiler::count("knnSearchCount");
#endif
  _pointSet->iterateClosestPoints(center, neighborCount, f, nullptr);
}

template <typename T, int D>
double
ParticleSystem<T, D>::sampleScalarProperty(size_t propertyId, double radius,
                                           const Point<T, D> &target) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  std::vector<Point<T, D>> points;
  std::vector<double> values;
  _pointSet->search(BBox<T, D>(target, radius), [&](size_t i) {
    points.emplace_back((*_pointSet)[i]);
    values.emplace_back(_scalarProperties[propertyId][i]);
  });
  return _interpolant->interpolateAt(target, points, values);
}

template <typename T, int D>
Vector<double, D> ParticleSystem<T, D>::sampleVectorFromProperties(
    const std::vector<size_t> &propertyIds, double radius,
    const Point<T, D> &target) const {
  THROW(propertyIds.size() == D, "ParticleSystem::sampleVectorFromProperties "
                                 "wrong number of property ids.");
  std::vector<Point<T, D>> points;
  std::vector<size_t> pointIds;
  _pointSet->search(BBox<T, D>(target, radius), [&](size_t i) {
    points.emplace_back((*_pointSet)[i]);
    pointIds.emplace_back(i);
  });
  Vector<double, D> v;
  std::vector<double> values(points.size(), 0.);
  for (int d = 0; d < D; d++) {
    THROW(propertyIds[d] < _scalarProperties.size(),
          "ParticleSystem::sampleVectorFromProperties invalid property id");
    for (size_t i = 0; i < points.size(); i++)
      values[i] = _scalarProperties[propertyIds[d]][pointIds[i]];
    v[d] = _interpolant->interpolateAt(target, points, values);
  }
  return v;
}

template <typename T, int D>
void ParticleSystem<T, D>::iterateParticles(
    const std::function<void(size_t, Point<T, D>)> &f) const {
  // TODO we should call update instead of size!
  _pointSet->size();
  _pointSet->iteratePoints(f);
}

#ifdef _USE_OPENMP
template <typename T, int D>
inline void
ParticleSystem<T, D>::iterateParticles_par(
  const std::function<void(size_t, Point<T, D>)> &f) const
{
  _pointSet->iteratePoints(f);
}
#endif // _USE_OPENMP

template <typename T, int D> void ParticleSystem<T, D>::remove(size_t i) {
  _pointSet->remove(i);
}
