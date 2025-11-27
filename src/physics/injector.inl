template <int D> void ParticleInjector<D>::setSpacing(double s) {
  _spacing = s;
}

template <int D>
size_t ParticleInjector<D>::setupBoxShape(BBox<double, D> region, int group, int groupPropertyId) {
  Point<double, D> position = region.lower();
  position += Vector<double, D>(_spacing / 2);
  size_t total = 0;
  while (position <= region.upper()) {
    auto id = _particleSystem->addParticle(position);
    _particleSystem->setScalarProperty(groupPropertyId, id, group);
    total++;
    position[0] += _spacing;
    for (int i = 0; i < D - 1; i++)
      if (position[i] > region.upper()[i]) {
        position[i] = region.lower()[i] + _spacing / 2;
        position[i + 1] += _spacing;
      }
  }
  return total;
}

template <int D> void ParticleInjector<D>::fromFile(std::string filename) {
  IO::loadFromPV(filename, _particleSystem);
}

template <int D>
size_t ParticleInjector<D>::setupCircleShape(Point<double, D> center,
                                             double radius) {
  Vector<double, D> delta(radius);
  BBox<double, D> region(center - radius, center + radius);
  Point<double, D> position = region.lower();
  size_t total = 0;
  while (position <= region.upper()) {
    if (distance2(position, center) <= radius * radius) {
      _particleSystem->addParticle(position);
      total++;
    }
    position[0] += _spacing;
    for (int i = 0; i < D - 1; i++)
      if (position[i] > region.upper()[i]) {
        position[i] = region.lower()[i];
        position[i + 1] += _spacing;
      }
  }
  return total;
}

template <int D> void ParticleInjector<D>::loadFromFile(const char *filename) {
  IO::loadFromPV(filename, _particleSystem);
}

template <int D>
ParticleInjector<D>::ParticleInjector(BBox<double, D> region,
                                      Transform<D> transform,
                                      Vector<double, D> velocity,
                                      double inflowRate,
                                      ParticleSystem<double, D> *particleSystem)
    : _injectionRegion(region), _transform(transform), _velocity(velocity),
      _inflowRate(inflowRate), _particleSystem(particleSystem),
      _elapsedTime(0.0), _startTime(0.0), _endTime(INFINITY) {
  _worldSpaceRegion = _transform(_injectionRegion);
}

template <int D>
ParticleInjector<D>::ParticleInjector(ParticleSystem<double, D> *particleSystem,
                                      double spacing)
    : _spacing(spacing), _particleSystem(particleSystem) {}

template <int D> size_t ParticleInjector<D>::inject(double dt) {
  this->_elapsedTime = this->_elapsedTime + dt;
  if (_elapsedTime >= _startTime && _elapsedTime <= _endTime) {
    size_t n = static_cast<size_t>(_inflowRate * dt);
    static BoxSampler sampler;
    for (size_t i = 0; i < n; i++) {
      auto pos = sampler.sample(_injectionRegion);
      size_t index = _particleSystem->addParticle(_transform(pos));
      for (int d = 0; d < D; d++)
        _particleSystem->setScalarProperty(d, index, _velocity[d]);
    }
    return n;
  }
  return 0;
}

template <int D> BBox<double, D> ParticleInjector<D>::region() const {
  return _worldSpaceRegion;
}

template <int D>
void ParticleInjector<D>::setStartEndTime(double startTime, double endTime) {
  _startTime = startTime;
  _endTime = endTime;
}
