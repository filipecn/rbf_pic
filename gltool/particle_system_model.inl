template <int D> ParticleSystemModel<D>::ParticleSystemModel() {
  _particlesSet.reset(new GraphicElementSet(
      ponos::create_icosphere_mesh(ponos::Point3(), 1.f, 3, false, false)));
  _velocitiesSet.reset(new GraphicElementSet(createArrow<D>()));
}

template <int D>
ParticleSystemModel<D>::ParticleSystemModel(
    furoo::ParticleSystem<double, D> *ps, double radius) {
  set(ps, radius);
}

template <int D>
void ParticleSystemModel<D>::set(furoo::ParticleSystem<double, D> *ps,
                                 double radius) {
  this->ps = ps;
  particleRadius = radius;
  update();
}

template <int D>
void ParticleSystemModel<D>::draw(const aergia::CameraInterface *camera,
                                  ponos::Transform t) {
  if (!ps)
    return;
  _particlesSet->draw(camera, t);
  if (showVelocities)
    _velocitiesSet->draw(camera, t);
}

template <int D>
bool ParticleSystemModel<D>::intersect(const ponos::Ray3 &r, float *t) {
  if (!ps)
    return false;
  UNUSED_VARIABLE(t);
  furoo::Point<double, D> c;
  c.setX(r.o.x);
  c.setY(r.o.y);
  if (D == 3)
    c.setZ(r.o.z);
  _activeParticle = -1;
  ps->iterateParticles(furoo::BBox<double, D>(c, 5 * particleRadius),
                       [&](size_t id, furoo::Point<double, D> p) {
                         if (p.distance(c) <= particleRadius * 2)
                           _activeParticle = static_cast<int>(id);
                       });
  return _activeParticle >= 0;
}

template <int D> void ParticleSystemModel<D>::update() {
  updateParticles();
  updateVelocities();
}

template <int D> void ParticleSystemModel<D>::updateVelocities() {
  if (!ps)
    return;
  _velocitiesSet->resize(ps->size());
  double minValue = 0., maxValue = 0.;
  // unactivate dead particles, unfortunatelly we have to kill them all first
  // and make active ones live again
  {
    auto a = _velocitiesSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _velocitiesSet->count());
  }
  // update positions
  ps->iterateParticles([&](size_t id, furoo::Point<double, D> pos) {
    auto m = _velocitiesSet->transformBuffer(id);
    furoo::Vector<double, D> size;
    if (_curPropertyGroup >= 0)
      for (int d = 0; d < D; d++)
        size[d] = velocityScale *
                  ps->getScalarProperty(vectorGroups[_curPropertyGroup][d], id);
    float t[16];
    (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 1)) *
     ponos::scale(size[0], size[1], (D == 3) ? size[3] : 0.f))
        .matrix()
        .column_major(t);
    for (size_t k = 0; k < 16; k++)
      m[k] = t[k];
    aergia::Color color = aergia::COLOR_BLACK;
    auto c = _velocitiesSet->colorBuffer(id);
    c[0] = color.r;
    c[1] = color.g;
    c[2] = color.b;
    c[3] = color.a;
    auto a = _velocitiesSet->activeBuffer(id);
    a[0] = 1;
  });
}

template <int D> void ParticleSystemModel<D>::updateParticles() {
  if (!ps)
    return;
  _particlesSet->resize(ps->size());
  double minValue = 0., maxValue = 0.;
  bool usingPalette = false;
  _paletteD[_curProperty] = aergia::HEAT_MATLAB_PALETTE;
  if (_curProperty >= 0 &&
      _paletteD.find(static_cast<size_t>(_curProperty)) != _paletteD.end()) {
    minValue = ps->minPropertyValue(_curProperty);
    maxValue = ps->maxPropertyValue(_curProperty);
    usingPalette = true;
  }
  // unactivate dead particles, unfortunatelly we have to kill them all first
  // and make active ones live again
  {
    auto a = _particlesSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _particlesSet->count());
  }
  // update positions
  ps->iterateParticles([&](size_t id, furoo::Point<double, D> pos) {
    auto m = _particlesSet->transformBuffer(id);
    float t[16];
    (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
     ponos::scale(particleRadius, particleRadius, particleRadius))
        .matrix()
        .column_major(t);
    for (size_t k = 0; k < 16; k++)
      m[k] = t[k];
    aergia::Color color = particleColor;
    if (usingPalette) {
      auto t = static_cast<float>(ponos::smoothStep(
          ps->getScalarProperty(_curProperty, id), minValue, maxValue));
      color = _paletteD[_curProperty](t, particleColor.a);
    }
    auto c = _particlesSet->colorBuffer(id);
    c[0] = color.r;
    c[1] = color.g;
    c[2] = color.b;
    c[3] = color.a;
    auto a = _particlesSet->activeBuffer(id);
    a[0] = 1;
  });
}
