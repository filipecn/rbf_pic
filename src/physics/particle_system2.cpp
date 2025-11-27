#include <physics/particle_system2.h>
#include <geometry/transform.h>

namespace furoo {

ParticleSystem2::ParticleSystem2() = default;

ParticleSystem2::ParticleSystem2(PointSet2Interface *ps,
                                 Interpolant<Point2d> *i) {
  _pointSet.reset(ps);
  _interpolant.reset(i);
}

ParticleSystem2::~ParticleSystem2() = default;

double ParticleSystem2::sampleScalarProperty(size_t propertyId, double r,
                                             const Point2d &target, const Vector2d& reflectionAxis) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  std::vector<Point2d> points;
  std::vector<double> values;
  _pointSet->search(BBox2D(target, r), [&](size_t i) {
    points.emplace_back((*_pointSet)[i]);
    values.emplace_back(_scalarProperties[propertyId][i]);
  });
  if(std::fabs(reflectionAxis.x()) > 0. || std::fabs(reflectionAxis.y()) > 0.) {
    // centroid
    Point2d centroid(0.,0.);
    for(auto& p : points)
      centroid += p;
    ASSERT_FATAL(points.size() > 0);
    centroid /= 1.0 * points.size();
    values.emplace_back(_interpolant->interpolateAt(centroid, points, values));
    auto r =
        Transform2D::applyReflection(centroid - target, reflectionAxis.normalized());
    points.emplace_back(target + r);
  }
  return _interpolant->interpolateAt(target, points, values);
}

double ParticleSystem2::sampleScalarProperty(size_t propertyId, const Point2d &p, size_t n) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  std::vector<Point2d> points;
  std::vector<double> values;
  _pointSet->iterateClosestPoints(p, n, [&](size_t i, Point2d pos){
    points.emplace_back(pos);
    values.emplace_back(_scalarProperties[propertyId][i]);
  }, nullptr);
  return _interpolant->interpolateAt(p, points, values);
}

size_t ParticleSystem2::addScalarProperty(double v) {
  _scalarProperties.emplace_back(_pointSet->size(), v);
  return _scalarProperties.size() - 1;
}

size_t ParticleSystem2::addParticle(Point2d position) {
  for (size_t i = 0; i < _scalarProperties.size(); i++)
    _scalarProperties[i].emplace_back(0.);
  return _pointSet->add(position);
}

void ParticleSystem2::setDomainRegion(const BBox2D &region) {
  _pointSet->setDomainRegion(region);
}

void ParticleSystem2::setScalarProperty(size_t propertyId, size_t particleId,
                                        double value) {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  ASSERT_FATAL(particleId < _pointSet->size());
  _scalarProperties[propertyId][particleId] = value;
}

double ParticleSystem2::getScalarProperty(size_t propertyId,
                                          size_t particleId) const {
  ASSERT_FATAL(propertyId < _scalarProperties.size());
  ASSERT_FATAL(particleId < _pointSet->size());
  return _scalarProperties[propertyId][particleId];
}

Point2d ParticleSystem2::operator[](size_t id) const {
  ASSERT_FATAL(id < _pointSet->size());
  return (*_pointSet)[id];
}

size_t ParticleSystem2::size() { return _pointSet->size(); }

void ParticleSystem2::setPosition(unsigned int i, Point2d p) {
  _pointSet->setPosition(i, p);
}

void ParticleSystem2::iterateParticles(
    const std::function<void(unsigned int, Point2d)> &f) const {
  _pointSet->iteratePoints(f);
}

void ParticleSystem2::iterateParticles(
    const BBox2D &b,
    const std::function<void(unsigned int, Point2d)> &f) const {
  _pointSet->search(b, [&](size_t i) { f(i, (*_pointSet)[i]); });
}

void ParticleSystem2::iterateClosestParticles(Point2d center,
                                              size_t n,
                                              const std::function<void(size_t, Point2d)> &f) const {
  _pointSet->iterateClosestPoints(center, n, f, nullptr);
}

} // namespace furoo
