#ifndef FUROO_PHYSICS_PARTICLE_SYSTEM_2_H
#define FUROO_PHYSICS_PARTICLE_SYSTEM_2_H

#include <blas/interpolation.h>
#include <blas/linear_vector.h>
#include <memory>
#include <structures/point_set.h>

namespace furoo {

class ParticleSystem2 {
public:
  ParticleSystem2();
  ParticleSystem2(PointSet2Interface *ps,
    Interpolant<Point2d> *i = defaultInterpolant());
  virtual ~ParticleSystem2();
  void setDomainRegion(const BBox2D &);
  double sampleScalarProperty(size_t propertyId, double r,
                              const Point2d &p, const Vector2d& reflectionAxis = Vector2d(0.,0.)) const;
  double sampleScalarProperty(size_t propertyId, const Point2d &p, size_t n) const;
  void setScalarProperty(size_t, size_t, double);
  double getScalarProperty(size_t, size_t) const;
  size_t addScalarProperty(double v = 0.);
  size_t addParticle(Point2d);
  Point2d operator[](size_t) const;
  void setPosition(unsigned int, Point2d);
  size_t size();
  void
  iterateParticles(const std::function<void(unsigned int, Point2d)> &) const;
  void
  iterateParticles(const BBox2D &,
                   const std::function<void(unsigned int, Point2d)> &) const;
  void iterateClosestParticles(Point2d, size_t, const std::function<void(size_t, Point2d)>&) const;

protected:
  std::vector<LinearVector> _scalarProperties;
  std::unique_ptr<PointSet2Interface> _pointSet;
  std::unique_ptr<Interpolant<Point2d>> _interpolant;

private:
  static Interpolant<Point2d>* defaultInterpolant()
  {
    static QuinticKernel<Point2d> kernel;
    return new RBFInterpolant<Point2d>{kernel};
  }
};

} // namespace furoo

#endif // FUROO_PHYSICS_PARTICLE_SET_2_H
