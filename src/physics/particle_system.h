#ifndef FUROO_PARTICLE_SYSTEM_H
#define FUROO_PARTICLE_SYSTEM_H

#include <blas/interpolation.h>
#include <common/code_profiler.h>
#include <memory>
#include <structures/point_set.h>

namespace furoo {

template<typename T, int D> class ParticleSystem {
public:
  ParticleSystem();

  /// \param pointSet point set object pointer
  /// \param interpolant (interpolation used for sampling particle properties)
  explicit ParticleSystem(PointSetInterface<T, D> *pointSet,
                          Interpolant<Point<T, D>> *interpolant = new MLS<D>());

  virtual ~ParticleSystem();

  /// \param region domain region where particles can exist
  void setDomainRegion(const BBox<T, D> &region);

  /// \param propertyId property field id
  /// \param radius search ball radius
  /// \param target search ball center, interpolation point
  /// \return interpolated property value at target
  double sampleScalarProperty(size_t propertyId, double radius,
                              const Point<T, D> &target) const;

  /// \param propertyIds properties corresponding to each vector component
  /// \param radius search ball radius
  /// \param target search ball center, interpolation point
  /// \return interpolated vector of property value of each component at target
  Vector<double, D>
  sampleVectorFromProperties(const std::vector<size_t> &propertyIds,
                             double radius, const Point<T, D> &target) const;

  /// \param propertyIds properties corresponding to each vector component
  /// \param target interpolation point
  /// \param n number of closest particles
  /// \return interpolated vector of property value of each component at target
  Vector<double, D>
  sampleVectorFromProperties(const std::vector<size_t> &propertyIds,
                             const Point<T, D> &target, size_t n) const;

  /// Interpolates a property field based on a set of closest points
  /// \param propertyId property field id
  /// \param p interpolation point
  /// \param n number of closest points
  /// \return interpolated property value at target
  double sampleScalarProperty(size_t propertyId, const Point<T, D> &p,
                              size_t n) const;

  /// \param propertyId property field id
  /// \param particleId particle id
  /// \param value new value
  void setScalarProperty(size_t propertyId, size_t particleId, double value);

  /// \param propertyId property field id
  /// \param particleId particle id
  /// \return scalar property value at particle id
  double getScalarProperty(size_t propertyId, size_t particleId) const;

  /// Builds a scalar vector from a set of D properties of the same type
  /// \param propertyIds list of properties
  /// \param particleId particle id
  /// \return the vector with the values of each property on its respective
  /// component
  Vector<T, D> getVector(const std::vector<size_t> &propertyIds,
                         size_t particleId) const;

  /// \param v initial value
  /// \return field id of new property field
  size_t addScalarProperty(double v = 0.);

  /// \param position particle's position
  /// \return particle id
  size_t addParticle(Point<T, D> position);

  /// \param id particle id
  /// \return particle's position
  Point<T, D> operator[](size_t id) const;

  /// \param propertyId property field id
  /// \return maximum value of the property field
  double maxPropertyValue(size_t propertyId) const;

  /// \param propertyId property field id
  /// \return minimum value of the property field
  double minPropertyValue(size_t propertyId) const;

  /// \param particleId particle id
  /// \param position new position
  void setPosition(size_t particleId, Point<T, D> position);

  /// \param i particle id to be removed
  void remove(size_t i);

  /// remove all particles (doesn't free memory)
  void clear();

  /// update underlying point set (for OpenMP, need to be fixed)
  void update();

  /// \return number of active particles
  size_t size();

  /// \param f callback for each particle
  void
  iterateParticles(const std::function<void(size_t, Point<T, D>)> &f) const;

#ifdef _USE_OPENMP
  /// Iterate over particles in parallel using OpenMP
  /// \param f callback (no race conditions) for each particle
  void iterateParticles_par(const std::function<void(size_t,
                                                     Point<T, D>)> &f) const;
#endif // _USE_OPENMP

  /// \param searchRegion bounding box shape search region
  /// \param f callback for each found particle
  void
  iterateParticles(const BBox<T, D> &searchRegion,
                   const std::function<void(size_t, Point<T, D>)> &f) const;

  void
  iterateParticlesZ(const BBox<T, D> &b, size_t level,
                    const std::function<void(size_t, Point<T, D>)> &f) const;

  /// Checks if there is at least one particle inside search region
  /// \param searchRegion search region
  /// \return true if particle exists
  bool containsParticle(const BBox<T, D> &searchRegion) const;
  bool containsParticleZ(const BBox<T, D> &searchRegion, size_t level) const;

  /// k-nn search
  /// \param center center point of search
  /// \param neighborCount number of nearest neighbors
  /// \param f callback to each nearest neighbor
  void iterateClosestParticles(
      Point<T, D> center, size_t neighborCount,
      const std::function<void(size_t, Point<T, D>)> &f) const;

  size_t scalarPropertiesCount() const { return _scalarProperties.size(); }

protected:
  std::vector<LinearVector> _scalarProperties;
  std::unique_ptr<PointSetInterface<T, D>> _pointSet;
  std::unique_ptr<Interpolant<Point<T, D>>> _interpolant;
};

#include "physics/particle_system.inl"

typedef ParticleSystem<double, 2> ParticleSystem2d;
typedef ParticleSystem<double, 3> ParticleSystem3d;

} // namespace furoo

#endif // FUROO_PARTICLE_SYSTEM_H
