#ifndef FUROO_INJECTOR_H
#define FUROO_INJECTOR_H

#include <common/io.h>
#include <physics/particle_system.h>

namespace furoo {

template <int D> class ParticleInjector {
public:
  /// \param particleSystem particle system point object
  /// \param spacing spacing between particles
  explicit ParticleInjector(ParticleSystem<double, D> *particleSystem,
                            double spacing = 0.03);

  /// \param region injector seeding region
  /// \param transform transform applied to injector region
  /// \param velocity new particle velocity
  /// \param inflowRate inflow rate (particles per second)
  /// \param particleSystem particle system raw pointer
  explicit ParticleInjector(BBox<double, D> region, Transform<D> transform,
                            Vector<double, D> velocity, double inflowRate,
                            ParticleSystem<double, D> *particleSystem);

  void fromFile(std::string filename);

  /// \param s particle spacing
  void setSpacing(double s);

  /// Sets starting time and ending time of partcile injection
  /// \param startTime start time
  /// \param endTime end time
  void setStartEndTime(double startTime, double endTime = INFINITY);

  /// Fill box region with particles
  /// \param region
  /// \return number of particles injected
  size_t setupBoxShape(BBox<double, D> region, int group = 0, int groupPropertyId = D);

  /// Fill ball region with particles
  /// \param center
  /// \param radius
  /// \return number of particles injected
  size_t setupCircleShape(Point<double, D> center, double radius);

  /// Load particles from file
  /// \param filename file path
  void loadFromFile(const char *filename);

  /// Computes the number of particles to be injected over the time step and add
  /// them to the simulation.
  /// \param dt time step
  /// \return Return the bounding box surronding the injector. If the injector
  /// is rotated, return the bounding box of the final rotated injector.
  size_t inject(double dt);

  /// Injection region bounding box
  /// \return world space region
  BBox<double, D> region() const;

private:
  BBox<double, D> _worldSpaceRegion;
  BBox<double, D> _injectionRegion;
  Transform<D> _transform;
  Vector<double, D> _velocity;
  double _inflowRate = 500;
  double _spacing = 0.03;
  double _elapsedTime;
  double _startTime;
  double _endTime;
  ParticleSystem<double, D> *_particleSystem;
};

#include <physics/injector.inl>

typedef ParticleInjector<2> Injector2;
typedef ParticleInjector<3> Injector3;

} // namespace furoo

#endif // FUROO_INJECTOR_H
