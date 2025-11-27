#ifndef FUROO_BLAS_TEMPORAL_INTEGRATOR_H
#define FUROO_BLAS_TEMPORAL_INTEGRATOR_H

#include <geometry/point.h>

namespace furoo {

/// Responsable ofr integrating particles with velocity over time steps
/// and computing their final positions
/// \tparam D space dimension
template<size_t D>
class TemporalIntegratorInterface {
 public:
  TemporalIntegratorInterface() = default;
  virtual ~TemporalIntegratorInterface() = default;

  /// Computes the end point of a particle moving with constant velocity over
  /// a time step period of time
  /// \param position particle's starting position
  /// \param velocity particle's velocity
  /// \param timeStep time step
  /// \return the final position after integration
  virtual Point<double, D> integrate(Point<double, D> position,
                                     Vector<double, D> velocity,
                                     double timeStep) const = 0;
};

/// Implements straightforward euler method of integration doing
/// finalPosition = startPosition + timeStep * velocity
/// \tparam D space dimension
template<size_t D>
class EulerTemporalIntegrator : public TemporalIntegratorInterface<D> {
 public:
  EulerTemporalIntegrator() = default;
  ~EulerTemporalIntegrator() = default;

  Point<double, D> integrate(Point<double, D>, Vector<double, D>,
                             double) const override;
};

/// Implements the 2nd order Runge Kutta integration by breaking timeStep
/// in 3 pieces and then computing the finalPosition as
/// startPosition + (2. / 9.) * dt * k1 + (3. / 9.) * dt * k2 + (4. / 9.) * dt * k3
/// where ki is an intermediary velocity
/// \tparam D space dimention
template<size_t D>
class RungeKutta2TemporalIntegrator : public TemporalIntegratorInterface<D> {
 public:
  RungeKutta2TemporalIntegrator() = default;
  ~RungeKutta2TemporalIntegrator() = default;

  Point<double, D> integrate(Point<double, D>, Vector<double, D>,
                             double) const override;
};

#include <blas/temporal_integrator.inl>

typedef EulerTemporalIntegrator<2> EulerTemporalIntegrator2;

} // namespace furoo

#endif
