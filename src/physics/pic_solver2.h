#ifndef FUROO_PHYSICS_PIC_SOLVER_2_H
#define FUROO_PHYSICS_PIC_SOLVER_2_H

#include <animation/physics_animation.h>
#include <blas/temporal_integrator.h>
#include <physics/particle_system2.h>
#include <physics/simulation_domain.h>

namespace furoo {

class PICSolver2 : public PhysicsAnimation {
public:
  PICSolver2();
  PICSolver2(DomainBoundaryHandler2 *, SimulationDomain2 *, ParticleSystem2 *,
             EulerTemporalIntegrator2 * = new EulerTemporalIntegrator2());
  virtual ~PICSolver2();

  void setBoundaryHandler(DomainBoundaryHandler2 *);
  void setSimulationDomain(SimulationDomain2 *);
  void setParticleSystem(ParticleSystem2 *);
  SimulationDomain2 *domain();
  ParticleSystem2 *particles();
  DomainBoundaryHandler2* boundaryHandler();
  void setGravity(Vector2d);

  LinearMatrix pressureMatrix;
  LinearVector pressureRhs;
  LinearVector pressureSolution;
  LinearVector pressureDivergent;

  void exportVelocityField(int frame) const;

protected:
  void onInitialize() override;
  void onAdvanceTimeStep(double) override;
  double cfl(double) const;
  unsigned int numberOfSubTimeSteps(double) const override;
  void applyBoundaryCondition();
  void beginAdvanceTimeStep(double);
  void computeExternalForces(double);
  void computeViscosity(double);
  void computePressure(double);
  void correctVelocities(double);
  void propagateVelocitiesToDomain();
  void computeAdvection(double);
  void endAdvanceTimeStep(double);
  void transferFromParticlesToGrid();
  void transferFromGridToParticles();

  unsigned int _velocityField;
  unsigned int _pressureField;
  unsigned int _particleVx;
  unsigned int _particleVy;
  DomainBoundaryHandler2* _boundaryHandler;
  std::shared_ptr<SimulationDomain2> _domain;
  ParticleSystem2* _particles;
  Vector2d _gravity = Vector2d(0, -9.8);
  double _maxCfl = 5.;
  // TODO always euler??
  std::unique_ptr<EulerTemporalIntegrator2> _temporalIntegrator;
};

} // namespace furoo

#endif // FUROO_PHYSICS_PIC_SOLVER_2_H
