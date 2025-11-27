#ifndef FUROO_SIMULATION_H
#define FUROO_SIMULATION_H

#include <blas/temporal_integrator.h>
#include <common/definitions.h>
#include <physics/particle_system.h>
#include <physics/solver.h>
#include <structures/structure_interface.h>
#include <memory>

namespace furoo {

/// Simulation Class
/// General class used for physics simulation. Must set a structure,
/// particle system, temporal integrator and solver.
/// Can be used for grid or hybrid simulation purposes.
/// \tparam D number of dimensions
template<int D>
class Simulation {
 public:
  Simulation() = default;
  /// \param structure data structure object pointer
  /// \param particles
  /// \param integrator
  /// \param solver
  explicit Simulation(StructureInterface<D> *structure,
                      ParticleSystem<double, D> *particles,
                      TemporalIntegratorInterface<D> *integrator = nullptr,
                      Solver<D> *solver = nullptr);
  ~Simulation();

  /// \param region define simulation domain region
  void setRegion(BBox<double, D> region);

  /// \param integrator (for advecting particles)
  void setIntegrator(TemporalIntegratorInterface<D> *integrator);

  /// \param solver (for pressure solve)
  void setSolver(Solver<D> *solver);

  /// \return particles system
  ParticleSystem<double, D> *particles();

  /// \return structure
  StructureInterface<D> *structure();

 protected:
  /// Respects the following order: LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
  Definitions::Boundary _domainBoundaries[D * 2];
  std::shared_ptr<StructureInterface<D>> _structure;
  std::shared_ptr<ParticleSystem<double, D>> _particles;
  std::shared_ptr<TemporalIntegratorInterface<D>> _integrator;
  std::shared_ptr<Solver<D>> _solver;
};

#include "physics/simulation.inl"

typedef Simulation<2> Simulation2;
typedef Simulation<3> Simulation3;
}

#endif //FUROO_SIMULATION_H
