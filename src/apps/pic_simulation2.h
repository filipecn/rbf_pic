#ifndef FUROO_PIC_SIMULATION_H
#define FUROO_PIC_SIMULATION_H

#include <physics/simulation.h>

namespace furoo {

class PicSimulation2 : public Simulation2 {
public:
  PicSimulation2() = default;
  /// \param structure
  /// \param particles
  /// \param integrator
  /// \param solver
  PicSimulation2(StructureInterface<2> *structure, ParticleSystem2d *particles,
                 TemporalIntegratorInterface<2> *integrator, Solver<2> *solver);
  /// Define material type for each cell
  void markCells();

  /// Define material/boundary type for each face
  void markFaces();

  /// Define a boundary condition value for some wall side.
  /// \param side side of the domain to receive the boundary condition
  /// \param vel velocity vector to be assigned to boundary side
  virtual void defineBoundaryConditions(Definitions::Side side, Vector2d vel);

  /// Apply external forces (ex: gravity) to face velocities
  /// \param dt time step (seconds)
  virtual void computeExternalForces(double dt);

  /// Make all neumann face velocities zero
  /// Also apply any Dirichlet velocity that was previously defined through
  /// defineBoundaryContidions(.) function
  virtual void applyBoundaryCondition();

  /// Interpolate particle velocities into grid faces
  virtual void transferFromParticlesToGrid();

  /// Interpolate grid face velocities into particles
  virtual void transferFromGridToParticles();

  /// Computes pressure and correct velocities with pressure gradient
  virtual void solvePressure(double dt);

  /// Update particles positions
  virtual void advectParticles(double dt);

  /// Computes the biggest time step size that respects the cfl condition
  virtual double cfl() const;

protected:
  size_t surfaceMaskFieldId;
  size_t cellMaterialTypeFieldId;
  size_t faceMaterialTypeFieldId;
  size_t faceBoundaryTypeFieldId;
  size_t faceVelocityFieldId;
  size_t _faceBoundaryVelocity[2];
  size_t _facePressureFieldId;
  size_t particleVelocityFieldIds[2];
  size_t pressureFieldId;
  size_t gradientFieldId;
  size_t divergenceFieldId;
  Vector2d _gravity;
};

} // namespace furoo

#endif // FUROO_PIC_SIMULATION_H
