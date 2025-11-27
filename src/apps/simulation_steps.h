#ifndef FUROO_SIMULATION_STEPS_H
#define FUROO_SIMULATION_STEPS_H

#include <blas/rbf.h>
#include <blas/temporal_integrator.h>
#include <common/distribution.h>
#include <common/random.h>
#include <geometry/bbox.h>
#include <physics/injector.h>
#include <physics/particle_system.h>
#include <physics/solid_box.h>
#include <physics/solid_interface.h>
#include <structures/cell_graph.h>
#ifdef _USE_OPENMP
#include <mutex>
#endif // _USE_OPENMP


namespace furoo {

/// Holds a collection of implemented steps that can be used to setup simulation
/// methods steps. It helps quick tests on various combinations of different
/// types of steps. For example, one could easily get a pic algorithm and plug
/// different methods for particle advection, particle grid transferences, and
/// so on.
class SimulationSteps {
public:
  template <int D>
  static void initImplicitScene(
      const std::vector<int> &fluidBoxes,
      const std::vector<std::shared_ptr<SolidInterface<D>>> &solids,
      size_t maxLevel, size_t minLevel, size_t cellDistanceFieldId,
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialFieldId, size_t cellSurfaceMaskFieldId, bool graded);

  template <int D>
  static void
  applyFunctionToParticles(ParticleSystem<double, D> *particles,
                           size_t testFunctionPropertyId,
                           const std::function<double(Point<double, D>)> &f);

  ///
  /// \tparam D
  /// \param structure
  /// \param testFunctionFieldId
  /// \param cellParticleCenterFieldId
  /// \param f
  template <int D>
  static void
  applyFunctionToCells(StructureInterface<D> *structure,
                       size_t testFunctionFieldId,
                       size_t cellParticleCenterFieldId,
                       const std::function<double(Point<double, D>)> &f);

  ///
  ///
  /// \tparam D
  /// \param dimension
  /// \param cellGraph
  /// \param functionFieldId
  /// \param gradientFieldId
  /// \param cellMaterialFieldId
  /// \param cellParticleCenterFieldId
  /// \param rbf
  template <int D>
  static void computeGradient(int dimension, StructureInterface<D> *cellGraph,
                              size_t functionFieldId, size_t gradientFieldId,
                              size_t cellMaterialFieldId,
                              size_t cellParticleCenterFieldId,
                              DifferentialRBF<D> *rbf);

  ///
  /// \tparam D
  /// \param cellGraph
  /// \param functionFieldId
  /// \param laplcianFieldId
  /// \param cellMaterialFieldId
  /// \param cellParticleCenterFieldId
  /// \param rbf
  template <int D>
  static void
  computeLaplacian(StructureInterface<D> *cellGraph, size_t functionFieldId,
                   size_t laplcianFieldId, size_t cellMaterialFieldId,
                   size_t cellParticleCenterFieldId, DifferentialRBF<D> *rbf);

  template <int D>
  static void computeError(StructureInterface<D> *structure, size_t resultId,
                           size_t solutionId);

  /// Computes the cfl condition value based on particles velocities
  /// \param particleSystem raw pointer to particles system data
  /// \param fieldIds velocity component field ids
  /// \param dx smallest spatial length a particle is allowed to cross
  /// \return dx divided by maximum particle velocity
  template <int D>
  static double cfl(const ParticleSystem<double, D> *particleSystem,
                    const std::vector<size_t> &fieldIds, double dx);

  /// Computes the cfl condition value based on grid velocities
  /// \param structure raw pointer to structure data
  /// \param xField velocity x component field id
  /// \param yField velocity y component field id
  /// \param dx smallest spatial length a particle is allowed to cross
  /// \return dx divided by maximum particle velocity
  template <int D>
  static double cfl(const StructureInterface<D> *structure,
                    const std::vector<size_t> &fieldIds, double dx);

  /// Updates positions of particles though integrating over a given time step
  /// in Lagrangian manner. Domain is used to query particles domain in order
  /// to avoid particles leaving simulation domain or penetrating in solid
  /// regions \param particleSystem raw pointer to particles system data
  /// \param velocityProperties velocity component field ids
  /// \param domain raw pointer to structure data
  /// \param solids solid list
  /// \param integrator temporal integrator
  /// \param dt time step
  template <int D>
  static void
  advectParticles(ParticleSystem<double, D> *particleSystem,
                  const std::vector<size_t>& velocityProperties,
                  const StructureInterface<D> *domain,
                  const std::vector<std::shared_ptr<SolidInterface<D>>> &solids,
                  const TemporalIntegratorInterface<D> *integrator, double dt);

  /// INPUT: a structure with cells already marked with material types
  ///        all required fields must exist
  /// OUTPUT: material and boundary fields of faces marked, each face
  ///         recieves its material type and boundary type
  /// \param structure raw pointer to structure data
  /// \param cellMaterialFieldId field that stores cell material types
  /// \param faceMaterialFieldId field that stores face material types
  /// \param faceBoundaryFieldId field that stores face boundary types
  /// \param domainBoundaries boundary types of domain extremes
  template <int D>
  static void markFaces(StructureInterface<D> *structure,
                        size_t cellMaterialFieldId, size_t faceMaterialFieldId,
                        size_t faceBoundaryFieldId,
                        Definitions::Boundary *domainBoundaries);

  /// Mark as FLUID cells all cells that contain particles
  /// \param structure raw pointer to structure data
  /// \param cellMaterialFieldId field that stores cell material types
  /// \param particles raw pointer to particles system data
  template <int D>
  static void markCells(StructureInterface<D> *structure,
                        size_t cellMaterialFieldId,
                        ParticleSystem<double, D> *particles);

  /// Uses the sampleScalarPropery method of the particle system to
  /// interpolate the particle property to the specified grid field \param
  /// particles raw pointer to particles system data \param propertyId source
  /// scalar property field \param n number of source particles for property
  /// interpolation \param grid raw pointer to grid data \param fieldId
  /// destination scalar field \param useElement returns true if point should
  /// be a destination (to avoid unecessary interpolations)
  template <int D>
  static void
  transferFromParticlesToGrid(const ParticleSystem<double, D> *particles,
                              size_t propertyId, size_t n,
                              StructureInterface<D> *grid, size_t fieldId,
                              const std::function<bool(size_t)> &useElement);

  /// The particles used on interpolation are given by region queries that
  /// depend on the field location
  /// - face center location: particles that lie inside the two adjacent cells
  /// are used
  /// - vertex location: particles that lie inside all the connecting cells
  /// are used
  /// - cell location: particles that lie inside the cell are used
  /// The interpolation is done through the interpolator object
  /// \param particles raw pointer to particles system data
  /// \param propertyId source scalar property field
  /// \param grid raw pointer to grid data
  /// \param fieldId destination scalar field
  /// \param useElement returns true if point should be a destination (to
  /// avoid unecessary interpolations)
  template <int D>
  static void
  transferFromParticlesToGrid(const ParticleSystem<double, D> *particles,
                              size_t propertyId, StructureInterface<D> *grid,
                              size_t fieldId,
                              const std::function<bool(size_t)> &useElement,
                              Interpolant<Point<double, D>> *interpolator);
  template <int D>
  static void transferFromParticlesToGridHandlingCellLevels(
      const ParticleSystem<double, D> *particles, size_t propertyId,
      StructureInterface<D> *grid, size_t cellMaterialFieldId, size_t fieldId,
      const std::function<bool(size_t)> &useElement,
      Interpolant<Point<double, D>> *interpolator);

  template <int D>
  static void transferVelocitiesFromParticlesToMidPointFaces(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      const std::vector<size_t>& particleVelocityPropertyIds,
      const std::vector<size_t>& faceVelocityFieldIds, size_t faceMaterialFieldId,
      size_t cellMaterialFieldId, size_t faceMidPointFieldId);

  /// Interpolates a field into particles property through sample*Field
  /// methods of structures
  /// \param structure raw pointer to grid data
  /// \param fieldId source grid field id
  /// \param particles raw pointer to particles  system data
  /// \param propertyId destination property id
  /// \param useFace returns true if face should be a source
  /// \param useCell returns true if cell should be a source
  /// \param boundaryVelocity how boundary velocities are handled
  template <int D>
  static void
  transferFromGridToParticles(const StructureInterface<D> *grid, size_t fieldId,
                              ParticleSystem<double, D> *particles,
                              size_t propertyId,
                              const std::function<bool(size_t)> &useFace,
                              const std::function<bool(size_t)> &useCell,
                              Definitions::BoundaryVelocity boundaryVelocity);

  template <int D>
  static void transferFromGridToParticles(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialTypeFieldId, size_t faceMaterialTypeFieldId,
      size_t surfaceCellMaskFieldId, const std::vector<size_t>& faceVelocityFieldIds,
      const std::vector<size_t>& particleVelocityFieldIds,
      Interpolant<Point<double, D>> *interpolator);

  /// Applies acceleration to a scalar field (ex: gravity)
  /// f = f + acceleration * dt
  /// \param structure raw pointer to structure data
  /// \param fieldId scalar field index
  /// \param acceleration acceleration value
  /// \param dt time step
  /// \param useElement returns true if element should processed
  template <int D>
  static void
  applyExternalForce(StructureInterface<D> *structure, size_t fieldId,
                     double acceleration, double dt,
                     const std::function<bool(size_t)> &useElement = nullptr);

  template <int D>
  static void applyExternalForceToMidPointFaces(
      StructureInterface<D> *structure,
      const std::vector<size_t>& faceVelocityFieldIds,
      Vector<double, D> acceleration, double dt, size_t faceMaterialFieldId,
      size_t cellMaterialFieldId);

  ///
  /// \tparam D
  /// \param particles
  /// \param propertyId
  /// \param acceleration
  /// \param dt
  template <int D>
  static void applyExternalForce(ParticleSystem<double, D> *particles,
                                 size_t propertyId, double acceleration,
                                 double dt);

  template <int D>
  static void
  computeDivergenceOnParticles(ParticleSystem<double, D> *particles,
                               const std::vector<size_t>& velocityPropertyIds,
                               size_t divergencePropertyId,
                               DifferentialRBF<D> *rbf);
  /// Computes the divergence field in cell centers from a staggered grid via
  /// RBF gradient operator
  /// \tparam D dimension
  /// \param structure raw pointer to structure data
  /// \param divergenceCellFieldId cell centered divergence field id
  /// \param velocityFaceFieldIds face centered velocity component field ids
  /// (x,y,..)
  /// \param useCell filter cells where the divergence must be
  /// calculated (true if yes)
  /// \param rbf differential rbf
  /// \param adaptiveRBFRadius adapt radius depending on rbf stencil
  template <int D>
  static void
  computeDivergenceField(StructureInterface<D> *structure,
                         size_t divergenceCellFieldId,
                         const std::vector<size_t>& velocityFaceFieldIds,
                         const std::function<bool(size_t)> &useCell,
                         DifferentialRBF<D> *rbf, bool adaptiveRBFRadius);

  template <int D>
  static void computeDivergenceFieldFromMidPointFaces(
      StructureInterface<D> *structure, size_t divergenceCellFieldId,
      size_t faceMidPointField, const std::vector<size_t>& velocityFaceFieldIds,
      const std::function<bool(size_t)> &useCell, DifferentialRBF<D> *rbf);

  template <int D>
  static void computePressureGradientOnMidPointFaces(
      StructureInterface<D> *structure, size_t cellMaterialfield,
      size_t cellPressureField, size_t faceMaterialFieldId,
      const std::vector<size_t>& facePressureGradientFieldIds,
      DifferentialRBF<D> *rbf);

  /// Compute the cell centered divergence field in a staggered grid via RBF
  /// gradient operator also using cell centered velocities
  /// \param structure raw pointer to structure data
  /// \param divergenceCellFieldId cell centered divergence field id
  /// \param VelocityFaceFieldIds face centered velocity field ids
  /// \param VelocityCellFieldIds cell centered velocity field ids
  /// \param useCell filter cells where the divergence must be calculated
  /// (true if yes) \param rbf differential rbf \param adaptiveRBFRadius adapt
  /// radius depending on rbf stencil
  template <int D>
  static void
  computeDivergenceField(StructureInterface<D> *structure,
                         size_t divergenceCellFieldId,
                         const std::vector<size_t> &velocityFaceFieldIds,
                         const std::vector<size_t> &velocityCellFieldIds,
                         const std::function<bool(size_t)> &useCell,
                         DifferentialRBF<D> *rbf, bool adaptiveRBFRadius);

  /// Compute the cell centered divergence field in a staggered grid via RBF
  /// using particles as stencil points
  /// \tparam D dimension
  /// \param structure raw pointer to structure data
  /// \param particles
  /// \param divergenceCellFieldId
  /// \param cellMaterialFieldId
  /// \param cellParticleCenterFieldId cell's particle id
  /// \param velocityParticleFieldIds
  /// \param useCell
  /// \param rbf differential rbf
  template <int D>
  static void computeDivergenceFieldFromParticles(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t divergenceCellFieldId, size_t cellMaterialFieldId,
      size_t cellParticleCenterFieldId,
      const std::vector<size_t> &velocityParticleFieldIds,
      const std::function<bool(size_t)> &useCell, DifferentialRBF<D> *rbf);

  /// Apply falloff function to particles near solids. This method creates a
  /// continuity near the walls, where normal velocity should be zero
  /// \tparam D dimension
  /// \param structure raw pointer to structure data
  /// \param particles pointer to particles data
  /// \param cellMaterialFieldId
  /// \param velocityParticleFieldIds
  template <int D>
  static void
  applyFalloffToParticles(StructureInterface<D> *structure,
                          ParticleSystem<double, D> *particles,
                          size_t cellMaterialFieldId,
                          const std::vector<size_t> &velocityParticleFieldIds);

  /// Given an cell, setups the topological neighborhood and computes the
  /// weights of this stencil based on a rbf
  /// \param cellGraph raw pointer to cell graph data
  /// \param particles
  /// \param cellId cell index
  /// \param cellMaterialFieldId  cell material field id
  /// \param cellParticleCenterFieldId cell's particle id
  /// \param rbf differential rbf
  /// \param adaptiveRBFRadius adapt radius depending on rbf stencil
  /// \return the list of weights of the stencil where the first element is
  /// the center weight
  template <int D>
  static std::vector<double> computePressureMatrixStencilWeights(
      const StructureInterface<D> *structure, size_t cellId,
      size_t cellMaterialFieldId, size_t cellParticleCenterFieldId,
      DifferentialRBF<D> *rbf, bool adaptiveRBFRadius);

  /// Computes the pressure field of a cell graph based structure
  /// INPUT: - the cell centered divergence field must be previously computed
  /// \param structure raw pointer to structure data
  /// \param pressureCellFieldId cell centered pressure field id
  /// \param divergenceCellFieldId cell centered divergence field id
  /// \param cellMaterialFieldId  cell material field id
  /// \param dt time step
  /// \param rbf differential rbf
  /// \param adaptiveRBFRadius adapt radius depending on rbf stencil
  template <int D>
  static void
  solvePressure(StructureInterface<D> *structure,
                const ParticleSystem<double, D> *particles,
                size_t pressureCellFieldId, size_t divergenceCellFieldId,
                size_t cellMaterialFieldId, size_t surfaceMaskFieldId,
                size_t cellParticleCenterFieldId, double dt,
                DifferentialRBF<D> *rbf, bool adaptiveRBFRadius);

  /// Given an cell, setups the topological neighborhood and computes the
  /// weights of this stencil based on a rbf
  /// \param cellGraph raw pointer to cell graph data
  /// \param cellId cell index
  /// \param cellMaterialFieldId  cell material field id
  /// \param rbf differential rbf
  /// \param adaptiveRBFRadius adapt radius depending on rbf stencil
  /// \return the list of weights of the stencil where the first element is
  /// the center weight
  template <int D>
  static std::vector<double>
  computePressureMatrixStencilWeights(const StructureInterface<D> *structure,
                                      size_t cellId, size_t cellMaterialFieldId,
                                      DifferentialRBF<D> *rbf,
                                      bool adaptiveRBFRadius);

  template <int D>
  static std::vector<double> computePressureRingMatrixStencilWeights(
      const StructureInterface<D> *structure, size_t cellId,
      DifferentialRBF<D> *rbf, const std::vector<size_t> &neighbors);

  /// Computes the pressure field of a cell graph based structure
  /// INPUT: - the cell centered divergence field must be previously computed
  /// \param structure raw pointer to structure data
  /// \param pressureCellFieldId cell centered pressure field id
  /// \param divergenceCellFieldId cell centered divergence field id
  /// \param cellMaterialFieldId  cell material field id
  /// \param dt time step
  /// \param rbf differential rbf
  /// \param adaptiveRBFRadius adapt radius depending on rbf stencil
  template <int D>
  static void
  solvePressure(StructureInterface<D> *structure, size_t pressureCellFieldId,
                size_t divergenceCellFieldId, size_t cellMaterialFieldId,
                double dt, DifferentialRBF<D> *rbf, bool adaptiveRBFRadius);

  template <int D>
  static void solvePressureRing(StructureInterface<D> *structure,
                                size_t pressureCellFieldId,
                                size_t divergenceCellFieldId,
                                size_t cellMaterialFieldId, double dt,
                                DifferentialRBF<D> *rbf);

  /// Computes pressure gradient at face centers. Only FLUID faces are
  /// considered INPUT: Pressure field must be previously computed on cell and
  /// face centers
  /// \param cellGraph raw pointer to cell graph data
  /// \param pressureCellFieldId cell centered pressure field id
  /// \param pressureFaceFieldId face centered pressure field id
  /// \param pressureGradientFaceFieldIds face centered pressure gradient
  /// field ids \param faceMaterialFieldId face centered material field id
  template <int D>
  static void
  computePressureGradient(StructureInterface<D> *cellGraph,
                          size_t pressureCellFieldId,
                          size_t pressureFaceFieldId,
                          const std::vector<size_t>& pressureGradientFaceFieldIds,
                          size_t faceMaterialFieldId, DifferentialRBF<D> *rbf);
  template <int D>
  static void computePressureGradient(
      StructureInterface<D> *cellGraph, size_t faceMaterialFieldId,
      size_t pressureCellFieldId, size_t cellMaterialFIeldId,
      const std::vector<size_t>& pressureGradientFaceFieldIds,
      DifferentialRBF<D> *rbf);

  /// Computes pressure gradient at cells centers. Only FLUID cells are
  /// considered INPUT: Pressure field must be previously computed on cell
  /// centers
  /// \tparam D dimensions
  /// \param cellGraph raw pointer to cell graph data
  /// \param pressureCellFieldId cell centered pressure field id
  /// \param pressureGradientCellFieldIds face centered pressure gradient
  /// field ids \param cellMaterialFieldId cell centered material field id
  /// \param rbf differential rbf
  template <int D>
  static void
  computePressureGradient(StructureInterface<D> *cellGraph,
                          size_t pressureCellFieldId,
                          const std::vector<size_t>& pressureGradientCellFieldIds,
                          size_t cellMaterialFieldId, DifferentialRBF<D> *rbf);

  template <int D>
  static void computePressureGradient(
      StructureInterface<D> *cellGraph,
      const ParticleSystem<double, D> *particles, size_t pressureCellFieldId,
      const std::vector<size_t>& pressureGradientCellFieldIds,
      size_t cellMaterialFieldId, size_t cellParticleCenterFieldId,
      DifferentialRBF<D> *rbf);

  template <int D>
  static void computePressureGradient(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialFieldId, size_t pressureCellFieldId,
      size_t cellParticleCenterFieldId,
      const std::vector<size_t>& particlePressureGradientFieldIds,
      DifferentialRBF<D> *rbf);

  /// Correct velocities at faces through pressure gradient at faces
  /// performs velocity -= dt * pressureGradient
  /// \param structure raw pointer to structure data
  /// \param pressureGradientFaceFieldIds pressure gradient face field ids
  /// \param velocityFaceFieldIds face centered velocity field ids
  /// \param faceMaterialFieldId face material field id
  /// \param dt time step
  template <int D>
  static void
  correctVelocities(StructureInterface<D> *structure,
                    const std::vector<size_t>& pressureGradientFaceFieldIds,
                    const std::vector<size_t>& velocityFaceFieldIds,
                    size_t faceMaterialFieldId, double dt);

  template <int D>
  static void correctVelocitiesOnMidPointFaces(
      StructureInterface<D> *structure,
      const std::vector<size_t>& pressureGradientFaceFieldIds,
      const std::vector<size_t>& velocityFaceFieldIds, size_t faceMaterialFieldId,
      double dt);

  ///
  /// Use mid point face flow to correct velocities in the particles
  /// \tparam D dimension
  /// \param structure
  /// \param particles
  /// \param faceFlowFieldId id of fluid flow values
  /// \param particleVelocityIds
  template <int D>
  static void transferMidPointFacesFlowToCellVelocities(
      StructureInterface<D> *structure, size_t cellPressureFieldId,
      size_t cellMaterialFieldId, size_t faceMaterialFieldId,
      size_t faceMidPointFieldId, const std::vector<size_t> &cellVelocityIds,
      const std::vector<size_t> &faceVelocityIds, DifferentialRBF<D> *rbf,
      double dt);

  /// Correct velocities at particles by interpolating pressure gradient from
  /// cells into particles position first
  /// \tparam D dimensions
  /// \param structure raw pointer to structure data
  /// \param pressureGradientCellFieldId pressure gradient cell field id
  /// \param particles particles system raw pointer
  /// \param propertyId particle's velocity component property index
  /// \param useElement returns true if the cell holds a value for pressure
  /// gradient
  /// \param dt time step
  template <int D>
  static void
  correctVelocities(StructureInterface<D> *structure,
                    size_t pressureGradientCellFieldId,
                    ParticleSystem<double, D> *particles, size_t propertyId,
                    const std::function<bool(size_t)> &useElement, double dt);

  /// \param surfaceMaskCellFieldId surface mask cell field id
  template <int D>
  static void correctVelocities(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      const std::vector<size_t>& pressureGradientCellFieldIds,
      size_t surfaceMaskCellFieldId, size_t cellMaterialFieldId,
      size_t cellParticleCenterFieldId, const std::vector<size_t>& particleVelocityIds,
      double dt, Interpolant<Point<double, D>> *interpolator);

  template <int D>
  static void
  correctVelocities(ParticleSystem<double, D> *particles,
                    const std::vector<size_t>& particlePressureGradientFieldIds,
                    const std::vector<size_t>& particleVelocityFieldIds, double dt);

  template <int D>
  static void correctVelocitiesOnSurface(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      const std::vector<size_t>& pressureGradientCellFieldIds,
      size_t surfaceMaskCellFieldId, size_t cellMaterialFieldId,
      size_t cellParticleCenterFieldId, const std::vector<size_t>& particleVelocityIds,
      double dt, Interpolant<Point<double, D>> *interpolator);

  /// Propagate faces velocities from surface cells to air cells near the
  /// surface
  /// \tparam D dimension
  /// \param structure raw pointer to structure data
  /// \param velocityFieldIds list of velocity fields
  /// \param cellMaterialFieldId cell material field id
  /// \param faceMaterialFieldId face material field id
  /// \param surfaceMaskCellFieldId surface mask cell field id
  template <int D>
  static void propagateVelocities(StructureInterface<D> *structure,
                                  const std::vector<size_t>& velocityFieldIds,
                                  size_t cellMaterialFieldId,
                                  size_t faceMaterialFieldId,
                                  size_t surfaceMaskCellFieldId);

  /// Reseeds particles inside fluid cells leaving surface particles unchanged
  /// NOTE: This method does not set the velocity to the new particles
  /// \param structure raw pointer to structure data
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param particles raw pointer to particles system data
  template <int D>
  static void reseedParticles(const StructureInterface<D> *structure,
                              size_t cellMaterialFieldId,
                              size_t surfaceMaskCellFieldId,
                              ParticleSystem<double, D> *particles, size_t particleGroupPropertyId = D);

  /// Reseeds particles inside fluid cells leaving surface particles unchanged
  /// \param structure raw pointer to structure data
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param particles raw pointer to particles system data
  /// \param particleVelocityFields particle velocity property ids
  template <int D>
  static void
  reseedParticles(const StructureInterface<D> *structure,
                  size_t cellMaterialFieldId, size_t surfaceMaskCellFieldId,
                  ParticleSystem<double, D> *particles,
                  const std::vector<size_t> &particleVelocityFields);

  /// Overwrites the surface mask field
  /// \param structure raw pointer to structure data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellMaterialFieldId cell material field id
  /// \return set of fluid cells that belong to the surface
  template <int D>
  static std::set<size_t> markFluidSurface(StructureInterface<D> *structure,
                                           size_t surfaceMaskCellFieldId,
                                           size_t cellMaterialFieldId);

  /// Marks the surface, particles need to be used for new surface cells that
  /// might appear
  /// \param particles raw pointer to particles system data
  /// \param cellGraph raw pointer to cell graph data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceLevel surface tree level (usually the max graph level)
  static void refineFluidSurface(ParticleSystem2d *particles,
                                 CellGraph2 *cellGraph,
                                 size_t surfaceMaskCellFieldId,
                                 size_t cellMaterialFieldId,
                                 size_t surfaceLevel);

  /// Refines the graph towards the fluid body
  /// INPUT: a cellgraph with only the root node (unrefined)
  /// \param particles raw pointer to particles system data
  /// \param cellGraph raw pointer to cell graph data
  /// \param cellMaterialFieldId cell material field id
  /// \param fluidLevel fluid level
  template <int D>
  static void refineFluid(ParticleSystem<double, D> *particles,
                          StructureInterface<D> *cellGraph,
                          size_t cellMaterialFieldId, size_t fluidLevel);

  /// This method performs the following tasks:
  /// 1 - queue for every air cell that intersects any injector
  /// 2 - detect new fluid cells and refine them if necessary
  /// 3 - update surface
  /// \tparam D dimensions
  /// \param cellGraph raw pointer to cell graph data
  /// \param cellMaterialFieldId  cell materia field id
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param surfaceLevel surface tree level (usually the max graph level)
  /// \param particles raw pointer to particles system data
  /// \param injectors list of praticle injectors
  template <int D>
  static void
  handleInjectors(StructureInterface<D> *cellGraph, size_t cellMaterialFieldId,
                  size_t surfaceMaskCellFieldId, size_t surfaceLevel,
                  ParticleSystem<double, D> *particles,
                  const std::vector<ParticleInjector<D>> &injectors);

  /// Given a cell graph with cell material types computed and surface mask,
  /// coarse all cells following the criteria providade by doCoarse function
  /// INPUT: - cell materials and surface mask previously computed
  /// \param graph raw pointer to cell graph data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellMaterialFieldId cell material field id
  /// \param doCoarse callback to each cell, true if cell can be coarsed
  template <int D>
  static void coarseGraph(
      StructureInterface<D> *graph, size_t surfaceCellFieldId,
      size_t cellMaterialFieldId,
      const std::function<bool(Definitions::Material, size_t)> &doCoarse);

  /// After advection, particles have new positions and the surface is not the
  /// same anymore. This function updates surface cell mask and cell material
  /// fields.
  /// \param particles raw pointer to particle system data
  /// \param structure raw pointer to structure data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceLevel fluid surface tree level \return queue of surface
  /// cells ids
  template <int D>
  static std::queue<size_t>
  updateSurfaceCells(ParticleSystem<double, D> *particles,
                     StructureInterface<D> *structure,
                     size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
                     size_t surfaceLevel);

  template <int D>
  static std::queue<size_t> updateSurfaceCellsCoarseBottom(
      ParticleSystem<double, D> *particles, StructureInterface<D> *structure,
      size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
      size_t surfaceLevel);

  /// We must keep a narrow band that encloses the surface with cells all
  /// refined. This function refines to surface level both fluid and air cells
  /// that are too close to surface cells
  /// \param graph raw pointer to cell graph data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellTDistanceSurfaceFieldId field containing the topological
  /// distance to surface cells
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceCells queue of surface cells
  /// \return queue of surface cells (an updated version of the input param
  /// queue)
  template <int D>
  static std::queue<size_t> refineNearSurface(
      StructureInterface<D> *structure, size_t surfaceMaskCellFieldId,
      size_t cellTDistanceSurfaceFieldId, size_t cellMaterialFieldId,
      std::queue<size_t> &surfaceCells);

  ///
  /// In order to achieve better precision near solids, we need the cells there
  /// to have the same refinement as the surface. Otherwise, the suavization
  /// caused by zero normal velocity affect different cells with different
  /// sizes.
  /// \tparam D
  /// \param structure
  /// \param cellMaterialField
  /// \param surfaceLevel
  template <int D>
  static void refineNearSolid(StructureInterface<D> *structure,
                              size_t cellMaterialField, size_t surfaceLevel);

  /// Fluid cells that share a vertex with air cells also need to be marked as
  /// surface cells.
  /// \param structure raw pointer to structure data
  /// \param surfaceMaskCellFieldId surface mask cell field id
  /// \param cellMaterialFieldId cell material field id
  /// \param surfaceCells queue of surface cells
  template <int D>
  static void
  fixSurface(StructureInterface<D> *structure, size_t surfaceMaskCellFieldId,
             size_t cellMaterialFieldId, std::queue<size_t> &surfaceCells);

  /// Refine fluid cells to make the interior graded
  /// \param structure raw pointer to structure data
  /// \param cellMaterialFieldId cell material field id
  template <int D>
  static void makeGraded(StructureInterface<D> *structure,
                         size_t cellMaterialFieldId);

  /// Updates the graph given a new configuration of particles. This function
  /// calls:
  /// 1 - updateSurfaceCells
  /// 2 - coarseGraph
  /// 3 - refineNearSurface
  /// 4 - makeGraded
  /// 5 - fixSurface
  /// \param particles raw pointer to particles system data
  /// \param structure raw pointer to structure data
  /// \param surfaceMaskFieldId cell centered surface mask field id
  /// \param cellTDistanceSurfaceFieldId field containing the topological
  /// distance to surface cells
  /// \param cellMaterialFieldId cell centered material field id
  /// \param minLevel min tree level
  /// \param surfaceLevel max tree level
  /// \param graded calls step 4 to make the tree graded
  template <int D>
  static void
  updateTree(ParticleSystem<double, D> *particles,
             StructureInterface<D> *structure, size_t surfaceMaskFieldId,
             size_t surfaceTDistanceFieldId, size_t cellMaterialFieldId,
             size_t minLevel, size_t surfaceLevel, bool graded = true);

  template <int D>
  static void updateTreeCoarseBottom(ParticleSystem<double, D> *particles,
                                     StructureInterface<D> *structure,
                                     size_t surfaceMaskFieldId,
                                     size_t surfaceTDistanceFieldId,
                                     size_t cellMaterialFieldId,
                                     size_t minLevel, size_t surfaceLevel,
                                     bool graded = true);

  /// Cells that are intersected by solids are marked as solid cells
  /// \tparam D dimensions
  /// \param structure structure raw pointer
  /// \param cellMaterialFieldId cell centered material field id
  /// \param solids list of solid object pointers
  template <int D>
  static void
  markSolids(StructureInterface<D> *structure, size_t cellMaterialFieldId,
             const std::vector<std::shared_ptr<SolidInterface<D>>> &solids);

  /// Interpolates a cell centered field into face centered positions field.
  /// For each face, it uses the topological cell neighbors of its adjacent
  /// cells into as interpolation points.
  /// \param cellGraph raw pointer to cell graph data
  /// \param cellFieldId cell centered field
  /// \param faceFieldId face centered field
  /// \param useCell true if cell can be used
  /// \param useFace true if face can be used
  template <int D>
  static void
  transferCellFieldToFaceField(StructureInterface<D> *cellGraph,
                               size_t cellFieldId, size_t faceFieldId,
                               const std::function<bool(size_t)> &useCell,
                               const std::function<bool(size_t)> &useFace);

  /// Interpolates a face centered field into vertex centered field.
  /// \param cellGraph raw pointer to cell graph data
  /// \param faceFieldId face centered field
  /// \param vertexFieldId vertex centered field
  /// \param useCell true if cell must be used
  /// \param useFace true if face must be used
  static void
  transferFaceFieldToVertexField(CellGraph2 *cellGraph, size_t faceFieldId,
                                 size_t vertexFieldId,
                                 const std::function<bool(size_t)> &useCell,
                                 const std::function<bool(size_t)> &useFace);

  template <int D>
  static void transferCellVelocitiesToParticles(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialFieldId,
      const std::vector<size_t> &cellVelocityFieldIds,
      const std::vector<size_t> &particleVelocityFieldIds,
      Interpolant<Point<double, D>> *interpolant);

  ///
  /// Transfer face velocities at mid point position to particles. No-slip
  /// velocity condition is used at the boundary
  /// \tparam D dimension
  /// \param structure raw pointer to the structure
  /// \param particles raw pointer to the particles
  /// \param cellMaterialTypeFieldId cell material fiel id
  /// \param faceVelocityFieldIds
  /// \param particleVelocityFieldIds
  /// \param interpolator pointer to the interpolator that will be used to
  /// interpolater the velocities
  template <int D>
  static void transferMidPointFaceVelocitiesToParticles(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialTypeFieldId, const std::vector<size_t>& faceVelocityFieldIds,
      const std::vector<size_t>& particleVelocityFieldIds,
      Interpolant<Point<double, D>> *interpolator);

  ///
  /// Transfer face velocities at mid point position to particles. Slip
  /// condition for boundary velocities have to be specified as a parameter
  /// \tparam D dimension
  /// \param structure raw pointer to the structure
  /// \param particles raw pointer to the particles
  /// \param cellMaterialTypeFieldId cell material fiel id
  /// \param faceMaterialTypeFieldId face material fiel id
  /// \param faceVelocityFieldIds
  /// \param particleVelocityFieldIds
  /// \param boundaryVelocityCondition either noslip or free slip condition
  /// \param interpolator pointer to the interpolator that will be used to
  /// interpolater the velocities
  template <int D>
  static void transferMidPointFaceVelocitiesToParticles(
      StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
      size_t cellMaterialTypeFieldId, size_t faceMaterialTypeFieldId,
      const std::vector<size_t>& faceVelocityFieldIds,
      const std::vector<size_t>& particleVelocityFieldIds,
      Definitions::BoundaryVelocity boundaryVelocityCondition,
      Interpolant<Point<double, D>> *interpolator);

  /// For each fluid cell, computes the centroid of its particles
  /// \tparam D dimensions
  /// \param structure structure raw pointer
  /// \param particles particle system raw pointer
  /// \param cellMaterialFieldId cell centered material field id
  /// \param cellParticleCenterFieldId cell's particle id
  template <int D>
  static void computeCellParticlesCentroid(StructureInterface<D> *structure,
                                           ParticleSystem<double, D> *particles,
                                           size_t surfaceMaskCellFieldId,
                                           size_t cellMaterialFieldId,
                                           size_t cellParticleCenterFieldId);

  template <int D>
  static void computeFaceMidPoint(StructureInterface<D> *structure,
                                  size_t faceMidPointFieldId,
                                  size_t faceMaterialField);
};

#include "simulation_steps.inl"
#include "simulation_steps_trask.inl"

} // namespace furoo

#endif // FUROO_SIMULATION_STEPS_H
