#include "apps/pic_simulation2.h"
#include <blas/rbf.h>
#include <map>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

namespace furoo {

PicSimulation2::PicSimulation2(StructureInterface<2> *structure,
                               ParticleSystem2d *particles,
                               TemporalIntegratorInterface<2> *integrator,
                               Solver<2> *solver) : Simulation2(structure,
                                                                particles,
                                                                integrator,
                                                                solver) {
  // setup fields
  // FACE BOUNDARY TYPE 0
  faceBoundaryTypeFieldId = _structure->addFaceField<Definitions::Boundary>(
      Definitions::Boundary::NONE);
  _structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId)->increaseSize(
      _structure->faceCount());
  // FACE MATERIAL TYPE 1
  faceMaterialTypeFieldId =
      _structure->addFaceField<Definitions::Material>(
          Definitions::Material::AIR);
  _structure->field<Definitions::Material>(faceMaterialTypeFieldId)->increaseSize(
      _structure->faceCount());
  // FACE VELOCITY FIELD 2
  faceVelocityFieldId = _structure->addFaceField<double>(0.);
  _structure->field<double>(faceVelocityFieldId)->increaseSize(
      _structure->faceCount());
  // CELL MATERIAL TYPE 3
  cellMaterialTypeFieldId = _structure->addCellField<Definitions::Material>(
      Definitions::Material::AIR);
  _structure->field<Definitions::Material>(cellMaterialTypeFieldId)->increaseSize(
      _structure->cellCount());
  // CELL SURFACE MASK 4
  surfaceMaskFieldId = _structure->addCellField<unsigned char>(0);
  _structure->field<unsigned char>(surfaceMaskFieldId)->increaseSize(
      _structure->cellCount());
  // CELL PRESSURE FIELD 5
  pressureFieldId = _structure->addCellField<double>(0.);
  _structure->field<double>(pressureFieldId)->increaseSize(_structure->cellCount());
  // FACE GRADIENT FIELD 6
  gradientFieldId = _structure->addFaceField<double>(0.);
  _structure->field<double>(gradientFieldId)->increaseSize(_structure->faceCount());
  // CELL DIVERGENCE FIELD 7
  divergenceFieldId = _structure->addCellField<double>(0.);
  _structure->field<double>(divergenceFieldId)->increaseSize(_structure->cellCount());
  for (size_t &particleVelocityField : particleVelocityFieldIds)
    particleVelocityField = _particles->addScalarProperty();
  for (auto &_domainBoundarie : _domainBoundaries)
    _domainBoundarie = Definitions::Boundary::NEUMANN;
}

void PicSimulation2::markCells() {
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  _structure->iterateCells([&](size_t cellId) {
    size_t count = 0;
    _particles->iterateParticles(_structure->cellRegion(cellId),
                                 [&](size_t id, Point2d p) {
                                   UNUSED_VARIABLE(id);
                                   UNUSED_VARIABLE(p);
                                   count++;
                                 });
    cellMaterialField[cellId] =
        (count) ? Definitions::Material::FLUID : Definitions::Material::AIR;
  });
}

void PicSimulation2::markFaces() {
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  auto &faceBoundaryField =
      *_structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId);
  _structure->iterateFaces([&](size_t faceId) {
    Definitions::Boundary faceBoundaryType = Definitions::Boundary::NONE;
    auto cells = _structure->faceCells(faceId);
    if (cells.size() == 1) {
      cells.push_back({-1, Definitions::oppositeSide(cells[0].side)});
    }
    THROW(cells.size() == 2u, "wrong number of face cells.");
    for (size_t n = 0; n < 2; n++)
      if (cells[n].id < 0)
        switch (cells[n].side) {
          case Definitions::Side::LEFT: faceBoundaryType = _domainBoundaries[0];
            break;
          case Definitions::Side::BOTTOM:
            faceBoundaryType = _domainBoundaries[1];
            break;
          case Definitions::Side::RIGHT:faceBoundaryType = _domainBoundaries[2];
            break;
          case Definitions::Side::TOP: faceBoundaryType = _domainBoundaries[3];
            break;
          default: NOT_IMPLEMENTED();
        }
    std::vector<Definitions::Material>
        cellType(2, Definitions::Material::SOLID);
    for (size_t i = 0; i < 2; i++)
      if (cells[i].id >= 0)
        cellType[i] = cellMaterialField[cells[i].id];
    if (cells[0].id >= 0 && cells[1].id >= 0)
      if ((cellType[0] == Definitions::Material::FLUID &&
          cellType[1] == Definitions::Material::AIR) ||
          (cellType[1] == Definitions::Material::FLUID &&
              cellType[0] == Definitions::Material::AIR))
        faceBoundaryType = Definitions::Boundary::DIRICHLET;
    faceBoundaryField[faceId] = faceBoundaryType;
    auto faceMaterialType = Definitions::Material::AIR;
    if (cells[0].id < 0)
      faceMaterialType =
          (faceBoundaryType == Definitions::Boundary::DIRICHLET) ? cellType[1]
                                                                 : Definitions::Material::SOLID;
    else if (cells[1].id < 0)
      faceMaterialType =
          (faceBoundaryType == Definitions::Boundary::DIRICHLET) ? cellType[0]
                                                                 : Definitions::Material::SOLID;
    else {
      if (cellType[0] == Definitions::Material::FLUID
          || cellType[1] == Definitions::Material::FLUID)
        faceMaterialType = Definitions::Material::FLUID;
      else if (cellType[0] == Definitions::Material::SOLID
          || cellType[1] == Definitions::Material::SOLID)
        faceMaterialType = Definitions::Material::SOLID;
      else
        faceMaterialType = Definitions::Material::AIR;
    }
    faceMaterialField[faceId] = faceMaterialType;
  });
}

void PicSimulation2::computeExternalForces(double dt) {
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  _structure->iterateFaces([&](size_t id) {
    if (_structure->faceOrientation(id) == Definitions::Orientation::HORIZONTAL)
      velocityField[id] -= 9.81 * dt;
    else
      velocityField[id] += 0;
  });
}

void PicSimulation2::applyBoundaryCondition() {
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  auto &boundaryField =
      *_structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId);
  _structure->iterateFaces([&](size_t id) {
    if (boundaryField[id] == Definitions::Boundary::NEUMANN)
      velocityField[id] = 0.;
  });
}

void PicSimulation2::transferFromParticlesToGrid() {
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  _structure->iterateFaces([&](size_t id) {
    velocityField[id] = 0.;
    if (faceMaterialField[id] == Definitions::Material::FLUID) {
      size_t fieldId = particleVelocityFieldIds[0];
      if (_structure->faceOrientation(id)
          == Definitions::Orientation::HORIZONTAL)
        fieldId = particleVelocityFieldIds[1];
      velocityField[id] = _particles->sampleScalarProperty(fieldId,
                                                           _structure->faceCenterPosition(
                                                               id), 13);
    }
  });
}

void PicSimulation2::transferFromGridToParticles() {
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  std::function<bool(size_t)> isValid = [&](size_t id) -> bool {
    return faceMaterialField[id] == Definitions::Material::FLUID;
  };
  _particles->iterateParticles([&](size_t id, Point2d p) {
    _particles->setScalarProperty(particleVelocityFieldIds[0], id,
                                  _structure->sampleFaceField(
                                      faceVelocityFieldId,
                                      p,
                                      Definitions::Orientation::VERTICAL,
                                      isValid));
    _particles->setScalarProperty(particleVelocityFieldIds[1], id,
                                  _structure->sampleFaceField(
                                      faceVelocityFieldId,
                                      p,
                                      Definitions::Orientation::HORIZONTAL,
                                      isValid));
  });
}

void PicSimulation2::solvePressure(double dt) {
  DifferentialRBF<2> rbf;
  rbf.setKernel(new GaussianKernel<Point2d>(1.));
  size_t fluidCount = 0;
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  auto &pressureField = *_structure->field<double>(pressureFieldId);
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      fluidCount++;
  });
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(fluidCount + 1);
//  std::cerr << "fluid count " << fluidCount << std::endl;
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  std::vector<Eigen::Triplet<double>> triplets;
  // structure index -> matrix index
  std::map<size_t, size_t> indexMap;
  std::function<size_t(size_t)> cellIndex = [&](size_t cellId) {
    auto it = indexMap.find(cellId);
    if (it == indexMap.end()) {
      size_t nextIndex = indexMap.size();
      indexMap[cellId] = nextIndex;
      return indexMap.size() - 1;
    }
    return it->second;
  };
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    auto cellRegionSize = _structure->cellRegion(cellId).size();
    auto center = _structure->cellCenterPosition(cellId);
    auto neighbors = _structure->cellNeighbors(cellId);
    auto faces = _structure->cellFaces(cellId);
    { // DIV
      double div = 0.;
      // compute divergence
      { // X direction
        std::vector<Point2d> points;
        std::vector<double> values;
        for (auto f : faces)
          if (_structure->faceOrientation(f)
              == Definitions::Orientation::HORIZONTAL) {
            values.emplace_back(velocityField[f]);
            points.emplace_back(_structure->faceCenterPosition(f));
          }
        div += rbf.gradientAt(center, 1, points, values);
//        std::cerr << "div at y" << div << std::endl;
      }
      { // Y direction
        std::vector<Point2d> points;
        std::vector<double> values;
        for (auto f : faces)
          if (_structure->faceOrientation(f)
              == Definitions::Orientation::VERTICAL) {
            values.emplace_back(velocityField[f]);
            points.emplace_back(_structure->faceCenterPosition(f));
          }
        div += rbf.gradientAt(center, 0, points, values);
//        std::cerr << "div at x" <<
//                  rbf.gradientAt(center, 0, points, values)
//                  << std::endl;
      }
      rhs[cellIndex(cellId)] = div / dt;
    }
    {
      std::vector<Point2d> points(neighbors.size() + 1);
      points[0] = center;
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0)
          points[i + 1] = _structure->cellCenterPosition(neighbors[i].id);
        else {
          switch (neighbors[i].side) {
            case Definitions::Side::LEFT:
              points[i + 1] = Point2d(center.x() - cellRegionSize.x(),
                                      center.y());
              break;
            case Definitions::Side::RIGHT:
              points[i + 1] = Point2d(center.x() + cellRegionSize.x(),
                                      center.y());
              break;
            case Definitions::Side::BOTTOM:
              points[i + 1] = Point2d(center.x(),
                                      center.y() - cellRegionSize.y());
              break;
            case Definitions::Side::TOP:
              points[i + 1] = Point2d(center.x(),
                                      center.y() + cellRegionSize.y());
              break;
            default: THROW(false, "PicSimulation2: invalid neighbor side");
          }
        }
      }
      std::vector<double> weights = rbf.laplacianWeights(points);
      double centerWeight = weights[0];
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0
            && cellMaterialField[neighbors[i].id] == Definitions::Material::AIR)
          continue;
        if (neighbors[i].id >= 0)
          triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
                                                       cellIndex(neighbors[i].id),
                                                       weights[i + 1]));
        else
          centerWeight += weights[i + 1];
      }
      triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
                                                   cellIndex(cellId),
                                                   centerWeight));
    }
  });
  triplets.emplace_back(Eigen::Triplet<double>(indexMap.size(),
                                               indexMap.size(),
                                               1.));
//  std::cerr << "RHS \n" << rhs << std::endl;
  Eigen::SparseMatrix<double> A(indexMap.size() + 1, indexMap.size() + 1);
  A.setFromTriplets(triplets.begin(), triplets.end());
//  std::cerr << A << std::endl;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver;
  solver.compute(A);
  Eigen::VectorXd pressure = solver.solve(rhs);
//  std::cerr << pressure << std::endl;
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  // correct velocities
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      pressureField[cellId] = pressure[cellIndex(cellId)];
  });
  _structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    auto cells = _structure->faceCells(faceId);
    std::vector<double> pressureValues(2, 0.);
    std::vector<Point2d> points(2);
    for (size_t i = 0; i < 2; i++) {
      THROW(cells[i].id >= 0,
            "PicSimulation2::computePressure invalid face cell neighbor.")
      points[i] = _structure->cellCenterPosition(cells[i].id);
      if (cellMaterialField[cells[i].id] == Definitions::Material::FLUID)
        pressureValues[i] = pressure[cellIndex(cells[i].id)];
    }
    double grad = 0.;
    switch (_structure->faceOrientation(faceId)) {
      case Definitions::Orientation::VERTICAL:
        grad = rbf.gradientAt(_structure->faceCenterPosition(faceId), 0,
                              points,
                              pressureValues);
//        std::cerr << "GRADx " << grad * dt << std::endl;
        break;
      case Definitions::Orientation::HORIZONTAL:
        grad = rbf.gradientAt(_structure->faceCenterPosition(faceId), 1,
                              points,
                              pressureValues);
//        std::cerr << "GRADy " << grad * dt << std::endl;
        break;
      default: THROW(false,
                     "PicSimulation2::computePressure wrong face orientation");
    }
//    if (_structure->faceOrientation(faceId)
//        == Definitions::Orientation::VERTICAL)
//      std::cerr << velocityField[faceId] << " -> ";
    velocityField[faceId] -= grad * dt;
//    if (_structure->faceOrientation(faceId)
//        == Definitions::Orientation::VERTICAL)
//      std::cerr << velocityField[faceId] << std::endl;
  });
}

void PicSimulation2::advectParticles(double dt) {
  std::cerr << "advect " << dt << std::endl;
  _particles->iterateParticles([&](unsigned int i, Point2d p) {
    // EULER
    auto v =
        Vector2d(_particles->getScalarProperty(particleVelocityFieldIds[0], i),
                 _particles->getScalarProperty(particleVelocityFieldIds[1], i));
    auto np = _integrator->integrate(p, v, dt);
    // Check if point is still inside domain
    auto region = _structure->domainRegion();
    if (np.x() > region.upper().x() || np.x() < region.lower().x()
        || np.y() > region.upper().y() || np.y() < region.lower().y())
      _particles->remove(i);
    else {
      np = Point2d(clamp(np.x(), 0.0, 0.99999), clamp(np.y(), 0.0, 0.99999));
      _particles->setPosition(i, np);
      THROW(_structure->domainRegion().contains(np),
            "PicSimulation2::advectParticles Particles outside domain!");
    }
  });
}

double PicSimulation2::cfl() const {
  NOT_IMPLEMENTED();
  return 0.;
}

void PicSimulation2::defineBoundaryConditions(Definitions::Side side,
                                              Vector2d vel) {

}

}