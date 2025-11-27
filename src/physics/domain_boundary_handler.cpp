#include <physics/domain_boundary_handler.h>
#include <physics/simulation_domain.h>

namespace furoo {

DomainBoundaryHandler2::DomainBoundaryHandler2()
    : DomainBoundaryHandler2(BoundaryType::CUSTOM) {}

DomainBoundaryHandler2::DomainBoundaryHandler2(BoundaryType bt) {
  for (size_t i = 0; i < 4; i++)
    _boundaries[i] = bt;
}

DomainBoundaryHandler2::DomainBoundaryHandler2(BoundaryType left,
                                               BoundaryType bottom,
                                               BoundaryType right,
                                               BoundaryType top) {
  _boundaries[0] = left;
  _boundaries[1] = bottom;
  _boundaries[2] = right;
  _boundaries[3] = top;
}

void DomainBoundaryHandler2::markFaces(const SimulationDomain2 *domain) {
  for (size_t f = 0; f < domain->faceCount(); f++) {
    _faceBoundaries[f] = BoundaryType::NONE;
    auto neighborhood = domain->faceNeighborhood(f);
    for (size_t n = 0; n < 2; n++)
      if (neighborhood[n].id < 0)
        switch (neighborhood[n].orientation) {
        case SimulationDomain2::Orientation::LEFT:_faceBoundaries[f] = _boundaries[0];
          break;
        case SimulationDomain2::Orientation::BOTTOM:_faceBoundaries[f] = _boundaries[1];
          break;
        case SimulationDomain2::Orientation::RIGHT:_faceBoundaries[f] = _boundaries[2];
          break;
        case SimulationDomain2::Orientation::TOP:_faceBoundaries[f] = _boundaries[3];
          break;
        default: NOT_IMPLEMENTED();
        }
    if (neighborhood[0].id >= 0 && neighborhood[1].id >= 0) {
      if ((cellType(neighborhood[0].id) == MaterialType::FLUID &&
          cellType(neighborhood[1].id) == MaterialType::AIR) ||
          (cellType(neighborhood[1].id) == MaterialType::FLUID &&
              cellType(neighborhood[0].id) == MaterialType::AIR)) {
        _faceBoundaries[f] = BoundaryType::DIRICHLET;
      }
    }
    if (neighborhood[0].id < 0)
      _faces[f] = (_faceBoundaries[f] == BoundaryType::DIRICHLET) ? cellType(neighborhood[1].id) : MaterialType::SOLID;
    else if (neighborhood[1].id < 0)
      _faces[f] = (_faceBoundaries[f] == BoundaryType::DIRICHLET) ? cellType(neighborhood[0].id) : MaterialType::SOLID;
    else {
      auto t1 = cellType(neighborhood[0].id);
      auto t2 = cellType(neighborhood[1].id);
      if (t1 == MaterialType::SOLID || t2 == MaterialType::SOLID)
        _faces[f] = MaterialType::SOLID;
      else if (t1 == MaterialType::FLUID || t2 == MaterialType::FLUID)
        _faces[f] = MaterialType::FLUID;
      else
        _faces[f] = MaterialType::AIR;
    }
  }
  // REMARK CELLS
  for (size_t i = 0; i < domain->cellCount(); i++) {
    if (_cells[i] == MaterialType::AIR) {
      auto faces = domain->cellFaces(i);
      bool allFluid = true;
      for (auto f : faces)
        if (_faces[f.id] == MaterialType::AIR) {
          allFluid = false;
          break;
        }
      if (allFluid)
        _cells[i] = MaterialType::FLUID;
    }
  }
}

void DomainBoundaryHandler2::markCells(const SimulationDomain2 *domain,
                                       const ParticleSystem2 *ps) {
  for (auto &t : _cells)
    t = MaterialType::AIR;
  ps->iterateParticles([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(i);
    int cid = domain->cellId(p);
    if (cid >= 0)
      _cells[cid] = MaterialType::FLUID;
  });
  // TODO iterate per cell instead ?
}

DomainBoundaryHandler2::Intersection
DomainBoundaryHandler2::intersect(const SimulationDomain2 *domain,
                                  Point2d initialPosition,
                                  Point2d finalPosition) const {
  Vector2d n;
  Point2d p;
  bool iv = false;
  auto r = domain->region();
  if (!r.inside(finalPosition)) {
    iv = true;
    r.intersect(initialPosition, (finalPosition - initialPosition).normalized(),
                p, n);
  }
  return Intersection(p, n, iv);
}

DomainBoundaryHandler2::~DomainBoundaryHandler2() {}

DomainBoundaryHandler2::BoundaryType
DomainBoundaryHandler2::faceBoundaryType(size_t faceId) const {
  ASSERT_FATAL(faceId < _faceBoundaries.size());
  return _faceBoundaries[faceId];
}

DomainBoundaryHandler2::MaterialType
DomainBoundaryHandler2::faceType(size_t faceId) const {
  ASSERT_FATAL(faceId < _faces.size());
  return _faces[faceId];
}

double DomainBoundaryHandler2::faceValue(size_t faceId) const {
  ASSERT_FATAL(faceId < _boundaryFacesValues.size());
  return _boundaryFacesValues[faceId];
}

void DomainBoundaryHandler2::setBoundaryFaceValue(size_t faceId, double v) {
  ASSERT_FATAL(faceId < _boundaryFacesValues.size());
  _boundaryFacesValues[faceId] = v;
}

DomainBoundaryHandler2::MaterialType
DomainBoundaryHandler2::cellType(size_t cellId) const {
  ASSERT_FATAL(cellId < _cells.size());
  return _cells[cellId];
}

void DomainBoundaryHandler2::resize(size_t cellCount, size_t faceCount) {
  _faceBoundaries.resize(faceCount, BoundaryType::NONE);
  _cells.resize(cellCount, MaterialType::CUSTOM);
  _faces.resize(faceCount, MaterialType::CUSTOM);
  _boundaryFacesValues.resize(faceCount, 0.);
}

size_t DomainBoundaryHandler2::dirichletBoundaryCount() const {
  size_t d = 0;
  for (auto bd : _faceBoundaries) {
    if (bd == BoundaryType::DIRICHLET)
      d++;
  }
  return d;
}

unsigned DomainBoundaryHandler2::fluidCellCount() const {
  unsigned count = 0;
  for (auto c : _cells)
    if (c == MaterialType::FLUID)
      count++;
  return count;
}

void DomainBoundaryHandler2::printFaces() {
  for (size_t i = 0; i < _faceBoundaries.size(); i++) {
    std::cout << "(" << i << ", " << _faceBoundaries[i] << ") ";
    /* code */
  }
  std::cout << '\n';
}

void DomainBoundaryHandler2::iterateBoundaryFaces(
    const std::function<void(size_t, BoundaryType)> &f) const {
  for (size_t i = 0; i < _faceBoundaries.size(); i++) {
    if (_faceBoundaries[i] != BoundaryType::NONE)
      f(i, _faceBoundaries[i]);
  }
}

void DomainBoundaryHandler2::iterateBoundaryCells(
    const SimulationDomain2 *domain, MaterialType t,
    const std::function<void(size_t)> &f) const {
  for (size_t i = 0; i < _cells.size(); i++) {
    if (_cells[i] != t)
      continue;
    auto neighborhood = domain->cellNeighborhood(i);
    for (auto n : neighborhood)
      if (n.id < 0 || _cells[n.id] != t) {
        f(i);
        break;
      }
  }
}

void DomainBoundaryHandler2::setFaceType(unsigned int faceId, BoundaryType t) {
  ASSERT_FATAL(_faceBoundaries.size() > faceId);
  _faceBoundaries[faceId] = t;
}

void DomainBoundaryHandler2::setCellType(unsigned int cellId, MaterialType t) {
  ASSERT_FATAL(_cells.size() > cellId);
  _cells[cellId] = t;
}

} // namespace furoo
