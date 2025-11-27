// DEPRECATED
#include <blas/interpolation.h>
#include <common/debug.h>
#include <geometry/numeric.h>
#include <limits>
#include <physics/simulation_regular_domain.h>
#include <queue>

namespace furoo {

SimRegularDomain2::SimRegularDomain2() : SimRegularDomain2(8) {}

SimRegularDomain2::SimRegularDomain2(size_t
                                     n) {
  double hx = 1.0 / n;

  _cellToGrid = Transform2D::scale(n - 1., n - 1.) *
      Transform2D(BBox2D(Point2d(hx / 2., hx / 2.),
                         Point2d(1. - hx / 2., 1 - hx / 2.)));
  _cellToWorld = _cellToGrid.inverse();
  _uFaceToGrid =
      Transform2D::scale(static_cast<double>(n), n - 1.) *
          Transform2D(BBox2D(Point2d(0., hx / 2.), Point2d(1., 1 - hx / 2.)));
  _uFaceToWorld = _uFaceToGrid.inverse();
  _vFaceToGrid =
      Transform2D::scale(n - 1., static_cast<double>(n)) *
          Transform2D(BBox2D(Point2d(hx / 2., 0.), Point2d(1 - hx / 2., 1.)));
  _vFaceToWorld = _vFaceToGrid.inverse();

  _grid = RegularGrid(n, _cellToWorld);
  this->
      _interpolant = SimulationDomain2::InterpolationMethod::BICUBIC;
}

RegularGrid SimRegularDomain2::grid() const { return _grid; }

unsigned SimRegularDomain2::uFacesCount() const {
  auto s = _grid.gridSize();
  return (s.x() + 1) * s.y();
}

int SimRegularDomain2::faceIdFromFaceId(
    unsigned faceId, SimulationDomain2::Orientation o) const {
  auto s = _grid.gridSize();
  unsigned nx = s.x(), ny = s.y();
  unsigned fid = faceId;
  if (faceId < uFacesCount())
    nx++;
  else {
    ny++;
    fid -= uFacesCount();
  }
  switch (o) {
    case SimulationDomain2::Orientation::LEFT:
      return (fid % nx > 0) ? faceId - 1 : -1;
      break;
    case SimulationDomain2::Orientation::RIGHT:
      return (fid % nx < nx - 1) ? faceId + 1 : -1;
      break;
    case SimulationDomain2::Orientation::BOTTOM:
      return (fid / nx > 0) ? faceId - nx : -1;
      break;
    case SimulationDomain2::Orientation::TOP:
      return (fid / nx < ny - 1) ? faceId + nx : -1;
      break;
    case SimulationDomain2::Orientation::OTHER: ASSERT_FATAL(false);
      break;
  }
  ASSERT_FATAL(false);
  return -1;
}

int SimRegularDomain2::faceIdFromCellId(
    unsigned cellId, SimulationDomain2::Orientation o) const {
  auto s = _grid.gridSize();
  int i = cellId / s.x();
  switch (o) {
    case SimulationDomain2::Orientation::LEFT:return cellId + i;
      break;
    case SimulationDomain2::Orientation::RIGHT:return cellId + i + 1;
      break;
    case SimulationDomain2::Orientation::BOTTOM:
      return cellId + (s.x() + 1) * s.y();
      break;
    case SimulationDomain2::Orientation::TOP:
      return cellId + (s.x() + 1) * s.y() + s.x();
      break;
    case SimulationDomain2::Orientation::OTHER: ASSERT_FATAL(false);
      break;
  }
  ASSERT_FATAL(false);
  return -1;
}

int SimRegularDomain2::cellIdFromFaceId(
    unsigned faceId, SimulationDomain2::Orientation o) const {
  auto s = _grid.gridSize();
  unsigned _uFacesCount = uFacesCount();
  if (faceId < _uFacesCount) {
    unsigned j = faceId % (s.x() + 1);
    unsigned i = faceId / (s.x() + 1);
    if (o == SimulationDomain2::Orientation::LEFT)
      return (j > 0) ? faceId - i - 1 : -1;
    else if (o == SimulationDomain2::Orientation::RIGHT)
      return (j < s.x()) ? faceId - i : -1;
    else {
      ASSERT_FATAL(false);
      return -1;
    }
  } else {
    faceId -= _uFacesCount;
    unsigned i = faceId / s.x();
    if (o == SimulationDomain2::Orientation::BOTTOM)
      return (i > 0) ? faceId - s.x() : -1;
    else if (o == SimulationDomain2::Orientation::TOP)
      return (i < s.y()) ? faceId : -1;
    else {
      ASSERT_FATAL(false);
      return -1;
    }
  }
  ASSERT_FATAL(false);
  return -1;
}

int SimRegularDomain2::cellIdFromCellId(
    unsigned cellId, SimulationDomain2::Orientation o) const {
  auto s = _grid.gridSize();
  switch (o) {
    case SimulationDomain2::Orientation::LEFT:
      return (cellId % s.x() > 0) ? cellId - 1 : -1;
      break;
    case SimulationDomain2::Orientation::RIGHT:
      return (cellId % s.x() < s.x() - 1) ? cellId + 1 : -1;
      break;
    case SimulationDomain2::Orientation::BOTTOM:
      return (cellId / s.x() > 0) ? cellId - s.x() : -1;
      break;
    case SimulationDomain2::Orientation::TOP:
      return (cellId / s.x() < s.y() - 1) ? cellId + s.x() : -1;
      break;
    case SimulationDomain2::Orientation::OTHER: ASSERT_FATAL(false);
      break;
  }
  ASSERT_FATAL(false);
  return -1;
}

unsigned SimRegularDomain2::cellIdFromIJ(Point2i ij) const {
  auto s = _grid.gridSize();
  return ij[1] * s.x() + ij[0];
}

unsigned SimRegularDomain2::uFaceIdFromIJ(Point2i ij) const {
  auto s = _grid.gridSize();
  return ij[1] * (s.x() + 1) + ij[0];
}

unsigned SimRegularDomain2::vFaceIdFromIJ(Point2i ij) const {
  auto s = _grid.gridSize();
  return uFacesCount() + ij[1] * s.y() + ij[0];
}

Point2i SimRegularDomain2::ijFromCellId(unsigned cellId) const {
  auto s = _grid.gridSize();
  int j = cellId / s.x();
  int i = cellId % s.x();
  return Point2i(i, j);
}

Point2i SimRegularDomain2::ijFromFaceId(unsigned faceId) const {
  auto s = _grid.gridSize();
  unsigned ux = s.x() + 1;
  unsigned uFaceCount = ux * s.y();
  if (faceId < uFaceCount)
    return Point2i(faceId % ux, faceId / ux);
  faceId -= uFaceCount;
  return Point2i(faceId % s.x(), faceId / s.x());
}

int SimRegularDomain2::cellId(Point2d p) const {
  auto gp = _cellToGrid(p);
  int i = round2Int(gp.x());
  int j = round2Int(gp.y());
  auto s = _grid.gridSize();
  if (i < 0 || i >= static_cast<int>(s.x()) || j < 0 ||
      j >= static_cast<int>(s.y()))
    return -1;
  return static_cast<int>(j * s.x() + i);
}

void SimRegularDomain2::setBoundaryHandler(DomainBoundaryHandler2 *bh) {
  this->_boundaryHandler = bh;
  bh->resize(cellCount(), faceCount());
}

void SimRegularDomain2::setInterpolationMethod(
    InterpolationMethod interpolant) {
  this->_interpolant = interpolant;
}

Vector2d SimRegularDomain2::smallestCellSize() const { return _grid.spacing(); }

BBox2D SimRegularDomain2::region() const {
  auto s = _grid.gridSize();
  return BBox2D(_cellToWorld(Point2d(-.5, -.5)),
                _cellToWorld(Point2d((s.x() - 1) + .5, (s.y() - 1) + .5)));
}

BBox2D SimRegularDomain2::cellRegion(size_t cellId) const {
  ASSERT_FATAL(cellId < _grid.cellCount());
  auto hh = .5 * _grid.spacing();
  Point2d c = cellCenterPosition(cellId);
  return BBox2D(c - hh, c + hh);
}

Point2d SimRegularDomain2::cellCenterPosition(size_t cellId) const {
  return _cellToWorld(ijFromCellId(cellId));
}

Point2d SimRegularDomain2::faceCenterPosition(size_t faceId) const {
  if (faceId < uFacesCount())
    return _uFaceToWorld(ijFromFaceId(faceId));
  return _vFaceToWorld(ijFromFaceId(faceId));
}

double &SimRegularDomain2::faceCenteredScalar(size_t fieldId, size_t faceId) {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  ASSERT_FATAL(faceId < faceCount());
  if (faceId < _uFaceFields[fieldId].size())
    return _uFaceFields[fieldId][faceId];
  return _vFaceFields[fieldId][faceId - _uFaceFields[fieldId].size()];
}

std::vector<SimulationDomain2::NeighborOrientation>
SimRegularDomain2::faceNeighborhood(size_t faceId) const {
  std::vector<SimulationDomain2::NeighborOrientation> r;
  if (faceId < uFacesCount()) {
    r.emplace_back(
        cellIdFromFaceId(faceId, SimulationDomain2::Orientation::LEFT),
        SimulationDomain2::Orientation::LEFT);
    r.emplace_back(
        cellIdFromFaceId(faceId, SimulationDomain2::Orientation::RIGHT),
        SimulationDomain2::Orientation::RIGHT);
  } else {
    r.emplace_back(
        cellIdFromFaceId(faceId, SimulationDomain2::Orientation::BOTTOM),
        SimulationDomain2::Orientation::BOTTOM);
    r.emplace_back(
        cellIdFromFaceId(faceId, SimulationDomain2::Orientation::TOP),
        SimulationDomain2::Orientation::TOP);
  }
  return r;
}

std::vector<SimulationDomain2::NeighborOrientation>
SimRegularDomain2::cellNeighborhood(size_t cellId) const {
  std::vector<SimulationDomain2::NeighborOrientation> r;
  r.emplace_back(cellIdFromCellId(cellId, SimulationDomain2::Orientation::LEFT),
                 SimulationDomain2::Orientation::LEFT);
  r.emplace_back(
      cellIdFromCellId(cellId, SimulationDomain2::Orientation::RIGHT),
      SimulationDomain2::Orientation::RIGHT);
  r.emplace_back(
      cellIdFromCellId(cellId, SimulationDomain2::Orientation::BOTTOM),
      SimulationDomain2::Orientation::BOTTOM);
  r.emplace_back(cellIdFromCellId(cellId, SimulationDomain2::Orientation::TOP),
                 SimulationDomain2::Orientation::TOP);
  return r;
}

std::vector<SimulationDomain2::NeighborOrientation> SimRegularDomain2::cellFaces(
    unsigned cellId) const {
  ASSERT_FATAL(cellId < cellCount());
  std::vector<SimulationDomain2::NeighborOrientation> r;
  r.emplace_back(faceIdFromCellId(cellId, SimulationDomain2::Orientation::LEFT),
                 SimulationDomain2::Orientation::LEFT);
  r.emplace_back(
      faceIdFromCellId(cellId, SimulationDomain2::Orientation::RIGHT),
      SimulationDomain2::Orientation::RIGHT);
  r.emplace_back(
      faceIdFromCellId(cellId, SimulationDomain2::Orientation::BOTTOM),
      SimulationDomain2::Orientation::BOTTOM);
  r.emplace_back(faceIdFromCellId(cellId, SimulationDomain2::Orientation::TOP),
                 SimulationDomain2::Orientation::TOP);
  return r;
}

void SimRegularDomain2::iterateStencilAtCellCenter(
    size_t id,
    const std::function<void(int nid, double w, Point2d p,
                             SimulationDomain2::Orientation o)> &f) const {
  Vector2d h = _grid.spacing();
  double h2 = h[0] * h[1];
  Point2i ij = ijFromCellId(id);
  Point2d nc;
  size_t neumannCount = 0, dirichletCount = 0;

  // LEFT =============================================
  int neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::LEFT);
  double neighW = 0.0;
  DomainBoundaryHandler2::BoundaryType neighBd =
      this->_boundaryHandler->faceBoundaryType(neighId);
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
    nc = _cellToWorld(Point2d((ij[0] - 0.5), ij[1]));
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
    nc = _cellToWorld(Point2d((ij[0] - 1), ij[1]));
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::LEFT),
    -neighW / h2,
    nc,
    SimulationDomain2::Orientation::LEFT);

  // BOTTOM =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::BOTTOM);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
    nc = _cellToWorld(Point2d(ij[0], (ij[1] - 0.5)));
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
    nc = _cellToWorld(Point2d(ij[0], (ij[1] - 1)));
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::BOTTOM), -neighW / h2,
    nc, SimulationDomain2::Orientation::BOTTOM);

  // RIGHT =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::RIGHT);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
    nc = _cellToWorld(Point2d((ij[0] + 0.5), ij[1]));
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
    nc = _cellToWorld(Point2d((ij[0] + 1), ij[1]));
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::RIGHT), -neighW / h2,
    nc, SimulationDomain2::Orientation::RIGHT);

  // TOP =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::TOP);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
    nc = _cellToWorld(Point2d(ij[0], (ij[1] + 0.5)));
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
    nc = _cellToWorld(Point2d(ij[0], (ij[1] + 1)));
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::TOP), -neighW / h2, nc,
    SimulationDomain2::Orientation::TOP);

  // CENTER
  f(id, (4.0 - neumannCount) / h2,
    _cellToWorld(Point2d((ij[0]), ij[1])),
    SimulationDomain2::Orientation::OTHER);
}

void SimRegularDomain2::iterateStencilAtCellCenter(
    size_t id,
    const std::function<void(int nid, double w,
                             DomainBoundaryHandler2::BoundaryType, double)> &f)
const {
  Vector2d h = _grid.spacing();
  double h2 = h[0] * h[1];
  size_t neumannCount = 0, dirichletCount = 0;

  // LEFT =============================================
  size_t neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::LEFT);
  double neighW = 0.0;
  DomainBoundaryHandler2::BoundaryType neighBd =
      this->_boundaryHandler->faceBoundaryType(neighId);
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::LEFT), -neighW / h2,
    neighBd, this->_boundaryHandler->faceValue(neighId));

  // BOTTOM =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::BOTTOM);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::BOTTOM), -neighW / h2,
    neighBd, this->_boundaryHandler->faceValue(neighId));

  // RIGHT =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::RIGHT);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::RIGHT), -neighW / h2,
    neighBd, this->_boundaryHandler->faceValue(neighId));

  // TOP =============================================
  neighId = faceIdFromCellId(id, SimulationDomain2::Orientation::TOP);
  neighBd = this->_boundaryHandler->faceBoundaryType(neighId);
  neighW = 0;
  if (neighBd == DomainBoundaryHandler2::BoundaryType::DIRICHLET) {
    neighW = 0;
    dirichletCount++;
  } else {
    if (neighBd == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      neumannCount++;
    else
      neighW = 1;
  }
  f(cellIdFromCellId(id, SimulationDomain2::Orientation::TOP), -neighW / h2,
    neighBd, this->_boundaryHandler->faceValue(neighId));

  // CENTER
  f(id, (4.0 - neumannCount) / h2,
    DomainBoundaryHandler2::BoundaryType::NONE, 0.);
}

void SimRegularDomain2::iterateDxScalar(
    size_t fieldId, const std::function<void(double &val, Point2d pos)> &f) {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  for (unsigned i = 0; i < _uFaceFields[fieldId].size(); i++)
    f(_uFaceFields[fieldId][i], _uFaceToWorld(ijFromFaceId(i)));
}

void SimRegularDomain2::iterateDxScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  for (size_t i = 0; i < _uFaceFields[fieldId].size(); i++)
    f(i, _uFaceFields[fieldId][i]);
}

void SimRegularDomain2::iterateDyScalar(
    size_t fieldId, const std::function<void(double &val, Point2d pos)> &f) {
  auto _uFacesCount = uFacesCount();
  ASSERT_FATAL(fieldId < this->_vFaceFields.size());
  for (unsigned i = 0; i < _vFaceFields[fieldId].size(); i++)
    f(this->_vFaceFields[fieldId][i],
      _vFaceToWorld(ijFromFaceId(i + _uFacesCount)));
}

void SimRegularDomain2::iterateDyScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < _vFaceFields.size());
  for (size_t i = 0; i < _vFaceFields[fieldId].size(); i++)
    f(i + _uFaceFields[fieldId].size(), _vFaceFields[fieldId][i]);
}

void SimRegularDomain2::iterateCellScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  ASSERT_FATAL(fieldId < this->_cellCenteredScalarFields.size());
  for (unsigned i = 0; i < cellCount(); i++)
    f(this->_cellCenteredScalarFields[fieldId][i],
      _cellToWorld(ijFromCellId(i)));
}

void SimRegularDomain2::iterateCellScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < this->_cellCenteredScalarFields.size());
  for (unsigned i = 0; i < cellCount(); i++)
    f(i, this->_cellCenteredScalarFields[fieldId][i]);
}

void SimRegularDomain2::iterateCellByte(
    size_t fieldId, const std::function<void(unsigned char &)> &f) {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(f);
  NOT_IMPLEMENTED();
}

void SimRegularDomain2::iterateVertexScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(f);
}

size_t SimRegularDomain2::addFaceCenteredScalarField(double v) {
  Vector<size_t, 2> _size = _grid.gridSize();
  _uFaceFields.emplace_back(
      std::vector<double>(_size.x() * (_size.y() + 1), v));
  _vFaceFields.emplace_back(
      std::vector<double>((_size.x() + 1) * _size.y(), v));
  return _uFaceFields.size() - 1;
}

size_t SimRegularDomain2::addCellCenteredScalarField(double v) {
  _cellCenteredScalarFields.emplace_back(
      std::vector<double>(_grid.cellCount(), v));
  return _cellCenteredScalarFields.size() - 1;
}

size_t SimRegularDomain2::addCellCenteredByteField(unsigned char v) {
  _cellCenteredByteFields.emplace_back(
      std::vector<unsigned char>(_grid.cellCount(), v));
  return _cellCenteredByteFields.size() - 1;
}

size_t SimRegularDomain2::addVertexCenteredScalarField(double v) {
  UNUSED_VARIABLE(v);
  NOT_IMPLEMENTED();
  return 0;
}

size_t SimRegularDomain2::cellCount() const { return _grid.cellCount(); }

size_t SimRegularDomain2::faceCount() const {
  Vector<size_t, 2> _size = _grid.gridSize();
  return _size.x() * (_size.y() + 1) + (_size.x() + 1) * _size.y();
}

size_t SimRegularDomain2::vertexCount() const {
  NOT_IMPLEMENTED();
  return 0;
}

double SimRegularDomain2::sampleCellCenteredScalar(size_t fieldId,
                                                   Point2d target) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: {
      unsigned cells[4];
      Point2d p = cellCenteredGridCell(target, &cells[0]);
      std::vector<double> f;
      for (int i = 0; i < 4; i++)
        f.emplace_back(_cellCenteredScalarFields[fieldId][cells[i]]);
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return BilinearInterpolant().interpolateAt(p, points, f);
    }
    case SimulationDomain2::InterpolationMethod::BICUBIC: {
      unsigned cells[4];
      Point2d gp = cellCenteredGridCell(target, &cells[0]);
      std::vector<double> fx, fy, fxy, f;
      for (int i = 0; i < 4; i++) {
        fx.emplace_back(gradientXAtCellCenteredScalar(fieldId, cells[i]));
        fy.emplace_back(gradientYAtCellCenteredScalar(fieldId, cells[i]));
        fxy.emplace_back(gradientXYAtCellCenteredScalar(fieldId, cells[i]));
        f.emplace_back(_cellCenteredScalarFields[fieldId][cells[i]]);
      }
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return _bicubicInterpolant.interpolateAt(gp, f, fx, fy, fxy,
                                               _grid.spacing());
    }
    case SimulationDomain2::InterpolationMethod::CATMULL: {
      int sizeX = static_cast<int>(_grid.gridSize().x());
      int sizeY = static_cast<int>(_grid.gridSize().y());

      int i, j;
      double fx, fy;
      MonotonicCatmullRom::getFraction(target.x(), 0, sizeX - 1, &i, &fx);
      MonotonicCatmullRom::getFraction(target.y(), 0, sizeY - 1, &j, &fy);

      int is[4] = {
          std::max(i - 1, 0),
          i,
          std::min(i + 1, sizeX - 1),
          std::min(i + 2, sizeX - 1)
      };

      int js[4] = {
          std::max(j - 1, 0),
          j,
          std::min(j + 1, sizeY - 1),
          std::min(j + 2, sizeY - 1)

      };

      double values[4];

      for (int k = 0; k < 4; ++k) {
        values[k] = MonotonicCatmullRom::interpolate(
            _cellCenteredScalarFields[fieldId][cellIdFromIJ(Point2i(is[0],
                                                                    js[k]))],
            _cellCenteredScalarFields[fieldId][cellIdFromIJ(Point2i(is[1],
                                                                    js[k]))],
            _cellCenteredScalarFields[fieldId][cellIdFromIJ(Point2i(is[2],
                                                                    js[k]))],
            _cellCenteredScalarFields[fieldId][cellIdFromIJ(Point2i(is[3],
                                                                    js[k]))],
            fx);
      }

      return MonotonicCatmullRom::interpolate(values[0],
                                              values[1],
                                              values[2],
                                              values[3],
                                              fy);

    }
    case SimulationDomain2::InterpolationMethod::OTHER: {
      NOT_IMPLEMENTED();
      return 0.;
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;
}

double SimRegularDomain2::gradientXAtCellCenteredScalar(size_t fieldId,
                                                        size_t cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  size_t _n = _grid.gridSize()[0];
  Vector2d h = _grid.spacing();
  size_t i = cellId % _n;
  size_t j = cellId / _n;
  ASSERT_FATAL(i > 0 || i < _n - 1);
  if (i > 0 && i < _n - 1) {
    return (_cellCenteredScalarFields[fieldId][j * _n + i + 1] -
        _cellCenteredScalarFields[fieldId][j * _n + i - 1]) /
        (2. * h[0]);
  }
  if (i == 0)
    return (_cellCenteredScalarFields[fieldId][j * _n + i + 1] -
        _cellCenteredScalarFields[fieldId][j * _n + i]) /
        (h[0]);
  return (_cellCenteredScalarFields[fieldId][j * _n + i] -
      _cellCenteredScalarFields[fieldId][j * _n + i - 1]) /
      (h[0]);
}

double SimRegularDomain2::gradientYAtCellCenteredScalar(size_t fieldId,
                                                        size_t cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  size_t _n = _grid.gridSize()[0];
  size_t _m = _grid.gridSize()[1];
  Vector2d h = _grid.spacing();
  size_t i = cellId % _n;
  size_t j = cellId / _n;
  ASSERT_FATAL(j > 0 || j < _m - 1);
  if (j > 0 && j < _m - 1)
    return (_cellCenteredScalarFields[fieldId][(j + 1) * _n + i] -
        _cellCenteredScalarFields[fieldId][(j - 1) * _n + i]) /
        (2. * h[1]);
  if (j == 0)
    return (_cellCenteredScalarFields[fieldId][(j + 1) * _n + i] -
        _cellCenteredScalarFields[fieldId][j * _n + i]) /
        (h[1]);
  return (_cellCenteredScalarFields[fieldId][j * _n + i] -
      _cellCenteredScalarFields[fieldId][(j - 1) * _n + i]) /
      (h[1]);
}

double SimRegularDomain2::gradientXYAtCellCenteredScalar(size_t fieldId,
                                                         size_t cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  size_t _n = _grid.gridSize()[0];
  size_t _m = _grid.gridSize()[1];
  Vector2d h = _grid.spacing();
  size_t i = cellId % _n;
  size_t j = cellId / _n;
  ASSERT_FATAL(j > 0 || j < _m - 1);
  double dyl = 0., dyr = 0.;
  double hx = h[0];
  if (i < _n - 1 && i > 0) {
    dyl = gradientYAtCellCenteredScalar(fieldId, j * _n + i - 1);
    dyr = gradientYAtCellCenteredScalar(fieldId, j * _n + i + 1);
    hx *= 2;
  } else if (i == 0) {
    dyl = gradientYAtCellCenteredScalar(fieldId, j * _n + i);
    dyr = gradientYAtCellCenteredScalar(fieldId, j * _n + i + 1);
  } else {
    dyl = gradientYAtCellCenteredScalar(fieldId, j * _n + i - 1);
    dyr = gradientYAtCellCenteredScalar(fieldId, j * _n + i);
  }
  return (dyr - dyl) / hx;
}

double SimRegularDomain2::faceDivergentAtCellCenter(unsigned int fieldId,
                                                    unsigned int cellId) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  ASSERT_FATAL(cellId < cellCount());
  double gx = 0., gy = 0.;
  Vector2d h = _grid.spacing();
  auto uCount = uFacesCount();

  gx = (_uFaceFields[fieldId][faceIdFromCellId(
      cellId, SimulationDomain2::Orientation::RIGHT)] -
      _uFaceFields[fieldId][faceIdFromCellId(
          cellId, SimulationDomain2::Orientation::LEFT)]) /
      (h[0]);
  gy = (_vFaceFields[fieldId][faceIdFromCellId(
      cellId, SimulationDomain2::Orientation::TOP) -
      uCount] -
      _vFaceFields[fieldId]
      [faceIdFromCellId(cellId,
                        SimulationDomain2::Orientation::BOTTOM) -
          uCount]) /
      (h[1]);

  return gx + gy;
}

double SimRegularDomain2::cellGradientAtFaceCenter(unsigned int fieldId,
                                                   unsigned int faceId) const {
  std::vector<SimulationDomain2::NeighborOrientation> cellNeighbors =
      faceNeighborhood(faceId);
//  if (cellNeighbors[0].id == -1 || cellNeighbors[1].id == -1)
//    return 0;
  double h;
  if (cellNeighbors[0].orientation == SimulationDomain2::Orientation::LEFT) {
    h = _grid.spacing()[0];
  } else {
    h = _grid.spacing()[1];
  }
  return (_cellCenteredScalarFields[fieldId][cellNeighbors[1].id] -
      _cellCenteredScalarFields[fieldId][cellNeighbors[0].id]) /
      h;
}

Vector2d SimRegularDomain2::sampleFaceCenteredScalar(size_t fieldId,
                                                     Point2d target) const {
  return Vector2d(sampleUFaceField(fieldId, target),
                  sampleVFaceField(fieldId, target));;
}

double SimRegularDomain2::sampleVFaceField(size_t fieldId,
                                           Point2d target) const {
  ASSERT_FATAL(fieldId < _vFaceFields.size());
  int uCount = uFacesCount();
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: {
      unsigned cells[4];
      std::vector<double> f(4, 0.0);
      Point2d p = vFaceCenteredGridCell(target, &cells[0], &f[0], fieldId);
      //std::vector<double> f;
      //for (int i = 0; i < 4; i++)
      //  f.emplace_back(_vFaceFields[fieldId][cells[i] - uCount]);
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return BilinearInterpolant().interpolateAt(p, points, f);
    }
    case SimulationDomain2::InterpolationMethod::BICUBIC: {
      unsigned cells[4];
      Point2d gp = vFaceCenteredGridCell(target, &cells[0]);
      std::vector<double> fx, fy, fxy, f;
      for (int i = 0; i < 4; i++) {
        fx.emplace_back(gradientXAtFaceCenteredScalar(fieldId, cells[i]));
        fy.emplace_back(gradientYAtFaceCenteredScalar(fieldId, cells[i]));
        fxy.emplace_back(gradientXYAtFaceCenteredScalar(fieldId, cells[i]));
        f.emplace_back(_vFaceFields[fieldId][cells[i] - uCount]);
      }
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return _bicubicInterpolant.interpolateAt(gp, f, fx, fy, fxy,
                                               _grid.spacing());
    }
    case SimulationDomain2::InterpolationMethod::RBF: {
      NOT_IMPLEMENTED();
      return 0.;
    }
    case SimulationDomain2::InterpolationMethod::OTHER: {
      NOT_IMPLEMENTED();
      return 0.;
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;
}

double SimRegularDomain2::sampleUFaceField(size_t fieldId,
                                           Point2d target) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: {
      unsigned cells[4];
      std::vector<double> f(4, 0.0);
      Point2d p = uFaceCenteredGridCell(target, &cells[0], &f[0], fieldId);
      //std::vector<double> f;
      //for (int i = 0; i < 4; i++)
      //  f.emplace_back(_uFaceFields[fieldId][cells[i]]);
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return BilinearInterpolant().interpolateAt(p, points, f);
    }
    case SimulationDomain2::InterpolationMethod::BICUBIC: {
      unsigned cells[4];
      Point2d gp = uFaceCenteredGridCell(target, &cells[0]);
      std::vector<double> fx, fy, fxy, f;
      for (int i = 0; i < 4; i++) {
        fx.emplace_back(gradientXAtFaceCenteredScalar(fieldId, cells[i]));
        fy.emplace_back(gradientYAtFaceCenteredScalar(fieldId, cells[i]));
        fxy.emplace_back(gradientXYAtFaceCenteredScalar(fieldId, cells[i]));
        f.emplace_back(_uFaceFields[fieldId][cells[i]]);
      }
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return _bicubicInterpolant.interpolateAt(gp, f, fx, fy, fxy,
                                               _grid.spacing());
    }
    case SimulationDomain2::InterpolationMethod::RBF: {
      NOT_IMPLEMENTED();
      return 0.;
    }
    case SimulationDomain2::InterpolationMethod::OTHER: {
      NOT_IMPLEMENTED();
      return 0.;
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;

  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(target);
  return 0.;
}

double SimRegularDomain2::gradientXAtFaceCenteredScalar(size_t fieldId,
                                                        size_t faceId) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size() && fieldId < _vFaceFields.size());
  auto size = _grid.gridSize();
  auto h = _grid.spacing();
  if (faceId < uFacesCount()) {
    size_t n = size.x() + 1;
    size_t i = faceId % n;
    size_t j = faceId / n;
    if (i > 0 && i < size.x())
      return (_uFaceFields[fieldId][j * n + i + 1] -
          _uFaceFields[fieldId][j * n + i - 1]) /
          (2. * h[0]);
    if (i == 0)
      return (_uFaceFields[fieldId][j * n + i + 1] -
          _uFaceFields[fieldId][j * n + i]) /
          (h[0]);
    return (_uFaceFields[fieldId][j * n + i] -
        _uFaceFields[fieldId][j * n + i - 1]) /
        (h[0]);
  } else {
    faceId -= uFacesCount();
    size_t n = size.x();
    size_t i = faceId % n;
    size_t j = faceId / n;
    if (i > 0 && i < size.x() - 1)
      return (_vFaceFields[fieldId][j * n + i + 1] -
          _vFaceFields[fieldId][j * n + i - 1]) /
          (2. * h[0]);
    if (i == 0)
      return (_vFaceFields[fieldId][j * n + i + 1] -
          _vFaceFields[fieldId][j * n + i]) /
          (h[0]);
    return (_vFaceFields[fieldId][j * n + i] -
        _vFaceFields[fieldId][j * n + i - 1]) /
        (h[0]);
  }
}

double SimRegularDomain2::gradientYAtFaceCenteredScalar(size_t fieldId,
                                                        size_t faceId) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size() && fieldId < _vFaceFields.size());
  auto size = _grid.gridSize();
  auto h = _grid.spacing();
  if (faceId < uFacesCount()) {
    size_t n = size.x() + 1;
    size_t _m = size.y();
    size_t i = faceId % n;
    size_t j = faceId / n;
    if (j > 0 && j < _m - 1)
      return (_uFaceFields[fieldId][(j + 1) * n + i] -
          _uFaceFields[fieldId][(j - 1) * n + i]) /
          (2. * h[1]);
    if (j == 0)
      return (_uFaceFields[fieldId][(j + 1) * n + i] -
          _uFaceFields[fieldId][j * n + i]) /
          (h[1]);
    return (_uFaceFields[fieldId][j * n + i] -
        _uFaceFields[fieldId][(j - 1) * n + i]) /
        (h[1]);
  }
  faceId -= uFacesCount();
  size_t n = size.x();
  size_t _m = size.y() + 1;
  size_t i = faceId % n;
  size_t j = faceId / n;
  if (j > 0 && j < _m - 1)
    return (_vFaceFields[fieldId][(j + 1) * n + i] -
        _vFaceFields[fieldId][(j - 1) * n + i]) /
        (2. * h[1]);
  if (j == 0)
    return (_vFaceFields[fieldId][(j + 1) * n + i] -
        _vFaceFields[fieldId][j * n + i]) /
        (h[1]);
  return (_vFaceFields[fieldId][j * n + i] -
      _vFaceFields[fieldId][(j - 1) * n + i]) /
      (h[1]);
}

double SimRegularDomain2::gradientXYAtFaceCenteredScalar(size_t fieldId,
                                                         size_t faceId) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size() && fieldId < _vFaceFields.size());
  auto size = _grid.gridSize();
  auto h = _grid.spacing();
  if (faceId < uFacesCount()) {
    size_t _n = size.x() + 1;
    size_t _m = size.y();
    size_t i = faceId % _n;
    size_t j = faceId / _n;
    ASSERT_FATAL(j > 0 || j < _m - 1);
    double dyl = 0., dyr = 0.;
    double hx = h[0];
    if (i < _n - 1 && i > 0) {
      dyl = gradientYAtFaceCenteredScalar(fieldId, j * _n + i - 1);
      dyr = gradientYAtFaceCenteredScalar(fieldId, j * _n + i + 1);
      hx *= 2;
    } else if (i == 0) {
      dyl = gradientYAtFaceCenteredScalar(fieldId, j * _n + i);
      dyr = gradientYAtFaceCenteredScalar(fieldId, j * _n + i + 1);
    } else {
      dyl = gradientYAtFaceCenteredScalar(fieldId, j * _n + i - 1);
      dyr = gradientYAtFaceCenteredScalar(fieldId, j * _n + i);
    }
    return (dyr - dyl) / hx;
  }
  size_t uCount = uFacesCount();
  faceId -= uCount;
  size_t _n = size.x();
  size_t _m = size.y() + 1;
  size_t i = faceId % _n;
  size_t j = faceId / _n;
  ASSERT_FATAL(j > 0 || j < _m - 1);
  double dyl = 0., dyr = 0.;
  double hx = h[0];
  if (i < _n - 1 && i > 0) {
    dyl = gradientYAtFaceCenteredScalar(fieldId, j * _n + i - 1 + uCount);
    dyr = gradientYAtFaceCenteredScalar(fieldId, j * _n + i + 1 + uCount);
    hx *= 2;
  } else if (i == 0) {
    dyl = gradientYAtFaceCenteredScalar(fieldId, uCount + j * _n + i);
    dyr = gradientYAtFaceCenteredScalar(fieldId, uCount + j * _n + i + 1);
  } else {
    dyl = gradientYAtFaceCenteredScalar(fieldId, uCount + j * _n + i - 1);
    dyr = gradientYAtFaceCenteredScalar(fieldId, uCount + j * _n + i);
  }
  return (dyr - dyl) / hx;
}

void SimRegularDomain2::sampleDxScalar(size_t fieldId,
                                       const std::vector<Point2d> &targets,
                                       LinearVector &values) const {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(targets);
  UNUSED_VARIABLE(values);
  NOT_IMPLEMENTED();
}

void SimRegularDomain2::sampleDyScalar(size_t fieldId,
                                       const std::vector<Point2d> &targets,
                                       LinearVector &values) const {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(targets);
  UNUSED_VARIABLE(values);
  NOT_IMPLEMENTED();
}

void SimRegularDomain2::extrapolateToDomain(size_t fieldId) {
  unsigned maxUnsigned = std::numeric_limits<unsigned>::max();
  std::vector<unsigned> marker(faceCount(), maxUnsigned);
  std::vector<SimulationDomain2::Orientation> norientation = {
      SimulationDomain2::Orientation::LEFT,
      SimulationDomain2::Orientation::RIGHT,
      SimulationDomain2::Orientation::BOTTOM,
      SimulationDomain2::Orientation::TOP};
  for (unsigned i = 0; i < cellCount(); i++)
    if (_boundaryHandler->cellType(i) ==
        DomainBoundaryHandler2::MaterialType::FLUID) {
      marker[faceIdFromCellId(i, norientation[0])] = 0;
      marker[faceIdFromCellId(i, norientation[1])] = 0;
      marker[faceIdFromCellId(i, norientation[2])] = 0;
      marker[faceIdFromCellId(i, norientation[3])] = 0;
    }
  std::queue<unsigned> uWavefront, vWavefront;
  _boundaryHandler->iterateBoundaryCells(
      this, DomainBoundaryHandler2::MaterialType::FLUID, [&](size_t id) {
        auto neighbors = cellNeighborhood(id);
        for (auto i : neighbors)
          if (i.id >= 0 &&
              _boundaryHandler->cellType(i.id) ==
                  DomainBoundaryHandler2::MaterialType::AIR) {
            switch (i.orientation) {
              case SimulationDomain2::Orientation::LEFT: {
                auto fid = faceIdFromCellId(i.id, norientation[0]); // LEFT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[2]); // BOTTOM
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[3]); // TOP
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
              }
                break;
              case SimulationDomain2::Orientation::RIGHT: {
                auto fid = faceIdFromCellId(i.id, norientation[1]); // RIGHT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[2]); // BOTTOM
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[3]); // TOP
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
              }
                break;
              case SimulationDomain2::Orientation::TOP: {
                auto fid = faceIdFromCellId(i.id, norientation[0]); // LEFT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[1]); // RIGHT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[3]); // TOP
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
              }
                break;
              case SimulationDomain2::Orientation::BOTTOM: {
                auto fid = faceIdFromCellId(i.id, norientation[0]); // LEFT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[1]); // RIGHT
                marker[fid] = std::min(1u, marker[fid]);
                uWavefront.push(fid);
                fid = faceIdFromCellId(i.id, norientation[2]); // BOTTOM
                marker[fid] = std::min(1u, marker[fid]);
                vWavefront.push(fid);
              }
                break;
              default: ASSERT_FATAL(false);
            }
          }
      });
  // EXTRAPOLATE U
  while (!uWavefront.empty()) {
    unsigned id = uWavefront.front();
    uWavefront.pop();
    double sum = 0.;
    int k = 0;
    for (auto o : norientation) {
      int nid = faceIdFromFaceId(id, o);
      if (nid < 0)
        continue;
      if (marker[nid] < marker[id]) {
        sum += scalarAtFace(fieldId, nid);
        k++;
      } else if (marker[nid] == maxUnsigned) {
        marker[nid] = marker[id] + 1;
        uWavefront.push(nid);
      }
    }
    if (k > 0)
      scalarAtFace(fieldId, id) = sum / k;
  }
  // EXTRAPOLATE V
  while (!vWavefront.empty()) {
    unsigned id = vWavefront.front();
    vWavefront.pop();
    double sum = 0.;
    int k = 0;
    for (auto o : norientation) {
      int nid = faceIdFromFaceId(id, o);
      if (nid < 0)
        continue;
      if (marker[nid] < marker[id]) {
        sum += scalarAtFace(fieldId, nid);
        k++;
      } else if (marker[nid] == maxUnsigned) {
        marker[nid] = marker[id] + 1;
        vWavefront.push(nid);
      }
    }
    if (k > 0)
      scalarAtFace(fieldId, id) = sum / k;
  }
}

Point2d SimRegularDomain2::cellCenteredGridCell(Point2d p, unsigned *v,
                                                Point2d *ps) const {
  auto s = _grid.gridSize();
  int n = s.x();
  int m = s.y();
  Point2d gp = p;
  Point2i lower(floor2Int(gp.x()), floor2Int(gp.y()));
  Point2i upper(ceil2Int(gp.x()), ceil2Int(gp.y()));

  Point2i safeLower(std::max(0, std::min(n - 1, lower.x())),
                    std::max(0, std::min(m - 1, lower.y())));
  Point2i safeUpper(std::max(0, std::min(n - 1, upper.x())),
                    std::max(0, std::min(m - 1, upper.y())));
  gp -= Vector2d(lower.x(), lower.y());
  std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                 Point2d(0., 1.), Point2d(1., 1.)};
  ASSERT_FATAL(safeLower.x() >= -1);
  ASSERT_FATAL(safeUpper.x() <= n);
  ASSERT_FATAL(safeLower.y() >= -1);
  ASSERT_FATAL(safeUpper.y() <= m);
  v[0] = cellIdFromIJ(safeLower);
  v[1] = cellIdFromIJ(Point2i(safeUpper.x(), safeLower.y()));
  v[2] = cellIdFromIJ(Point2i(safeLower.x(), safeUpper.y()));
  v[3] = cellIdFromIJ(safeUpper);
  if (ps)
    for (int i = 0; i < 4; i++)
      ps[i] = cellCenterPosition(v[i]);
  return gp;
}

Point2d SimRegularDomain2::uFaceCenteredGridCell(Point2d p,
                                                 unsigned *v,
                                                 double *f,
                                                 unsigned fieldId,
                                                 Point2d *ps) const {
  auto s = _grid.gridSize();
  int n = s.x() + 1;
  int m = s.y();
  Point2d gp = p;
  Point2i lower(floor2Int(gp.x()), floor2Int(gp.y()));
  Point2i upper(ceil2Int(gp.x()), ceil2Int(gp.y()));
  Point2i safeLower(std::max(0, std::min(n - 1, lower.x())),
                    std::max(0, std::min(m - 1, lower.y())));
  Point2i safeUpper(std::max(0, std::min(n - 1, upper.x())),
                    std::max(0, std::min(m - 1, upper.y())));
  gp -= Vector2d(lower.x(), lower.y());
  std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                 Point2d(0., 1.), Point2d(1., 1.)};
  ASSERT(safeLower.x() >= -1);
  ASSERT(safeUpper.x() <= n);
  ASSERT(safeLower.y() >= -1);
  ASSERT(safeUpper.y() <= m);
  v[0] = uFaceIdFromIJ(safeLower);
  v[1] = uFaceIdFromIJ(Point2i(safeUpper.x(), safeLower.y()));
  v[2] = uFaceIdFromIJ(Point2i(safeLower.x(), safeUpper.y()));
  v[3] = uFaceIdFromIJ(safeUpper);

  if (f) {
    f[0] = _uFaceFields[fieldId][v[0]];
    f[1] = _uFaceFields[fieldId][v[1]];
    f[2] = _uFaceFields[fieldId][v[2]];
    f[3] = _uFaceFields[fieldId][v[3]];
  }
  if (ps)
    for (int i = 0; i < 4; i++)
      ps[i] = faceCenterPosition(v[i]);
  return gp;
}

Point2d SimRegularDomain2::vFaceCenteredGridCell(Point2d p,
                                                 unsigned *v,
                                                 double *f,
                                                 unsigned fieldId,
                                                 Point2d *ps) const {
  auto s = _grid.gridSize();
  int n = static_cast<int>(s.x());
  int m = static_cast<int>(s.y() + 1);
  Point2d gp = p;
  Point2i lower(floor2Int(gp.x()), floor2Int(gp.y()));
  Point2i upper(ceil2Int(gp.x()), ceil2Int(gp.y()));
  Point2i safeLower(std::max(0, std::min(n - 1, lower.x())),
                    std::max(0, std::min(m - 1, lower.y())));
  Point2i safeUpper(std::max(0, std::min(n - 1, upper.x())),
                    std::max(0, std::min(m - 1, upper.y())));

  gp -= Vector2d(lower.x(), lower.y());

  std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                 Point2d(0., 1.), Point2d(1., 1.)};
  ASSERT(safeLower.x() >= -1);
  ASSERT(safeUpper.x() <= n);
  ASSERT(safeLower.y() >= -1);
  ASSERT(safeUpper.y() <= m);
  v[0] = vFaceIdFromIJ(safeLower);
  v[1] = vFaceIdFromIJ(Point2i(safeUpper.x(), safeLower.y()));
  v[2] = vFaceIdFromIJ(Point2i(safeLower.x(), safeUpper.y()));
  v[3] = vFaceIdFromIJ(safeUpper);

  if (f) {
    f[0] = _vFaceFields[fieldId][v[0] - uFacesCount()];
    f[1] = _vFaceFields[fieldId][v[1] - uFacesCount()];
    f[2] = _vFaceFields[fieldId][v[2] - uFacesCount()];
    f[3] = _vFaceFields[fieldId][v[3] - uFacesCount()];
  }
  if (ps)
    for (int i = 0; i < 4; i++)
      ps[i] = faceCenterPosition(v[i]);
  return gp;
}

double SimRegularDomain2::scalarAtCell(unsigned fieldId,
                                       unsigned cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  return _cellCenteredScalarFields[fieldId][cellId];
}

double &SimRegularDomain2::scalarAtCell(unsigned fieldId, unsigned cellId) {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  return _cellCenteredScalarFields[fieldId][cellId];
}

void SimRegularDomain2::setScalarAtCell(unsigned fieldId,
                                        unsigned cellId,
                                        double value) {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());

  _cellCenteredScalarFields[fieldId][cellId] = value;
}

double SimRegularDomain2::scalarAtFace(unsigned fieldId,
                                       unsigned faceId) const {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  auto _uFacesCount = uFacesCount();
  if (faceId < _uFacesCount) {
    ASSERT_FATAL(faceId < _uFaceFields[fieldId].size());
    return _uFaceFields[fieldId][faceId];
  }
  faceId -= _uFacesCount;
  ASSERT_FATAL(faceId < _vFaceFields[fieldId].size());
  return _vFaceFields[fieldId][faceId];
}

double &SimRegularDomain2::scalarAtFace(unsigned fieldId, unsigned faceId) {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  auto _uFacesCount = uFacesCount();
  if (faceId < _uFacesCount) {
    ASSERT_FATAL(faceId < _uFaceFields[fieldId].size());
    return _uFaceFields[fieldId][faceId];
  }
  faceId -= _uFacesCount;
  ASSERT_FATAL(faceId < _vFaceFields[fieldId].size());
  return _vFaceFields[fieldId][faceId];
}

unsigned char SimRegularDomain2::byteAtCell(unsigned fieldId,
                                            unsigned cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredByteFields.size());
  ASSERT_FATAL(cellId < _cellCenteredByteFields[fieldId].size());
  return _cellCenteredByteFields[fieldId][cellId];
}

unsigned char &SimRegularDomain2::byteAtCell(unsigned fieldId,
                                             unsigned cellId) {
  ASSERT_FATAL(fieldId < _cellCenteredByteFields.size());
  ASSERT_FATAL(cellId < _cellCenteredByteFields[fieldId].size());
  return _cellCenteredByteFields[fieldId][cellId];
}

Vector2d SimRegularDomain2::sampleFaceFieldAtCellCenter(size_t fieldId,
                                                        size_t cellId) const {
  Point2d samplePoint = this->cellCenterPosition(cellId);
  return Vector2d(sampleUFaceField(fieldId, samplePoint),
                  sampleVFaceField(fieldId, samplePoint));
}

void SimRegularDomain2::setScalarAtFace(unsigned fieldId,
                                        unsigned faceId,
                                        double value) {
  ASSERT_FATAL(fieldId < _uFaceFields.size());
  auto _uFacesCount = uFacesCount();
  if (faceId < _uFacesCount) {
    ASSERT_FATAL(faceId < _uFaceFields[fieldId].size());
    _uFaceFields[fieldId][faceId] = value;
  } else {
    faceId -= _uFacesCount;
    ASSERT_FATAL(faceId < _vFaceFields[fieldId].size());
    _vFaceFields[fieldId][faceId] = value;
  }
}

} // namespace furoo
