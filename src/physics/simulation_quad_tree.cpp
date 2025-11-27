#include <blas/interpolation.h>
#include <blas/rbf.h>
#include <map>
#include <physics/simulation_quad_tree.h>
#include <set>

namespace furoo {

SimQuadTree::SimQuadTree() {}

SimQuadTree::~SimQuadTree() {}

int SimQuadTree::cellId(Point2d p) const { return _qt.intersect(p); }

void SimQuadTree::setPointSetStructures(PointSet2Interface *cs,
                                        PointSet2Interface *fsx,
                                        PointSet2Interface *fsy) {
  _cellCenters.reset(cs);
  auto leafs = _qt.leafs();
  for (auto &l : leafs)
    _cellCenters->add(l.node->region().center());
  _uFaceCenters.reset(fsx);
  _vFaceCenters.reset(fsy);
  auto edges = _qt.leafEdges();
  for (auto i : _vEdgeIndices)
    _vFaceCenters->add(edges[i].faceCenterPosition);
  for (auto i : _uEdgeIndices)
    _uFaceCenters->add(edges[i].faceCenterPosition);
}

void SimQuadTree::setBoundaryHandler(DomainBoundaryHandler2 *bh) {
  this->_boundaryHandler = bh;
  bh->resize(cellCount(), faceCount());
}

void SimQuadTree::setInterpolationMethod(InterpolationMethod interpolant) {
  this->_interpolant = interpolant;
}

Vector2d SimQuadTree::smallestCellSize() const {
  double invsize = 1. / std::pow(2., _qt.height());
  return Vector2d(_qt.domainRegion().size(0) * invsize,
                  _qt.domainRegion().size(1) * invsize);
}

BBox2D SimQuadTree::region() const { return _qt.domainRegion(); }

BBox2D SimQuadTree::cellRegion(size_t cellId) const {
  if (!(cellId < _qt.leafCount())) {
    std::cerr << "Cell region" << '\n';
    throw ("cellRegion");
  }
  auto &leafs = _qt.leafs();
  return leafs[cellId].node->region();
}

double &SimQuadTree::faceCenteredScalar(size_t fieldId, size_t faceId) {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  ASSERT_FATAL(faceId < faceCount());
  return _faceCenteredScalarFields[fieldId][faceId];
}

void SimQuadTree::iterateCells(Point2d center, double radius,
                               const std::function<void(size_t)> &f) const {
  ASSERT_FATAL(_cellCenters);
  BBox2D b(center - Vector2d(radius), center + Vector2d(radius));
  _cellCenters->search(b, f);
}

void SimQuadTree::iterateVFaces(Point2d center, double radius,
                                const std::function<void(size_t)> &f) const {
  ASSERT_FATAL(_vFaceCenters);
  BBox2D b(center - Vector2d(radius), center + Vector2d(radius));
  _vFaceCenters->search(b, [&](size_t i) {
    ASSERT_FATAL(i < _vEdgeIndices.size());
    f(_vEdgeIndices[i]);
  });
}

void SimQuadTree::iterateUFaces(Point2d center, double radius,
                                const std::function<void(size_t)> &f) const {
  ASSERT_FATAL(_uFaceCenters);
  BBox2D b(center - Vector2d(radius), center + Vector2d(radius));
  _uFaceCenters->search(b, [&](size_t i) {
    ASSERT_FATAL(i < _uEdgeIndices.size());
    f(_uEdgeIndices[i]);
  });
}

void SimQuadTree::iterateVFaces(Point2d center,
                                const std::function<void(size_t)> &f,
                                size_t neighborsStar) const {
  neighborsStar++;
  auto cid = cellId(center);
  if (cid < 0)
    return;
  ASSERT_FATAL(cid >= 0);
  // cellId ---> node distance
  std::map<size_t, size_t> visited;
  std::queue<size_t> q;
  q.push(static_cast<size_t>(cid));
  visited[static_cast<size_t>(cid)] = 0;
  std::set<size_t> facesIds;
  while (!q.empty()) {
    size_t id = q.front();
    q.pop();
    size_t level = visited[id];
    if (level >= neighborsStar)
      continue;
    auto faces = cellFaces(static_cast<unsigned int>(id));
    for (auto &face : faces)
      if (face.orientation == SimulationDomain2::Orientation::BOTTOM ||
          face.orientation == SimulationDomain2::Orientation::TOP)
        facesIds.insert(static_cast<unsigned long &&>(face.id));
    auto neighbors = cellNeighborhood(static_cast<size_t>(id));
    for (auto &n : neighbors) {
      if (n.id < 0 || visited.find(static_cast<const unsigned long &>(n.id)) !=
          visited.end())
        continue;
      visited[n.id] = level + 1;
      q.push(n.id);
    }
  }
  for (auto &i : facesIds)
    f(i);
}

void SimQuadTree::iterateUFaces(Point2d center,
                                const std::function<void(size_t)> &f,
                                size_t neighborsStar) const {
  neighborsStar++;
  auto cid = cellId(center);
  if (cid < 0)
    return;
  ASSERT_FATAL(cid >= 0);
  // cellId ---> node distance
  std::map<size_t, size_t> visited;
  std::queue<size_t> q;
  q.push(static_cast<size_t>(cid));
  visited[static_cast<size_t>(cid)] = 0;
  std::set<size_t> facesIds;
  while (!q.empty()) {
    size_t id = q.front();
    q.pop();
    size_t level = visited[id];
    if (level >= neighborsStar)
      continue;
    auto faces = cellFaces(static_cast<unsigned int>(id));
    for (auto &face : faces)
      if (face.orientation == SimulationDomain2::Orientation::LEFT ||
          face.orientation == SimulationDomain2::Orientation::RIGHT)
        facesIds.insert(static_cast<unsigned long &&>(face.id));
    auto neighbors = cellNeighborhood(static_cast<size_t>(id));
    for (auto &n : neighbors) {
      if (n.id < 0 || visited.find(static_cast<const unsigned long &>(n.id)) !=
          visited.end())
        continue;
      visited[n.id] = level + 1;
      q.push(n.id);
    }
  }
  for (auto &i : facesIds)
    f(i);
  /*
  auto neighbors = cellNeighborhood(static_cast<size_t>(cid));
  for (auto &n : neighbors) {
    if (n.id < 0)
      continue;
    auto faces = cellFaces(static_cast<unsigned int>(n.id));
    for (auto &face: faces)
      if (face.orientation == SimulationDomain2::Orientation::RIGHT ||
          face.orientation == SimulationDomain2::Orientation::LEFT)
        f(face.id);
  }*/
}

Point2d SimQuadTree::cellCenterPosition(size_t cellId) const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(cellId < leafs.size());
  return leafs[cellId].node->region().center();
}

Point2d SimQuadTree::faceCenterPosition(size_t faceId) const {
  auto &edges = _qt.leafEdges();
  ASSERT_FATAL(faceId < edges.size());
  auto &edge = edges[faceId];
  return edge.faceCenterPosition;
}

void SimQuadTree::iterateStencilAtCellCenter(
    size_t id,
    const std::function<void(int nid, double w, Point2d p,
                             SimulationDomain2::Orientation o)> &f) const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(id < leafs.size());
  auto &leaf = leafs[id];
  auto node = *leaf.node;
  auto center = node.region().center();
  std::vector<Point2d> points(1, center);
  // compute shape parameter
  double shapeParameter = 0.;
  for (size_t e = 0; e < leaf.edges.size(); e++) {
    auto nc = neighborCellCenter(id, e);
    shapeParameter = std::max(shapeParameter, distance(nc, center));
    points.emplace_back(nc);
  }
  RBF<Point2d> rbf(
      new GaussianKernel<Point2d>(std::max(0., 1. / sqrtf(shapeParameter))));
  auto x = rbf.laplacian(points);
  auto &edges = _qt.leafEdges();
  // count neumann boundary faces
  for (size_t e = 0; e < leaf.edges.size(); e++)
    if (this->_boundaryHandler->faceBoundaryType(leaf.edges[e]) ==
        DomainBoundaryHandler2::BoundaryType::NEUMANN)
      x[0] += x[e];
  f(id, x[0], points[0], SimulationDomain2::Orientation::OTHER);
  for (size_t e = 0; e < leaf.edges.size(); e++) {
    SimulationDomain2::Orientation o;
    switch (edges[leaf.edges[e]].neighborOrientation(id)) {
      case QuadTree::Edge::NeighborOrientation::LEFT:
        o = SimulationDomain2::Orientation::LEFT;
        break;
      case QuadTree::Edge::NeighborOrientation::RIGHT:
        o = SimulationDomain2::Orientation::RIGHT;
        break;
      case QuadTree::Edge::NeighborOrientation::TOP:
        o = SimulationDomain2::Orientation::TOP;
        break;
      case QuadTree::Edge::NeighborOrientation::BOTTOM:
        o = SimulationDomain2::Orientation::BOTTOM;
        break;
    }
    f(edges[leaf.edges[e]].neighborIndex(id), x[e + 1], points[e + 1], o);
  }
}

void SimQuadTree::iterateStencilAtCellCenter(
    size_t id,
    const std::function<void(int nid, double w,
                             DomainBoundaryHandler2::BoundaryType, double)> &f)
const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(id < leafs.size());
  auto &leaf = leafs[id];
  auto node = *leaf.node;
  auto center = node.region().center();
  std::vector<Point2d> points(1, center);
  // compute shape parameter
  double shapeParameter = 0.;
  for (size_t e = 0; e < leaf.edges.size(); e++) {
    auto nc = neighborCellCenter(id, e);
    shapeParameter = std::max(shapeParameter, distance(nc, center));
    points.emplace_back(nc);
  }
  RBF<Point2d> rbf(new GaussianKernel<Point2d>(
      0.13)); // std::max(0., 1. / sqrtf(shapeParameter))));
  auto x = rbf.laplacian(points);
  auto &edges = _qt.leafEdges();
  // count neumann boundary faces
  for (size_t e = 0; e < leaf.edges.size(); e++)
    if (this->_boundaryHandler->faceBoundaryType(leaf.edges[e]) ==
        DomainBoundaryHandler2::BoundaryType::NEUMANN)
      x[0] += x[e];
  f(id, x[0], DomainBoundaryHandler2::BoundaryType::NONE, 0.);
  for (size_t e = 0; e < leaf.edges.size(); e++) {
    f(edges[leaf.edges[e]].neighborIndex(id), x[e + 1],
      this->_boundaryHandler->faceBoundaryType(leaf.edges[e]),
      this->_boundaryHandler->faceValue(leaf.edges[e]));
  }
}

std::vector<SimulationDomain2::NeighborOrientation>
SimQuadTree::faceNeighborhood(size_t faceId) const {
  auto &edges = _qt.leafEdges();
  ASSERT_FATAL(faceId < edges.size());
  auto &edge = edges[faceId];
  std::vector<SimulationDomain2::NeighborOrientation> r;
  for (size_t i = 0; i < 2; i++)
    switch (edge.np[i]) {
      case QuadTree::Edge::NeighborOrientation::LEFT:
        r.emplace_back(edge.v[i],
                       SimulationDomain2::Orientation::LEFT);
        break;
      case QuadTree::Edge::NeighborOrientation::RIGHT:
        r.emplace_back(edge.v[i],
                       SimulationDomain2::Orientation::RIGHT);
        break;
      case QuadTree::Edge::NeighborOrientation::TOP:
        r.emplace_back(edge.v[i],
                       SimulationDomain2::Orientation::TOP);
        break;
      case QuadTree::Edge::NeighborOrientation::BOTTOM:
        r.emplace_back(edge.v[i],
                       SimulationDomain2::Orientation::BOTTOM);
        break;
    }
  return r;
}

std::vector<SimulationDomain2::NeighborOrientation>
SimQuadTree::cellNeighborhood(size_t cellId) const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(cellId < leafs.size());
  auto &edges = _qt.leafEdges();
  std::vector<SimulationDomain2::NeighborOrientation> r;
  for (size_t e = 0; e < leafs[cellId].edges.size(); e++) {
    auto &edge = edges[leafs[cellId].edges[e]];
    switch (edge.neighborOrientation(cellId)) {
      case QuadTree::Edge::NeighborOrientation::LEFT:
        r.emplace_back(edge.neighborIndex(cellId),
                       SimulationDomain2::Orientation::LEFT);
        break;
      case QuadTree::Edge::NeighborOrientation::RIGHT:
        r.emplace_back(edge.neighborIndex(cellId),
                       SimulationDomain2::Orientation::RIGHT);
        break;
      case QuadTree::Edge::NeighborOrientation::TOP:
        r.emplace_back(edge.neighborIndex(cellId),
                       SimulationDomain2::Orientation::TOP);
        break;
      case QuadTree::Edge::NeighborOrientation::BOTTOM:
        r.emplace_back(edge.neighborIndex(cellId),
                       SimulationDomain2::Orientation::BOTTOM);
        break;
    }
  }
  return r;
}

void SimQuadTree::iterateDxScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  auto &edges = _qt.leafEdges();
  for (auto &i : _uEdgeIndices)
    f(_faceCenteredScalarFields[fieldId][i], edges[i].faceCenterPosition);
}

void SimQuadTree::iterateDxScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  for (auto i : _uEdgeIndices)
    f(i, _faceCenteredScalarFields[fieldId][i]);
}

void SimQuadTree::iterateDyScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  auto &edges = _qt.leafEdges();
  for (auto &i : _vEdgeIndices)
    f(_faceCenteredScalarFields[fieldId][i], edges[i].faceCenterPosition);
}

void SimQuadTree::iterateDyScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  for (auto i : _vEdgeIndices)
    f(i, _faceCenteredScalarFields[fieldId][i]);
}

void SimQuadTree::iterateCellScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(_cellCenteredScalarFields[fieldId].size() == leafs.size());
  size_t size = leafs.size();
  for (size_t i = 0; i < size; i++)
    f(_cellCenteredScalarFields[fieldId][i], leafs[i].node->region().center());
}

void SimQuadTree::iterateCellScalar(
    size_t fieldId, const std::function<void(size_t, double &)> &f) {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(_cellCenteredScalarFields[fieldId].size() == leafs.size());
  size_t size = leafs.size();
  for (size_t i = 0; i < size; i++)
    f(i, _cellCenteredScalarFields[fieldId][i]);
}

void SimQuadTree::iterateCellByte(
    size_t fieldId, const std::function<void(unsigned char &)> &f) {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(f);
}

void SimQuadTree::iterateVertexScalar(
    size_t fieldId, const std::function<void(double &, Point2d)> &f) {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(f);
}

size_t SimQuadTree::addFaceCenteredScalarField(double v) {
  _faceCenteredScalarFields.emplace_back(faceCount(), v);
  return _faceCenteredScalarFields.size() - 1;
}

size_t SimQuadTree::addCellCenteredScalarField(double v) {
  _cellCenteredScalarFields.emplace_back(_qt.leafCount(), v);
  return _cellCenteredScalarFields.size() - 1;
}

size_t SimQuadTree::addCellCenteredByteField(unsigned char v) {
  _cellCenteredByteFields.emplace_back(_qt.leafCount(), v);
  return _cellCenteredByteFields.size() - 1;
}

size_t SimQuadTree::addVertexCenteredScalarField(double v) {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(v);
  return 0;
}

void SimQuadTree::updateEdgeIndices() {
  _uEdgeIndices.clear();
  _vEdgeIndices.clear();
  auto &edges = _qt.leafEdges();
  for (size_t i = 0; i < edges.size(); i++) {
    if (edges[i].np[0] == QuadTree::Edge::NeighborOrientation::LEFT ||
        edges[i].np[0] == QuadTree::Edge::NeighborOrientation::RIGHT)
      _uEdgeIndices.emplace_back(i);
    else
      _vEdgeIndices.emplace_back(i);
  }
}

Point2d SimQuadTree::neighborCellCenter(size_t i, size_t e) const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(i < leafs.size());
  auto &edges = _qt.leafEdges();
  ASSERT_FATAL(e < leafs[i].edges.size());
  auto &leaf = leafs[i];
  auto &edge = edges[leaf.edges[e]];
  if (edge.neighborIndex(i) < 0) {
    bool isD = this->_boundaryHandler->faceBoundaryType(leaf.edges[e]) ==
        DomainBoundaryHandler2::BoundaryType::DIRICHLET;
    switch (edge.neighborOrientation(i)) {
      case QuadTree::Edge::NeighborOrientation::LEFT:
        // std::cout << "left\n";
        return leaf.node->region().center() +
            Vector2d(-leaf.node->region().size(0), 0.f) * ((isD) ? 0.5f : 1.f);
      case QuadTree::Edge::NeighborOrientation::RIGHT:
        // std::cout << "right\n";
        return leaf.node->region().center() +
            Vector2d(leaf.node->region().size(0), 0.f) * ((isD) ? 0.5f : 1.f);
      case QuadTree::Edge::NeighborOrientation::TOP:
        // std::cout << "top\n";
        return leaf.node->region().center() +
            Vector2d(0.f, leaf.node->region().size(1)) * ((isD) ? 0.5f : 1.f);
      case QuadTree::Edge::NeighborOrientation::BOTTOM:
        // std::cout << "bottom\n";
        return leaf.node->region().center() +
            Vector2d(0.f, -leaf.node->region().size(1)) * ((isD) ? 0.5f : 1.f);
    }
  }
  return leafs[edge.neighborIndex(i)].node->region().center();
}

double SimQuadTree::sampleCellCenteredScalar(size_t fieldId,
                                             Point2d target) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  UNUSED_VARIABLE(target);
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: NOT_IMPLEMENTED();
      ASSERT_FATAL(false);
      return 0.;
    case SimulationDomain2::InterpolationMethod::BICUBIC: NOT_IMPLEMENTED();
      return 0.;
    case SimulationDomain2::InterpolationMethod::RBF: {
      std::vector<double> f;
      std::vector<Point2d> points;
      iterateCells(target, 0.05, [&](size_t id) {
        points.emplace_back(cellCenterPosition(id));
        f.emplace_back(_cellCenteredScalarFields[fieldId][id]);
      });
      RBFInterpolant<Point2d> interpolant;
      return interpolant.interpolateAt(target, points, f);
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;
}

double SimQuadTree::gradientXAtCellCenteredScalar(size_t fieldId,
                                                  size_t cellId) const {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(cellId);
  return 0.;
}

double SimQuadTree::gradientYAtCellCenteredScalar(size_t fieldId,
                                                  size_t cellId) const {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(cellId);
  return 0.;
}

double SimQuadTree::gradientXAtFaceCenteredScalar(size_t fieldId,
                                                  size_t faceId) const {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(faceId);
  return 0.;
}

double SimQuadTree::gradientYAtFaceCenteredScalar(size_t fieldId,
                                                  size_t faceId) const {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(faceId);
  return 0.;
}

double SimQuadTree::faceDivergentAtCellCenter(unsigned int fieldId,
                                              unsigned int cellId) const {

  auto &leafs = _qt.leafs();
  ASSERT_FATAL(cellId < leafs.size());
  double divergence = 0.;
  auto &edges = _qt.leafEdges();
  //
  RBF<Point2d> interpolant;
  std::vector<double> f;
  std::vector<Point2d> points;
  std::set<size_t> horizontalFaces, verticalFaces;
  std::set<size_t> neighborCells;

  // Get all cell neighbors
  auto neighbors = cellNeighborhood(cellId);
  neighborCells.insert(cellId);
  for (auto n : neighbors) {
    if (n.id >= 0)
      neighborCells.insert(static_cast<unsigned long &&>(n.id));
  }

  // From neighbors, add all their faces to correspondent set
  for (auto cell : neighborCells) {
    auto &leaf = leafs[cell];
    for (size_t e = 0; e < leaf.edges.size(); e++) {
      auto face = leaf.edges[e];
      switch (edges[face].neighborOrientation(cell)) {
        case QuadTree::Edge::NeighborOrientation::LEFT:
        case QuadTree::Edge::NeighborOrientation::RIGHT:
          horizontalFaces.insert(face);
          break;
        case QuadTree::Edge::NeighborOrientation::TOP:
        case QuadTree::Edge::NeighborOrientation::BOTTOM:
          verticalFaces.insert(face);
          break;
      }
    }
  }

  // Prepare data for interpolation
  // Horizontal faces
  for (auto face : horizontalFaces) {
    points.emplace_back(faceCenterPosition(face));
    f.emplace_back(scalarAtFace(fieldId, face));
  }
  try {
    auto target = cellCenterPosition(cellId);
    auto weights = interpolant.gradientXAt(target, points);
    for (size_t i = 0; i < weights.size(); i++) {
      divergence += weights[i] * f[i];
    }
  } catch (const char *e) {
    std::cerr << "Thrown from " << e << '\n';
    throw ("faceDivergentAtCellCenter");
  }

  // Vertical faces
  points.clear();
  f.clear();
  for (auto face : verticalFaces) {
    points.emplace_back(faceCenterPosition(face));
    f.emplace_back(scalarAtFace(fieldId, face));
  }
  try {
    auto target = cellCenterPosition(cellId);
    auto weights = interpolant.gradientYAt(target, points);
    for (size_t i = 0; i < weights.size(); i++) {
      divergence += weights[i] * f[i];
    }
  } catch (const char *e) {
    std::cerr << "Thrown from " << e << '\n';
    throw ("faceDivergentAtCellCenter");
  }

  return divergence;
  // ============================
  // FINITE DIFFERENCE DIVERGENCE
  // ============================
  // auto &leaf = leafs[cellId];
  // for (size_t e = 0; e < leaf.edges.size(); e++) {
  //   switch (edges[leaf.edges[e]].neighborOrientation(cellId)) {
  //   case QuadTree::Edge::NeighborOrientation::LEFT:
  //     divergence -= scalarAtFace(fieldId, leaf.edges[e]);
  //     break;
  //   case QuadTree::Edge::NeighborOrientation::RIGHT:
  //     divergence += scalarAtFace(fieldId, leaf.edges[e]);
  //     break;
  //   case QuadTree::Edge::NeighborOrientation::TOP:
  //     divergence += scalarAtFace(fieldId, leaf.edges[e]);
  //     break;
  //   case QuadTree::Edge::NeighborOrientation::BOTTOM:
  //     divergence -= scalarAtFace(fieldId, leaf.edges[e]);
  //     break;
  //   }
  // }
  // return divergence / cellRegion(cellId).size(0);
}

double SimQuadTree::cellGradientAtFaceCenter(unsigned int fieldId,
                                             unsigned int faceId) const {
  double gradient = 0.;
  auto &edges = _qt.leafEdges();

  RBF<Point2d> interpolant;
  std::vector<double> f;
  std::vector<Point2d> points;
  auto target = edges[faceId].faceCenterPosition;
  std::set<size_t> neighborCells;
  if (edges[faceId].v[0] >= 0) {
    auto neighbors = cellNeighborhood(static_cast<size_t>(edges[faceId].v[0]));
    for (auto n : neighbors)
      if (n.id >= 0)
        neighborCells.insert(static_cast<unsigned long &&>(n.id));
  }
  if (edges[faceId].v[1] >= 0) {
    auto neighbors = cellNeighborhood(static_cast<size_t>(edges[faceId].v[1]));
    for (auto n : neighbors)
      if (n.id >= 0)
        neighborCells.insert(static_cast<unsigned long &&>(n.id));
  }
  for (auto cell : neighborCells) {
    points.emplace_back(cellCenterPosition(cell));
    f.emplace_back(_cellCenteredScalarFields[fieldId][cell]);
  }

  auto orientation = edges[faceId].neighborOrientation(edges[faceId].v[0]);

  if (orientation == QuadTree::Edge::NeighborOrientation::LEFT ||
      orientation == QuadTree::Edge::NeighborOrientation::RIGHT) {
    try {
      auto weights = interpolant.gradientXAt(target, points);
      for (size_t i = 0; i < weights.size(); i++) {
        gradient += weights[i] * f[i];
      }
    } catch (const char *e) {
      std::cerr << "Thrown from " << e << '\n';
      throw ("cellGradientAtFaceCenter::Horizontal Gradient Interpolation");
    }
  } else {
    try {
      auto weights = interpolant.gradientYAt(target, points);
      for (size_t i = 0; i < weights.size(); i++) {
        gradient += weights[i] * f[i];
      }
    } catch (const char *e) {
      std::cerr << "Thrown from " << e << '\n';
      throw ("cellGradientAtFaceCenter::Vertical Gradient Interpolation");
    }
  }

  return gradient;
  /*
  iterateCells(target, 0.005, [&](size_t id) {
    points.emplace_back(cellCenterPosition(id));
    f.emplace_back(_cellCenteredScalarFields[fieldId][id]);
    std::cout << f[f.size() - 1] << " ";
  })
  std::cout << std::endl;
  auto &edges = _qt.leafEdges();*/
  // get any cell

  // ==========================
  // FINITE DIFFERENCE GRADIENT
  // ==========================
  // int cellA = edges[faceId].v[0];
  // int cellB = edges[faceId].v[1];
  // ASSERT_FATAL(cellA >= 0 && cellB >= 0);
  //
  // int cellId = cellA; // >= 0 ? edges[faceId].v[0] : edges[faceId].v[1];
  // auto o = edges[faceId].neighborOrientation(cellId);
  // double h = distance(cellCenterPosition(cellA), cellCenterPosition(cellB));
  //
  // // Horizontal gradient
  // if (o == QuadTree::Edge::NeighborOrientation::LEFT ||
  //     o == QuadTree::Edge::NeighborOrientation::RIGHT) {
  //   if (o == QuadTree::Edge::NeighborOrientation::RIGHT)
  //     gradient = _cellCenteredScalarFields[fieldId][cellB] -
  //                _cellCenteredScalarFields[fieldId][cellA];
  //   else
  //     gradient = _cellCenteredScalarFields[fieldId][cellA] -
  //                _cellCenteredScalarFields[fieldId][cellB];
  // } else {
  //   // Vertical gradient
  //   if (o == QuadTree::Edge::NeighborOrientation::TOP)
  //     gradient = _cellCenteredScalarFields[fieldId][cellB] -
  //                _cellCenteredScalarFields[fieldId][cellA];
  //   else
  //     gradient = _cellCenteredScalarFields[fieldId][cellA] -
  //                _cellCenteredScalarFields[fieldId][cellB];
  // }
  // gradient /= h;
  //
  // return gradient;
}

Vector2d SimQuadTree::sampleFaceCenteredScalar(size_t fieldId,
                                               Point2d target) const {
  return Vector2d(sampleUFaceField(fieldId, target),
                  sampleVFaceField(fieldId, target));
}

double SimQuadTree::sampleVFaceField(size_t fieldId, Point2d target) const {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  UNUSED_VARIABLE(target);
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: {
      std::vector<double> f(4, 0.0);
      // find which cell the particle is in
      auto cell = cellId(target);
      // project it to cell coordinates
      auto p = cellRegion(cell).pointCoordinates(target);
      // find out which deslocated cell it belongs
      auto neighbors = cellNeighborhood(cell);
      {
        auto faces = cellFaces(cell);
        for (auto face : faces)
          if (face.orientation == SimulationDomain2::Orientation::BOTTOM)
            f[0] = f[1] = scalarAtFace(fieldId, face.id);
          else if (face.orientation == SimulationDomain2::Orientation::TOP)
            f[2] = f[3] = scalarAtFace(fieldId, face.id);
      }
      if (p.x() < 0.5) { // left
        p += Vector2d(0.5, 0.);
        for (auto n : neighbors)
          if (n.orientation == SimulationDomain2::Orientation::LEFT &&
              n.id >= 0) {
            // get horizontal faces information
            auto faces = cellFaces(n.id);
            for (auto face : faces)
              if (face.orientation == SimulationDomain2::Orientation::BOTTOM)
                f[0] = scalarAtFace(fieldId, face.id);
              else if (face.orientation == SimulationDomain2::Orientation::TOP)
                f[2] = scalarAtFace(fieldId, face.id);
          }
      } else { // right
        p -= Vector2d(0.5, 0.);
        for (auto n : neighbors)
          if (n.orientation == SimulationDomain2::Orientation::RIGHT &&
              n.id >= 0) {
            // get horizontal faces information
            auto faces = cellFaces(n.id);
            for (auto face : faces)
              if (face.orientation == SimulationDomain2::Orientation::BOTTOM)
                f[1] = scalarAtFace(fieldId, face.id);
              else if (face.orientation == SimulationDomain2::Orientation::TOP)
                f[3] = scalarAtFace(fieldId, face.id);
          }
      }
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return BilinearInterpolant().interpolateAt(p, points, f);
    }
    case SimulationDomain2::InterpolationMethod::BICUBIC: NOT_IMPLEMENTED();
      return 0.;
    case SimulationDomain2::InterpolationMethod::RBF: {
      // IF the cell is in domain boundary, use bilinear
      // find which cell the particle is in
      auto cell = cellId(target);
      // check if boundary
      auto neighbors = cellNeighborhood(cell);
      // LEFT RIGHT
      int neighborsIds[2] = {-1, -1};
      for (auto n : neighbors)
        if (n.orientation == SimulationDomain2::Orientation::LEFT)
          neighborsIds[0] = n.id;
        else if (n.orientation == SimulationDomain2::Orientation::RIGHT)
          neighborsIds[1] = n.id;
      // project it to cell coordinates
      auto p = cellRegion(cell).pointCoordinates(target);
      if ((p.x() < 0.5 && neighborsIds[0] < 0) ||
          (p.x() >= 0.5 && neighborsIds[1] < 0)) {
        auto faces = cellFaces(cell);
        std::vector<double> f(4, 0.0);
        for (auto face : faces)
          if (face.orientation == SimulationDomain2::Orientation::BOTTOM)
            f[0] = f[1] = scalarAtFace(fieldId, face.id);
          else if (face.orientation == SimulationDomain2::Orientation::TOP)
            f[2] = f[3] = scalarAtFace(fieldId, face.id);
        std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                       Point2d(0., 1.), Point2d(1., 1.)};
        return BilinearInterpolant().interpolateAt(p, points, f);
        // for (auto face : cellFaces(cell))
        //   if (face.orientation == SimulationDomain2::Orientation::TOP
        //       || face.orientation == SimulationDomain2::Orientation::BOTTOM) {
        //     points.emplace_back(faceCenterPosition(face.id) - Vector2d(h, 0));
        //     f.emplace_back(scalarAtFace(fieldId, face.id));
        //     std::cerr << "REFLECTING FACE " << face.id << " " <<
        //     scalarAtFace(fieldId, face.id)
        //               << faceCenterPosition(face.id) << std::endl;
        //   }
      }
      // else if (p.x() >= 0.5 && neighborsIds[1] < 0) {
      //   for (auto face : faces)
      //     if (face.orientation == SimulationDomain2::Orientation::BOTTOM)
      //       f[0] = f[1] = scalarAtFace(fieldId, face.id);
      //     else if (face.orientation == SimulationDomain2::Orientation::TOP)
      //       f[2] = f[3] = scalarAtFace(fieldId, face.id);
      // }

      // NOT in boundary. Using RBF
      RBFInterpolant<Point2d> interpolant{};
      std::vector<double> f;
      std::vector<Point2d> points;
      std::vector<unsigned> fluidFaces;
      std::set<unsigned> fluidFacesIds;
      iterateVFaces(target, /*0.01,*/
                    [&](size_t id) {
                      if (fluidFacesIds.find(id) != fluidFacesIds.end())
                        return;
                      double
                          scalarValue = _faceCenteredScalarFields[fieldId][id];
                      // std::cerr << "ID " << id << '\n';
                      // std::cerr << "pos " << faceCenterPosition(id) << '\n';
                      // std::cerr << "scalar " << scalarValue << '\n';
                      // if (_boundaryHandler->faceType(id) ==
                      //     DomainBoundaryHandler2::MaterialType::AIR) {
                      //   // DomainBoundaryHandler2::MaterialType::FLUID) {
                      //   std::vector<double> values;
                      //   std::vector<Point2d> fluidPoints;
                      //   _vFaceCenters->iterateClosestPoints(
                      //       faceCenterPosition(id), 2,
                      //       [&](size_t i, Point2d p) {
                      //         values.emplace_back(
                      //             _faceCenteredScalarFields[fieldId][_vEdgeIndices[i]]);
                      //         fluidPoints.emplace_back(p);
                      //       },
                      //       [&](size_t i) -> bool {
                      //         return
                      //         _boundaryHandler->faceType(_vEdgeIndices[i]) ==
                      //                DomainBoundaryHandler2::MaterialType::FLUID;
                      //       });
                      //   try {
                      //     if ((target - fluidPoints[0]).length() ==
                      //         (target - fluidPoints[1]).length())
                      //       scalarValue = (values[0] + values[1]) / 2;
                      //     else
                      //       scalarValue = values[0];
                      //     // ASSERT_FATAL(values.size() == 1);
                      //     // scalarValue = values[0];
                      //     // scalarValue =
                      //     interpolant.interpolateAt(faceCenterPosition(id),
                      //     // fluidPoints, values);
                      //   } catch (const char *e) {
                      //     for (auto p : fluidPoints)
                      //       std::cerr << p << " ";
                      //   }
                      // } else if (_boundaryHandler->faceType(id) ==
                      //            DomainBoundaryHandler2::MaterialType::SOLID) {
                      //   scalarValue = 0.;
                      // }
                      fluidFacesIds.insert(id);
                      fluidFaces.emplace_back(f.size());
                      points.emplace_back(faceCenterPosition(id));
                      f.emplace_back(scalarValue);
                    },
                    2);

      // // Copy values to the outside of domain
      // auto cell = cellId(target);
      // auto neighbors = cellNeighborhood(cell);
      // auto faces = cellFaces(cell);
      // auto h = cellRegion(cell).size(0);
      // // We are in V faces
      // // Check if either left or right are solid
      // for (auto face : faces) {
      //   if (face.orientation == SimulationDomain2::Orientation::LEFT &&
      //       _boundaryHandler->faceType(face.id) ==
      //           DomainBoundaryHandler2::MaterialType::SOLID) {
      //     // If left is solid:
      //     for (auto face_ : faces) {
      //       // Copy values to the interior of the domain
      //       if (face_.orientation == SimulationDomain2::Orientation::TOP ||
      //           face_.orientation == SimulationDomain2::Orientation::BOTTOM) {
      //         auto fPosition = faceCenterPosition(face_.id);
      //         points.emplace_back(Point2d(fPosition.x() - h, fPosition.y()));
      //         // Check if corner cell
      //         if (_boundaryHandler->faceType(face_.id) ==
      //             DomainBoundaryHandler2::MaterialType::SOLID) {
      //           f.emplace_back(0.0);
      //         } else {
      //           f.emplace_back(scalarAtFace(fieldId, face_.id));
      //         }
      //       }
      //     }
      //   } else if (face.orientation == SimulationDomain2::Orientation::RIGHT &&
      //              _boundaryHandler->faceType(face.id) ==
      //                  DomainBoundaryHandler2::MaterialType::SOLID) {
      //     // If right face is solid
      //     for (auto face_ : faces) {
      //       if (face_.orientation == SimulationDomain2::Orientation::TOP ||
      //           face_.orientation == SimulationDomain2::Orientation::BOTTOM) {
      //         auto fPosition = faceCenterPosition(face_.id);
      //         points.emplace_back(Point2d(fPosition.x() + h, fPosition.y()));
      //         if (_boundaryHandler->faceType(face_.id) ==
      //             DomainBoundaryHandler2::MaterialType::SOLID)
      //           f.emplace_back(0.0);
      //         else
      //           f.emplace_back(scalarAtFace(fieldId, face_.id));
      //       }
      //     }
      //   }
      // }
      // ShepardInterpolant2 interpolant;
      try {
        return interpolant.interpolateAt(target, points, f);
      } catch (const char *e) {
        std::cerr << "Thrown from " << e << '\n';
        throw ("sampleVFaceField: InterpolationMethod::RBF");
      }
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;
}

double SimQuadTree::sampleUFaceField(size_t fieldId, Point2d target) const {
  ASSERT_FATAL(fieldId < _faceCenteredScalarFields.size());
  UNUSED_VARIABLE(target);
  switch (this->_interpolant) {
    case SimulationDomain2::InterpolationMethod::BILINEAR: {
      std::vector<double> f(4, 0.0);
      // find which cell the particle is in
      auto cell = cellId(target);
      // project it to cell coordinates
      auto p = cellRegion(cell).pointCoordinates(target);
      // find out which deslocated cell it belongs
      auto neighbors = cellNeighborhood(cell);
      {
        auto faces = cellFaces(cell);
        for (auto face : faces)
          if (face.orientation == SimulationDomain2::Orientation::LEFT) {
            // std::cerr << "face " << face.id << '\n';
            // std::cerr << faceCenterPosition(face.id) << '\n';
            f[0] = f[2] = scalarAtFace(fieldId, face.id);
          } else if (face.orientation
              == SimulationDomain2::Orientation::RIGHT) {
            f[1] = f[3] = scalarAtFace(fieldId, face.id);
            // std::cerr << "face " << face.id << '\n';
            // std::cerr << faceCenterPosition(face.id) << '\n';
          }
      }
      if (p.y() < 0.5) { // bottom
        p += Vector2d(0., 0.5);
        for (auto n : neighbors)
          if (n.orientation == SimulationDomain2::Orientation::BOTTOM &&
              n.id >= 0) {
            // get horizontal faces information
            auto faces = cellFaces(n.id);
            for (auto face : faces) {
              if (face.orientation == SimulationDomain2::Orientation::LEFT)
                f[0] = scalarAtFace(fieldId, face.id);
              else if (face.orientation
                  == SimulationDomain2::Orientation::RIGHT)
                f[1] = scalarAtFace(fieldId, face.id);
            }
          }
      } else { // top
        p -= Vector2d(0., 0.5);
        for (auto n : neighbors)
          if (n.orientation == SimulationDomain2::Orientation::TOP
              && n.id >= 0) {
            // get horizontal faces information
            auto faces = cellFaces(n.id);
            for (auto face : faces)
              if (face.orientation == SimulationDomain2::Orientation::LEFT)
                f[2] = scalarAtFace(fieldId, face.id);
              else if (face.orientation
                  == SimulationDomain2::Orientation::RIGHT)
                f[3] = scalarAtFace(fieldId, face.id);
          }
      }
      std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                     Point2d(0., 1.), Point2d(1., 1.)};
      return BilinearInterpolant().interpolateAt(p, points, f);
    }
    case SimulationDomain2::InterpolationMethod::BICUBIC: NOT_IMPLEMENTED();
      return 0.;
    case SimulationDomain2::InterpolationMethod::RBF: {
      // IF the cell is in domain boundary, use bilinear
      // find which cell the particle is in
      auto cell = cellId(target);
      // check if boundary
      auto neighbors = cellNeighborhood(cell);
      // LEFT RIGHT
      int neighborsIds[2] = {-1, -1};
      for (auto n : neighbors)
        if (n.orientation == SimulationDomain2::Orientation::BOTTOM)
          neighborsIds[0] = n.id;
        else if (n.orientation == SimulationDomain2::Orientation::TOP)
          neighborsIds[1] = n.id;
      // project it to cell coordinates
      auto p = cellRegion(cell).pointCoordinates(target);
      if ((p.y() < 0.5 && neighborsIds[0] < 0) ||
          (p.y() >= 0.5 && neighborsIds[1] < 0)) {
        std::vector<double> f(4, 0.0);
        auto faces = cellFaces(cell);
        for (auto face : faces)
          if (face.orientation == SimulationDomain2::Orientation::LEFT)
            f[0] = f[2] = scalarAtFace(fieldId, face.id);
          else if (face.orientation == SimulationDomain2::Orientation::RIGHT)
            f[1] = f[3] = scalarAtFace(fieldId, face.id);
        std::vector<Point2d> points = {Point2d(0., 0.), Point2d(1., 0.),
                                       Point2d(0., 1.), Point2d(1., 1.)};
        return BilinearInterpolant().interpolateAt(p, points, f);
      }

      RBFInterpolant<Point2d> interpolant{};
      std::vector<double> f;
      std::vector<Point2d> points;
      std::vector<unsigned> fluidFaces;
      std::set<unsigned> fluidFacesIds;
      iterateUFaces(target, /*0.01,*/
                    [&](size_t id) {
                      if (fluidFacesIds.find(id) != fluidFacesIds.end())
                        return;
                      double
                          scalarValue = _faceCenteredScalarFields[fieldId][id];
                      // std::cerr << "ID " << id << '\n';
                      // std::cerr << "pos " << faceCenterPosition(id) << '\n';
                      // std::cerr << "scalar " << scalarValue << '\n';
                      // if (_boundaryHandler->faceType(id) ==
                      //     DomainBoundaryHandler2::MaterialType::AIR) {
                      //   // DomainBoundaryHandler2::MaterialType::FLUID) {
                      //   std::vector<double> values;
                      //   std::vector<Point2d> fluidPoints;
                      //   _uFaceCenters->iterateClosestPoints(
                      //       faceCenterPosition(id), 2,
                      //       [&](size_t i, Point2d p) {
                      //         values.emplace_back(
                      //             _faceCenteredScalarFields[fieldId][_uEdgeIndices[i]]);
                      //         fluidPoints.emplace_back(p);
                      //       },
                      //       [&](size_t i) -> bool {
                      //         return
                      //         _boundaryHandler->faceType(_uEdgeIndices[i]) ==
                      //                DomainBoundaryHandler2::MaterialType::FLUID;
                      //       });
                      //   try {
                      //     if ((target - fluidPoints[0]).length() ==
                      //         (target - fluidPoints[1]).length())
                      //       scalarValue = (values[0] + values[1]) / 2;
                      //     else
                      //       scalarValue = values[0];
                      //     // ASSERT_FATAL(values.size() == 1);
                      //     // scalarValue =
                      //     //
                      //     interpolant.interpolateAt(faceCenterPosition(id),
                      //     // fluidPoints, values);
                      //     // scalarValue = values[0];
                      //   } catch (const char *e) {
                      //     for (auto p : fluidPoints)
                      //       std::cerr << p << " ";
                      //     std::cerr << "Thrown from  " << e << std::endl;
                      //     throw("sampleUFaceField:
                      //     InterpolationMethod::RBF");
                      //   }
                      // } else if (_boundaryHandler->faceType(id) ==
                      //            DomainBoundaryHandler2::MaterialType::SOLID)
                      //            {
                      //   scalarValue = 0.;
                      // }
                      fluidFacesIds.insert(id);
                      fluidFaces.emplace_back(f.size());
                      points.emplace_back(faceCenterPosition(id));
                      f.emplace_back(scalarValue);
                    },
                    2);

      // ShepardInterpolant2 interpolant;
      //
      // auto cell = cellId(target);
      // auto neighbors = cellNeighborhood(cell);
      // auto faces = cellFaces(cell);
      // auto h = cellRegion(cell).size(0);
      // for (auto face : faces) {
      //   if (face.orientation == SimulationDomain2::Orientation::BOTTOM &&
      //       _boundaryHandler->faceType(face.id) ==
      //           DomainBoundaryHandler2::MaterialType::SOLID) {
      //     for (auto face_ : faces) {
      //       if (face_.orientation == SimulationDomain2::Orientation::LEFT ||
      //           face_.orientation == SimulationDomain2::Orientation::RIGHT) {
      //         auto fPosition = faceCenterPosition(face_.id);
      //         points.emplace_back(Point2d(fPosition.x(), fPosition.y() - h));
      //         if (_boundaryHandler->faceType(face_.id) ==
      //             DomainBoundaryHandler2::MaterialType::SOLID)
      //           f.emplace_back(0.0);
      //         else
      //           f.emplace_back(scalarAtFace(fieldId, face_.id));
      //       }
      //     }
      //   } else if (face.orientation == SimulationDomain2::Orientation::TOP &&
      //              _boundaryHandler->faceType(face.id) ==
      //                  DomainBoundaryHandler2::MaterialType::SOLID) {
      //     for (auto face_ : faces) {
      //       if (face_.orientation == SimulationDomain2::Orientation::LEFT ||
      //           face_.orientation == SimulationDomain2::Orientation::RIGHT) {
      //         auto fPosition = faceCenterPosition(face_.id);
      //         points.emplace_back(Point2d(fPosition.x(), fPosition.y() + h));
      //         if (_boundaryHandler->faceType(face_.id) ==
      //             DomainBoundaryHandler2::MaterialType::SOLID)
      //           f.emplace_back(0.0);
      //         else
      //           f.emplace_back(scalarAtFace(fieldId, face_.id));
      //       }
      //     }
      //   }
      // }
      try {
        return interpolant.interpolateAt(target, points, f);
      } catch (const char *e) {
        std::cerr << "Thrown from " << e << '\n';
        throw ("sampleUFaceField: InterpolationMethod::RBF");
      }
    }
    default: ASSERT_FATAL(false);
  }
  ASSERT_FATAL(false);
  return 0.;
}

void SimQuadTree::sampleDxScalar(size_t fieldId,
                                 const std::vector<Point2d> &targets,
                                 LinearVector &values) const {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(targets);
  UNUSED_VARIABLE(values);
}

void SimQuadTree::sampleDyScalar(size_t fieldId,
                                 const std::vector<Point2d> &targets,
                                 LinearVector &values) const {
  UNUSED_VARIABLE(fieldId);
  UNUSED_VARIABLE(targets);
  UNUSED_VARIABLE(values);
}

double SimQuadTree::scalarAtCell(unsigned fieldId, unsigned cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  return _cellCenteredScalarFields[fieldId][cellId];
}

double &SimQuadTree::scalarAtCell(unsigned fieldId, unsigned cellId) {
  ASSERT_FATAL(fieldId < _cellCenteredScalarFields.size());
  ASSERT_FATAL(cellId < _cellCenteredScalarFields[fieldId].size());
  return _cellCenteredScalarFields[fieldId][cellId];
}

double SimQuadTree::scalarAtFace(unsigned fieldId, unsigned faceId) const {
  ASSERT_FATAL(faceId < _faceCenteredScalarFields[fieldId].size());
  return _faceCenteredScalarFields[fieldId][faceId];
}

double &SimQuadTree::scalarAtFace(unsigned fieldId, unsigned faceId) {
  ASSERT_FATAL(faceId < _faceCenteredScalarFields[fieldId].size());
  return _faceCenteredScalarFields[fieldId][faceId];
}

unsigned char SimQuadTree::byteAtCell(unsigned fieldId, unsigned cellId) const {
  ASSERT_FATAL(fieldId < _cellCenteredByteFields.size());
  ASSERT_FATAL(cellId < _cellCenteredByteFields[fieldId].size());
  return _cellCenteredByteFields[fieldId][cellId];
}

unsigned char &SimQuadTree::byteAtCell(unsigned fieldId, unsigned cellId) {
  ASSERT_FATAL(fieldId < _cellCenteredByteFields.size());
  ASSERT_FATAL(cellId < _cellCenteredByteFields[fieldId].size());
  return _cellCenteredByteFields[fieldId][cellId];
}

void SimQuadTree::extrapolateToDomain(size_t fieldId) {
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
      for (auto f : cellFaces(i))
        marker[f.id] = 0;
    }
  std::queue<unsigned> uWavefront, vWavefront;
  _boundaryHandler->iterateBoundaryCells(
      this, DomainBoundaryHandler2::MaterialType::FLUID, [&](size_t id) {
        auto neighbors = cellNeighborhood(id);
        for (auto i : neighbors)
          if (i.id >= 0 && _boundaryHandler->cellType(i.id) ==
              DomainBoundaryHandler2::MaterialType::AIR) {
            switch (i.orientation) {
              case SimulationDomain2::Orientation::LEFT: {
                for (auto f : cellFaces(i.id))
                  if (f.orientation == norientation[0]) { // LEFT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[2]) { // BOTTOM
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  } else if (f.orientation == norientation[3]) { // TOP
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  }
              }
                break;
              case SimulationDomain2::Orientation::RIGHT: {
                for (auto f : cellFaces(i.id))
                  if (f.orientation == norientation[1]) { // RIGHT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[2]) { // BOTTOM
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  } else if (f.orientation == norientation[3]) { // TOP
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  }
              }
                break;
              case SimulationDomain2::Orientation::TOP: {
                for (auto f : cellFaces(i.id))
                  if (f.orientation == norientation[0]) { // LEFT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[1]) { // RIGHT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[3]) { // TOP
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  }
              }
                break;
              case SimulationDomain2::Orientation::BOTTOM: {
                for (auto f : cellFaces(i.id))
                  if (f.orientation == norientation[0]) { // LEFT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[1]) { // RIGHT
                    marker[f.id] = std::min(1u, marker[f.id]);
                    uWavefront.push(f.id);
                  } else if (f.orientation == norientation[2]) { // BOTTOM
                    marker[f.id] = std::min(1u, marker[f.id]);
                    vWavefront.push(f.id);
                  }
              }
                break;
              default: ASSERT_FATAL(false);
            }
          }
      });
  auto faceIdFromCellId = [&](int cellId,
                              SimulationDomain2::Orientation p) -> int {
    if (cellId < 0)
      return -1;
    for (auto n : cellFaces(cellId))
      if (n.orientation == p)
        return n.id;
    return -1;
  };
  auto faceIdFromFaceId = [&](size_t faceId,
                              SimulationDomain2::Orientation o) -> int {
    auto fn = faceNeighborhood(faceId);
    auto cells = faceNeighborhood(faceId);
    for (auto cell : cells)
      if (cell.orientation == o)
        return faceIdFromCellId(cell.id, o);

    if (cells[0].orientation == SimulationDomain2::Orientation::LEFT ||
        cells[0].orientation ==
            SimulationDomain2::Orientation::RIGHT) { // it is a uface
      // get first valid cell and jump to neighbor
      for (auto cell : cells) {
        if (cell.id < 0)
          continue;
        // get left neighbor from the neighbor cell
        for (auto n : cellNeighborhood(cell.id))
          if (n.orientation == o)
            return faceIdFromCellId(
                n.id, (cell.orientation == SimulationDomain2::Orientation::LEFT)
                      ? SimulationDomain2::Orientation::RIGHT
                      : SimulationDomain2::Orientation::LEFT);
      }
    } else {
      // get first valid cell and jump to neighbor
      for (auto cell : cells) {
        if (cell.id < 0)
          continue;
        // get left neighbor from the neighbor cell
        for (auto n : cellNeighborhood(cell.id))
          if (n.orientation == o)
            return faceIdFromCellId(
                n.id, (cell.orientation == SimulationDomain2::Orientation::TOP)
                      ? SimulationDomain2::Orientation::BOTTOM
                      : SimulationDomain2::Orientation::TOP);
      }
    }
    ASSERT_FATAL(false);
    return -1;
  };
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

const QuadTree &SimQuadTree::tree() const { return _qt; }

void SimQuadTree::setTree(const QuadTree &tree) { _qt = tree; }

std::vector<SimulationDomain2::NeighborOrientation>
SimQuadTree::cellFaces(unsigned cellId) const {
  auto &leafs = _qt.leafs();
  ASSERT_FATAL(cellId < leafs.size());
  auto &leaf = leafs[cellId];
  auto &edges = _qt.leafEdges();
  std::vector<SimulationDomain2::NeighborOrientation> faces;
  for (size_t e = 0; e < leaf.edges.size(); e++) {
    switch (edges[leaf.edges[e]].neighborOrientation(cellId)) {
      case QuadTree::Edge::NeighborOrientation::LEFT:
        faces.emplace_back(leaf.edges[e],
                           SimulationDomain2::Orientation::LEFT);
        break;
      case QuadTree::Edge::NeighborOrientation::RIGHT:
        faces.emplace_back(leaf.edges[e],
                           SimulationDomain2::Orientation::RIGHT);
        break;
      case QuadTree::Edge::NeighborOrientation::TOP:
        faces.emplace_back(leaf.edges[e],
                           SimulationDomain2::Orientation::TOP);
        break;
      case QuadTree::Edge::NeighborOrientation::BOTTOM:
        faces.emplace_back(leaf.edges[e],
                           SimulationDomain2::Orientation::BOTTOM);
        break;
    }
  }
  return faces;
}

} // namespace furoo
