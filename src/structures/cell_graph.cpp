#include "cell_graph.h"

namespace furoo {

static auto defaultCellGraph2Interpolant()
{
  static CubicKernel<Point<double, 2>> kernel;
  return new RBFInterpolant<Point<double, 2>>{kernel};
}

CellGraph2::~CellGraph2() = default;

CellGraph2::CellGraph2(BBox2d region, size_t maxLevel)
    : _maxLevel(maxLevel), _domainRegion(region),
      _interpolator(defaultCellGraph2Interpolant()) {
  // create vertices
  newVertex(Point2d(region.lower().x(), region.lower().y()));
  newVertex(Point2d(region.upper().x(), region.lower().y()));
  newVertex(Point2d(region.upper().x(), region.upper().y()));
  newVertex(Point2d(region.lower().x(), region.upper().y()));
  size_t vertices[4] = {0, 1, 2, 3};
  newNode(vertices, 0);
  _nodes[0].setAddressCode(Point2u(2 << maxLevel, 2 << maxLevel));
  newNode(vertices, 0);
  _nodes[1].setAddressCode(Point2u(0, 0));
  auto regionCenter = region.center();
  // create edges and connect them to root node
  newEdge(Definitions::Side::LEFT, 0, 1,
          Point2d(region.lower().x(), regionCenter.y()));
  newEdge(Definitions::Side::RIGHT, 0, 1,
          Point2d(region.upper().x(), regionCenter.y()));
  newEdge(Definitions::Side::BOTTOM, 0, 1,
          Point2d(regionCenter.x(), region.lower().y()));
  newEdge(Definitions::Side::TOP, 0, 1,
          Point2d(regionCenter.x(), region.upper().y()));
  // connect vertices our cell
  // right now we have the following scheme
  //           e3
  //   v3______|______v2
  //    |             |
  // e0--|             |--e1
  //    |_____________|
  //   v0     |       v1
  //          e2
}

CellGraph2::CellGraph2(const BBox2d &region,
                       const std::function<bool(const CellGraph2::Node &)> &f,
                       size_t maxLevel)
    : CellGraph2(region, maxLevel) {
  refine(1, f);
}

std::vector<size_t> CellGraph2::refine(size_t nodeId) {
  THROW(nodeId < _nodes.size(), "invalid node id.");
  // This code is quite delicate. We need to track update neighbors, edges and
  // vertices! Handling nodes:
  //    - Remove the parent node and replace it by 4 children as in a quadtree
  // Handling edges:
  //    - All edges are removed, however we need to save the neighboring nodes
  //    first
  //      4 new edges are created to connect the sibling nodes that were just
  //      created,
  //    - Edges are created to connect children to the parent node neighbors,
  //      care must be taken in order connect each child to its respective
  //      neighbor
  // Handling vertices:
  //    - When refining no vertex needs to be removed
  //    - A center vertex is created and middle face vertices too (only if
  //    needed)
  // get parent node information and neighbors list
  auto nodeNeighbors = neighbors(nodeId);
  auto children = createChildren(nodeId);
  // assign address codes
  Point2u parentAddressCode = _nodes[nodeId].addressCode();
  auto parentLevel = _nodes[nodeId].level();
  size_t indexStep = 1u << (_maxLevel - (parentLevel + 1));
  _nodes[children[0]].setAddressCode(
      Point2u(parentAddressCode.x(), parentAddressCode.y() + indexStep));
  _nodes[children[1]].setAddressCode(Point2u(
      parentAddressCode.x() + indexStep, parentAddressCode.y() + indexStep));
  _nodes[children[2]].setAddressCode(
      Point2u(parentAddressCode.x(), parentAddressCode.y()));
  _nodes[children[3]].setAddressCode(
      Point2u(parentAddressCode.x() + indexStep, parentAddressCode.y()));
  // create edges between children
  interconnectChildren(children);
  // we now destroy this node to replace it by the children nodes
  destroyNode(nodeId);
  // for each neighbor, connect it to the correct child
  // add current node to avaliable Ids
  for (auto neighbor : nodeNeighbors) {
    if (_nodes[neighbor.id].level() <= parentLevel)
      connectChildrenToBigNeighbor(children, neighbor);
    else
      connectChildrenToSmallNeighbor(children, neighbor);
  }
  return children;
}

void CellGraph2::refine(size_t nodeId,
                        std::function<bool(const CellGraph2::Node &)> f) {

  std::queue<size_t> nodes;
  nodes.push(nodeId);
  while (!nodes.empty()) {
    auto curNode = nodes.front();
    nodes.pop();
    if (f(_nodes[curNode])) {
      auto children = refine(curNode);
      for (auto child : children)
        nodes.push(child);
    }
  }
}

void CellGraph2::refineAll() {
  std::vector<size_t> allIds;
  iterateCells([&](size_t cellId) { allIds.emplace_back(cellId); });
  for (size_t i = 0; i < allIds.size(); i++) {
    refine(allIds[i]);
  }
}

std::vector<Definitions::Neighbor> CellGraph2::neighbors(size_t nodeId) const {
  THROW(nodeId < _nodes.size(), "invalid node id.");
  auto edges = _nodes[nodeId].edges();
  std::vector<Definitions::Neighbor> _neighbors;
  for (auto e : edges)
    _neighbors.emplace_back(_edges[e].nodeNeighbor(nodeId));
  return _neighbors;
}

std::vector<size_t> CellGraph2::neighbors(size_t nodeId,
                                          Definitions::Side side) const {
  THROW(nodeId < _nodes.size(), "invalid node id.");
  auto edges = _nodes[nodeId].edges();
  std::vector<size_t> _neighbors;
  for (auto e : edges) {
    auto neighbor = _edges[e].nodeNeighbor(nodeId);
    if (neighbor.side == side)
      _neighbors.emplace_back(neighbor.id);
  }
  return _neighbors;
}

size_t CellGraph2::newNode(const size_t vertices[], size_t level) {
  size_t id = _nodes.size();
  BBox2d r(_vertices[vertices[0]].position(),
           _vertices[vertices[2]].position());
  if (!_availableNodeIds.empty()) {
    id = _availableNodeIds.front();
    _availableNodeIds.pop();
    _nodes[id].set(vertices, r, level);
    for (auto fieldId : this->_cellFieldIds)
      this->_fields[fieldId]->reset(id);
  } else {
    _nodes.emplace_back(vertices, r, id, level);
    for (auto fieldId : this->_cellFieldIds)
      this->_fields[fieldId]->increaseSize(1);
  }
  // increment cell counter of each vertex
  for (size_t i = 0; i < 4; i++)
    _vertices[vertices[i]].addCell();
  return id;
}

size_t CellGraph2::newVertex(const Point2d &position) {
  size_t id = _vertices.size();
  if (!_availableVertexIds.empty()) {
    id = _availableVertexIds.front();
    _availableVertexIds.pop();
    _vertices[id].set(position, id);
    for (auto fieldId : this->_vertexFieldIds)
      this->_fields[fieldId]->reset(id);
  } else {
    _vertices.emplace_back(position, id);
    for (auto fieldId : this->_vertexFieldIds)
      this->_fields[fieldId]->increaseSize(1);
  }
  return id;
}

void CellGraph2::destroyVertex(size_t vertexId) {
  THROW(vertexId < _vertices.size(),
        "CellGraph2::destroyVertex _vertexSize < vertex id")
  THROW(_vertices[vertexId].isValid(),
        "CellGraph2::destroyVertex invalid vertex.");
  _availableVertexIds.push(vertexId);
  _vertices[vertexId].destroy();
  // TODO: destroying a vertex should also destroy its nodes and edges...
}

void CellGraph2::destroyNode(size_t nodeId) {
  THROW(nodeId < _nodes.size(), "CellGraph2::destroyNode _nodesSize < node id")
  THROW(_nodes[nodeId].isValid(), "CellGraph2::destroyNode invalid node.");
  _availableNodeIds.push(nodeId);
  auto edges = _nodes[nodeId].edgesCopy();
  for (auto e : edges)
    destroyEdge(e);
  // decrement cell count of each vertex
  for (size_t i = 0; i < 4; i++) {
    _vertices[_nodes[nodeId].vertex(i)].removeCell();
    if (!_vertices[_nodes[nodeId].vertex(i)].cellCount())
      destroyVertex(_nodes[nodeId].vertex(i));
  }
  _nodes[nodeId].destroy();
}

size_t CellGraph2::newEdge(Definitions::Side aSide, size_t vertexA,
                           size_t vertexB, Point2d center) {
  // if (vertexA > 0 && vertexB > 0)
  //   center = (_nodes[vertexA].region().center() +
  //             _nodes[vertexB].region().center()) *
  //            0.5;
  size_t id = _edges.size();
  if (!_availableEdgeIds.empty()) {
    id = _availableEdgeIds.front();
    _availableEdgeIds.pop();
    _edges[id].set(aSide, vertexA, vertexB, center);
    for (auto fieldId : this->_faceFieldIds)
      this->_fields[fieldId]->reset(id);
  } else {
    _edges.emplace_back(aSide, vertexA, vertexB, center, id);
    for (auto fieldId : this->_faceFieldIds)
      this->_fields[fieldId]->increaseSize(1);
  }
  if (_nodes[vertexA].isValid())
    _nodes[vertexA].addEdge(id);
  if (_nodes[vertexB].isValid())
    _nodes[vertexB].addEdge(id);
  return id;
}

void CellGraph2::destroyEdge(size_t edgeId) {
  THROW(edgeId < _edges.size() && _edges[edgeId].isValid(),
        "invalid edge id or trying to destroy an invalid edge.");
  _availableEdgeIds.push(edgeId);
  auto nodes = _edges[edgeId].nodes();
  _nodes[nodes[0]].removeEdge(edgeId);
  _nodes[nodes[1]].removeEdge(edgeId);
  _edges[edgeId].destroy();
}

bool CellGraph2::canCoarse(std::vector<size_t> nodesId) const {
  if (nodesId.size() != 4)
    return false;
  size_t level = _nodes[nodesId[0]].level();
  size_t indexStep = 1u << (_maxLevel - (level - 1u));
  Point2u nodeFamily = _nodes[nodesId[0]].addressCode() / indexStep;
  for (auto id : nodesId) {
    if (_nodes[id].level() != level)
      return false;
    if (nodeFamily != _nodes[id].addressCode() / indexStep)
      return false;
  }
  return true;
}

std::vector<size_t> CellGraph2::nodeSiblings(size_t nodeId) const {
  size_t level = _nodes[nodeId].level();
  size_t indexStep = 1u << (_maxLevel - (level - 1u));
  Point2u nodeFamily = _nodes[nodeId].addressCode() / indexStep;
  std::vector<Definitions::Neighbor> siblings;
  for (auto n : neighbors(nodeId))
    if (level == _nodes[n.id].level() &&
        nodeFamily == _nodes[n.id].addressCode() / indexStep)
      siblings.emplace_back(n);
  if (siblings.size() == 1u) {
    std::vector<Definitions::Side> sides =
        Definitions::orientationSides(Definitions::orthogonalOrientation(
            Definitions::sideOrientation(siblings[0].side),
            Definitions::Orientation::DEPTH));
    for (auto side : sides)
      for (auto neighbor : neighbors(siblings[0].id, side))
        if (level == _nodes[neighbor].level() &&
            nodeFamily == _nodes[neighbor].addressCode() / indexStep) {
          Definitions::Neighbor neighborData{static_cast<int>(neighbor), side};
          siblings.emplace_back(neighborData);
        }
  } else if (siblings.size() > 1) {
    for (auto neighbor : neighbors(siblings[0].id, siblings[1].side))
      if (level == _nodes[neighbor].level() &&
          nodeFamily == _nodes[neighbor].addressCode() / indexStep) {
        Definitions::Neighbor neighborData{static_cast<int>(neighbor),
                                           siblings[1].side};
        siblings.emplace_back(neighborData);
      }
  }
  std::vector<size_t> siblingsId(1, nodeId);
  for (auto s : siblings)
    siblingsId.emplace_back(s.id);
  return siblingsId;
}

size_t CellGraph2::coarse(size_t nodeId) {
  // To coarse a node we need to first find out its siblings.
  // A node can be coarsed only if we can find a group of 4 with same parent's
  // addres code
  auto siblings = nodeSiblings(nodeId);
  if (siblings.size() != 4)
    return 0;
  // Before making any removal, we need to find the vertices
  // for convenience, lets sort siblings so we have
  //  1  3
  //  0  2
  {
    auto compare = [&](size_t a, size_t b) -> bool {
      if (_nodes[a].addressCode().x() == _nodes[b].addressCode().x())
        return _nodes[a].addressCode().y() < _nodes[b].addressCode().y();
      return _nodes[a].addressCode().x() < _nodes[b].addressCode().x();
    };
    std::sort(siblings.begin(), siblings.end(), compare);
  }
  // now we may compute the vertices
  size_t vertices[4] = {
      _nodes[siblings[0]].vertex(0), _nodes[siblings[2]].vertex(1),
      _nodes[siblings[3]].vertex(2), _nodes[siblings[1]].vertex(3)};
  // with the sibling being found, we need to delete all them and create a new
  // node that covers their summed region, which is equivalent to go up one
  // level on the tree.
  auto compare = [](const Definitions::Neighbor &a,
                    const Definitions::Neighbor &b) -> bool {
    return a.id < b.id || (a.id == b.id && a.side < b.side);
  };
  // we also need to save all neighbors to connect them later to the new node
  std::set<Definitions::Neighbor, decltype(compare)> _neighbors(compare);
  for (auto s : siblings)
    for (auto n : neighbors(s))
      if (n.id != static_cast<int>(siblings[0]) &&
          n.id != static_cast<int>(siblings[1]) &&
          n.id != static_cast<int>(siblings[2]) &&
          n.id != static_cast<int>(siblings[3]))
        _neighbors.insert(n);
  size_t level = _nodes[siblings[0]].level();
  BBox2d parentRegion =
      BBox2d::combine(_nodes[siblings[0]].region(),
                      BBox2d::combine(_nodes[siblings[1]].region(),
                                      _nodes[siblings[2]].region()));
  Point2u addressCode(std::min(_nodes[siblings[0]].addressCode().x(),
                               _nodes[siblings[1]].addressCode().x()),
                      std::min(_nodes[siblings[0]].addressCode().y(),
                               _nodes[siblings[1]].addressCode().y()));
  addressCode =
      Point2u(std::min(addressCode.x(), _nodes[siblings[2]].addressCode().x()),
              std::min(addressCode.y(), _nodes[siblings[2]].addressCode().y()));
  for (auto s : siblings)
    destroyNode(s);
  size_t parent = newNode(vertices, level - 1);
  _nodes[parent].setAddressCode(addressCode);
  for (auto n : _neighbors) {
    BBox2d region = parentRegion;
    Definitions::Side side = Definitions::oppositeSide(n.side);
    if (n.id && _nodes[n.id].level() > level - 1) {
      region = _nodes[n.id].region();
      side = n.side;
    }
    Point2d center = region.center();
    switch (side) {
    case Definitions::Side::TOP:
      center = Point2d(center.x(), region.lower().y());
      break;
    case Definitions::Side::BOTTOM:
      center = Point2d(center.x(), region.upper().y());
      break;
    case Definitions::Side::LEFT:
      center = Point2d(region.upper().x(), center.y());
      break;
    case Definitions::Side::RIGHT:
      center = Point2d(region.lower().x(), center.y());
      break;
    default:
      THROW(false, "invalid neighbor side");
    }
    newEdge(n.side, n.id, parent, center);
  }
  return parent;
}

size_t CellGraph2::nodeCount() const {
  return _nodes.size() - 1 - _availableNodeIds.size();
}

size_t CellGraph2::edgeCount() const {
  return _edges.size() - _availableEdgeIds.size();
}

const CellGraph2::Node &CellGraph2::node(size_t nodeId) const {
  return _nodes[nodeId];
}

const CellGraph2::Edge &CellGraph2::edge(size_t edgeId) const {
  return _edges[edgeId];
}

const CellGraph2::Vertex &CellGraph2::vertex(size_t vertexId) const {
  return _vertices[vertexId];
}

BBox2d CellGraph2::cellRegion(size_t cellId) const {
  THROW(cellId < _nodes.size() && _nodes[cellId].isValid(),
        "CellGraph2::cellRegion invalid cellId");
  return _nodes[cellId].region();
}

void CellGraph2::iterateCells(const std::function<void(size_t)> &f) const {
  for (size_t i = 1; i < _nodes.size(); i++)
    if (_nodes[i].isValid())
      f(_nodes[i].id());
}

void CellGraph2::iterateFaces(const std::function<void(size_t)> &f) const {
  for (size_t i = 0; i < _edges.size(); i++)
    if (_edges[i].isValid())
      f(i);
}

void CellGraph2::iterateVertices(const std::function<void(size_t)> &f) const {
  for (size_t i = 0; i < _vertices.size(); i++)
    if (_vertices[i].isValid())
      f(i);
}

#ifdef _USE_OPENMP
void CellGraph2::iterateCells_par(const std::function<void(size_t)> &f) const
{
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_nodes.size();

#pragma omp parallel for
  for (index_type i = 1; i < n; i++)
    if (_nodes[i].isValid())
      f(_nodes[i].id());
}

void CellGraph2::iterateFaces_par(const std::function<void(size_t)> &f) const
{
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_edges.size();

#pragma omp parallel for
  for (index_type i = 0; i < n; i++)
    if (_edges[i].isValid())
      f((size_t)i);
}

void CellGraph2::iterateVertices_par(const std::function<void(size_t)> &f) const
{
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_vertices.size();

#pragma omp parallel for
  for (index_type i = 0; i < n; i++)
    if (_vertices[i].isValid())
      f((size_t)i);
}
#endif // _USE_OPENMP

int CellGraph2::cellId(Point2d target) const {
  for (size_t curNode = 1; curNode < _nodes.size(); curNode++)
    if (_nodes[curNode].isValid() && _nodes[curNode].region().contains(target))
      return curNode;
  return -1;
  /*
// get any node connected to node 0
auto randomNodes = _edges[*_nodes[0].edges().begin()].nodes();
size_t curNode = (randomNodes[0]) ? randomNodes[0] : randomNodes[1];
while (!_nodes[curNode].region().contains(target)) {
  std::cerr << "curNode " << curNode << " " << _nodes[curNode].region() << " ---
" << target << std::endl; Definitions::Side targetDirection =
Definitions::Side::RIGHT; auto region = _nodes[curNode].region(); if (target.x()
< region.lower().x()) targetDirection = Definitions::Side::LEFT; else if
(target.x() > region.upper().x()) targetDirection = Definitions::Side::RIGHT;
  else if (target.y() > region.upper().y())
    targetDirection = Definitions::Side::TOP;
  else if (target.y() < region.lower().y())
    targetDirection = Definitions::Side::BOTTOM;
  else THROW(false, "CellGraph2::cellId: target not inside cell!");
  auto _neighbors = neighbors(curNode, targetDirection);
  if (!_neighbors[0])
    return -1;
  curNode = _neighbors[0];
}
return static_cast<int>(curNode);*/
}

std::vector<Definitions::Neighbor>
CellGraph2::cellNeighbors(size_t cellId) const {
  THROW(cellId < _nodes.size(), "invalid node id.");
  std::vector<Definitions::Neighbor> _neighbors;
  for (auto &n : neighbors(cellId)) {
    if (!n.id)
      n.id = -1;
    _neighbors.emplace_back(n);
  }
  return _neighbors;
}

std::vector<Definitions::Neighbor> CellGraph2::faceCells(size_t faceId) const {
  std::vector<Definitions::Neighbor> cells;
  auto nodes = _edges[faceId].nodes();
  for (auto n : nodes) {
    auto neighbor = _edges[faceId].nodeNeighbor(n);
    if (neighbor.id)
      cells.emplace_back(neighbor);
  }
  return cells;
}

Definitions::Orientation CellGraph2::faceOrientation(size_t faceId) const {
  return Definitions::orthogonalOrientation(_edges[faceId].orientation(),
                                            Definitions::Orientation::DEPTH);
}

size_t CellGraph2::cellCount() const { return nodeCount(); }

size_t CellGraph2::faceCount() const { return edgeCount(); }

size_t CellGraph2::vertexCount() const {
  return _vertices.size() - _availableVertexIds.size();
}

Point2d CellGraph2::cellCenterPosition(size_t cellId) const {
  return _nodes[cellId].region().center();
}

Point2d CellGraph2::faceCenterPosition(size_t faceId) const {
  return _edges[faceId].center();
}

Point2d CellGraph2::vertexPosition(size_t vertexId) const {
  return _vertices[vertexId].position();
}

void CellGraph2::iterateNeighborFaces(
    Point2d center, size_t star, Definitions::Orientation orientationFilter,
    const std::function<void(size_t)> &f) const {
  int centerCell = 0;
  try {
    centerCell = cellId(center);
  } catch (const char *e) {
    std::cerr << e << std::endl;
    throw "CellGraph2::neighborFaces: center outside domain!";
  }
  std::set<size_t> faces;
  std::set<size_t> visitedCells;
  struct CellData {
    size_t id, starLevel;
  };
  std::queue<CellData> q;
  q.push({static_cast<size_t>(centerCell), 0});
  visitedCells.insert(centerCell);
  while (!q.empty()) {
    auto curCell = q.front();
    q.pop();
    auto curNode = node(curCell.id);
    for (auto e : curNode.edges()) {
      if (faceOrientation(e) == orientationFilter)
        faces.insert(e);
      int neighborCell = edge(e).nodeNeighbor(curCell.id).id;
      if (curCell.starLevel + 1 <= star && neighborCell > 0 &&
          visitedCells.find(neighborCell) == visitedCells.end()) {
        visitedCells.insert(neighborCell);
        q.push({static_cast<size_t>(neighborCell), curCell.starLevel + 1});
      }
    }
  }
  for (auto face : faces)
    f(face);
}

double CellGraph2::sampleFaceField(
    size_t fieldId, Point2d target, Definitions::Orientation orientationFilter,
    const std::function<bool(size_t)> &faceIsValid,
    const std::function<bool(size_t)> &cellIsValid,
    Definitions::BoundaryVelocity boundaryVelocity) const {
  std::vector<Point2d> points;
  std::vector<double> f;
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  iterateNeighborFaces(target, 1, orientationFilter, [&](size_t id) {
    if (faceIsValid && !faceIsValid(id))
      return;
    f.emplace_back(field[id]);
    points.emplace_back(faceCenterPosition(id));
  });
  if (boundaryVelocity != Definitions::BoundaryVelocity::NONE) {
    bool useVertex[4] = {false, false, false, false};
    int centerId = cellId(target);
    THROW(centerId > 0, "CellGraph<>::sampleFaceField target outside domain");
    double vertexValue[4] = {0., 0., 0., 0.};
    //   0    1     2     3
    // LEFT RIGHT BOTTOM TOP
    double faceValues[4];
    if (_nodes[centerId].edges().size() == 4) {
      for (auto edgeId : _nodes[centerId].edges())
        switch (_edges[edgeId].nodeNeighbor(centerId).side) {
        case Definitions::Side::LEFT:
          faceValues[1] = field[edgeId];
          break;
        case Definitions::Side::RIGHT:
          faceValues[0] = field[edgeId];
          break;
        case Definitions::Side::BOTTOM:
          faceValues[3] = field[edgeId];
          break;
        case Definitions::Side::TOP:
          faceValues[2] = field[edgeId];
          break;
        }
      auto _neighbors = neighbors(cellId(target));
      for (auto neighbor : _neighbors) {
        if (neighbor.id < 0 || !cellIsValid(neighbor.id)) {
          switch (neighbor.side) {
          case Definitions::Side::LEFT:
            useVertex[0] = useVertex[3] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[0] = faceValues[2];
              vertexValue[3] = faceValues[3];
            }
            break;
          case Definitions::Side::RIGHT:
            useVertex[1] = useVertex[2] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[1] = faceValues[2];
              vertexValue[2] = faceValues[3];
            }
            break;
          case Definitions::Side::BOTTOM:
            useVertex[0] = useVertex[1] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[0] = faceValues[0];
              vertexValue[1] = faceValues[1];
            }
            break;
          case Definitions::Side::TOP:
            useVertex[2] = useVertex[3] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[2] = faceValues[0];
              vertexValue[3] = faceValues[1];
            }
            break;
          }
        }
      }
      for (size_t i = 0; i < 4; i++)
        if (useVertex[i]) {
          points.emplace_back(vertexPosition(_nodes[centerId].vertex(i)));
          f.emplace_back(vertexValue[i]);
        }
    }
  }

  if (points.size() == 0)
    std::cerr << "target " << target << std::endl;
  THROW(points.size() != 0,
        "CellGraph2::sampleFaceField No points to inteporlate");

  return _interpolator->interpolateAt(target, points, f);
}

double
CellGraph2::sampleCellField(size_t fieldId, Point2d target,
                            const std::function<bool(size_t)> &isValid) const {
  std::vector<Point2d> points;
  std::vector<double> f;
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  int cellIndex = cellId(target);
  THROW(cellIndex > 0, "CellGraph2::sampleCellField target outside domain");
  iterateNeighborCells(cellIndex, 1, [&](size_t id) {
    if (isValid && !isValid(id))
      return;
    f.emplace_back(field[id]);
    points.emplace_back(cellCenterPosition(id));
  });
  auto _neighbors = neighbors(cellIndex);
  for (auto neighbor : _neighbors) {
    if (!neighbor.id || (isValid && !isValid(neighbor.id))) {
      f.emplace_back(0.);
      switch (neighbor.side) {
      case Definitions::Side::LEFT:
        points.emplace_back(
            0.5 * (_nodes[cellIndex].vertex(0) + _nodes[cellIndex].vertex(3)));
        break;
      case Definitions::Side::RIGHT:
        points.emplace_back(
            0.5 * (_nodes[cellIndex].vertex(1) + _nodes[cellIndex].vertex(2)));
        break;
      case Definitions::Side::BOTTOM:
        points.emplace_back(
            0.5 * (_nodes[cellIndex].vertex(0) + _nodes[cellIndex].vertex(1)));
        break;
      case Definitions::Side::TOP:
        points.emplace_back(
            0.5 * (_nodes[cellIndex].vertex(2) + _nodes[cellIndex].vertex(3)));
        break;
      default:
        THROW(false, "CellGraph2::sampleCellField wrong side neighbor");
      }
    }
  }
  return _interpolator->interpolateAt(target, points, f);
}

double CellGraph2::sampleVertexField(size_t fieldId, Point2d target) const {
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  int nodeId = cellId(target);
  THROW(nodeId > 0, "CellGraph2::sampleVertexField() point outside domain");
  auto &node = _nodes[nodeId];
  std::vector<furoo::Point2d> points;
  /// points are expected in the order:
  /// 2 ---- 3
  /// |      |
  /// 0 ---- 1
  points.emplace_back(_vertices[node.vertex(0)].position());
  points.emplace_back(_vertices[node.vertex(1)].position());
  points.emplace_back(_vertices[node.vertex(3)].position());
  points.emplace_back(_vertices[node.vertex(2)].position());
  std::vector<double> f;
  f.emplace_back(field[node.vertex(0)]);
  f.emplace_back(field[node.vertex(1)]);
  f.emplace_back(field[node.vertex(3)]);
  f.emplace_back(field[node.vertex(2)]);
  return BilinearInterpolant().interpolateAt(target, points, f);
}

std::vector<size_t> CellGraph2::cellFaces(size_t cellId) const {
  std::vector<size_t> faces;
  for (auto e : _nodes[cellId].edges())
    faces.emplace_back(e);
  return faces;
}

std::vector<size_t> CellGraph2::cellVertices(size_t cellId) const {
  std::vector<size_t> vertices;
  for (auto v : _nodes[cellId].vertices())
    vertices.emplace_back(v);
  return vertices;
}

BBox2d CellGraph2::domainRegion() const { return _domainRegion; }

void CellGraph2::iterateNeighborCells(
    size_t cellId, size_t star, const std::function<void(size_t)> &f) const {
  std::set<size_t> visitedCells;
  struct CellData {
    size_t id, starLevel;
  };
  std::queue<CellData> q;
  q.push({cellId, 0});
  visitedCells.insert(cellId);
  while (!q.empty()) {
    auto curCell = q.front();
    f(curCell.id);
    q.pop();
    auto curNode = node(curCell.id);
    for (auto e : curNode.edges()) {
      int neighborCell = edge(e).nodeNeighbor(curCell.id).id;
      if (curCell.starLevel + 1 <= star && neighborCell > 0 &&
          visitedCells.find(neighborCell) == visitedCells.end()) {
        visitedCells.insert(neighborCell);
        q.push({static_cast<size_t>(neighborCell), curCell.starLevel + 1});
      }
    }
  }
}

void CellGraph2::setGraph(const std::vector<CellGraph2::Vertex> &vertices,
                          const std::vector<CellGraph2::Node> &nodes,
                          const std::vector<CellGraph2::Edge> &edges) {
  _availableVertexIds = std::queue<size_t>();
  _availableNodeIds = std::queue<size_t>();
  _availableEdgeIds = std::queue<size_t>();
  _vertices = vertices;
  _nodes = nodes;
  _edges = edges;
  for (size_t i = 0; i < _nodes.size(); i++) {
    if (!_nodes[i].isValid())
      _availableNodeIds.push(i);
    else
      for (size_t j = 0; j < 4; j++)
        _vertices[_nodes[i].vertex(j)].addCell();
  }
  for (size_t i = 0; i < _vertices.size(); i++)
    if (!_vertices[i].isValid())
      _availableVertexIds.push(i);
  for (size_t i = 0; i < _edges.size(); i++)
    if (!_edges[i].isValid())
      _availableEdgeIds.push(i);
    else {
      if (_edges[i].nodes()[0] == 0)
        _nodes[0].addEdge(i);
      else if (_edges[i].nodes()[1] == 0)
        _nodes[0].addEdge(i);
    }
}

void CellGraph2::middleNodeFaceVertices(size_t node, int *vertices) {
  // vertices are in order LEFT RIGHT BOTTOM TOP
  // remember: vertex stored in a node are
  //  3 - 2
  //  0 - 1
  auto nodeCenter = _nodes[node].region().center();
  auto &edges = _nodes[node].edges();
  for (size_t i = 0; i < 4u; i++)
    vertices[i] = -1;
  for (auto e : edges) {
    auto neighbor = _edges[e].nodeNeighbor(node);
    if (_nodes[neighbor.id].level() > _nodes[node].level()) {
      switch (neighbor.side) {
      case Definitions::Side::LEFT: {
        size_t neighborVertices[2] = {_nodes[neighbor.id].vertex(1),
                                      _nodes[neighbor.id].vertex(2)};
        if (vertices[0] < 0)
          vertices[0] = static_cast<int>(neighborVertices[0]);
        for (size_t neighborVertex : neighborVertices)
          if (std::fabs(_vertices[vertices[0]].position().y() -
                        nodeCenter.y()) >
              std::fabs(_vertices[neighborVertex].position().y() -
                        nodeCenter.y()))
            vertices[0] = static_cast<int>(neighborVertex);
      } break;
      case Definitions::Side::RIGHT: {
        size_t neighborVertices[2] = {_nodes[neighbor.id].vertex(3),
                                      _nodes[neighbor.id].vertex(0)};
        if (vertices[1] < 0)
          vertices[1] = static_cast<int>(neighborVertices[0]);
        for (size_t neighborVertex : neighborVertices)
          if (std::fabs(_vertices[vertices[1]].position().y() -
                        nodeCenter.y()) >
              std::fabs(_vertices[neighborVertex].position().y() -
                        nodeCenter.y()))
            vertices[1] = static_cast<int>(neighborVertex);
      } break;
      case Definitions::Side::BOTTOM: {
        size_t neighborVertices[2] = {_nodes[neighbor.id].vertex(3),
                                      _nodes[neighbor.id].vertex(2)};
        if (vertices[2] < 0)
          vertices[2] = static_cast<int>(neighborVertices[0]);
        for (size_t neighborVertex : neighborVertices)
          if (std::fabs(_vertices[vertices[2]].position().x() -
                        nodeCenter.x()) >
              std::fabs(_vertices[neighborVertex].position().x() -
                        nodeCenter.x()))
            vertices[2] = static_cast<int>(neighborVertex);
      } break;
      case Definitions::Side::TOP: {
        size_t neighborVertices[2] = {_nodes[neighbor.id].vertex(0),
                                      _nodes[neighbor.id].vertex(1)};
        if (vertices[3] < 0)
          vertices[3] = static_cast<int>(neighborVertices[0]);
        for (size_t neighborVertex : neighborVertices)
          if (std::fabs(_vertices[vertices[3]].position().x() -
                        nodeCenter.x()) >
              std::fabs(_vertices[neighborVertex].position().x() -
                        nodeCenter.x()))
            vertices[3] = static_cast<int>(neighborVertex);
      } break;
      default:
        THROW(false, "invalid neighbor side.");
      }
    }
  }
}

std::vector<size_t> CellGraph2::createChildren(size_t nodeId) {
  auto parentLevel = _nodes[nodeId].level();
  auto parentRegionCenter = _nodes[nodeId].region().center();
  // compute children regions
  auto childrenRegions = BBox2d::halfSplit(_nodes[nodeId].region());
  // TODO: fix children indices so these lines are not necessary
  std::swap(childrenRegions[0], childrenRegions[2]);
  std::swap(childrenRegions[1], childrenRegions[3]);
  // create vertices
  //                     topVertex
  //            p3 ------------------ p2
  //              |   c0    |    c1  |
  //              |         |        |
  // leftVertex    |---centerVertex---|  rightVertex
  //              |         |        |
  //              |   c2    |   c3   |
  //              --------------------
  //             p0   bottomVertex   p1
  int faceVertices[4] = {-1, -1, -1, -1};
  middleNodeFaceVertices(nodeId, faceVertices);
  size_t centerVertex = newVertex(parentRegionCenter);
  size_t topVertex = (faceVertices[3] < 0)
                         ? newVertex(childrenRegions[0].upper())
                         : static_cast<size_t>(faceVertices[3]);
  size_t leftVertex = (faceVertices[0] < 0)
                          ? newVertex(childrenRegions[0].lower())
                          : static_cast<size_t>(faceVertices[0]);
  size_t rightVertex = (faceVertices[1] < 0)
                           ? newVertex(childrenRegions[3].upper())
                           : static_cast<size_t>(faceVertices[1]);
  size_t bottomVertex = (faceVertices[2] < 0)
                            ? newVertex(childrenRegions[3].lower())
                            : static_cast<size_t>(faceVertices[2]);
  std::vector<size_t> children;
  {
    size_t vertices[4] = {leftVertex, centerVertex, topVertex,
                          _nodes[nodeId].vertex(3)};
    children.emplace_back(newNode(vertices, parentLevel + 1));
  }
  {
    size_t vertices[4] = {centerVertex, rightVertex, _nodes[nodeId].vertex(2),
                          topVertex};
    children.emplace_back(newNode(vertices, parentLevel + 1));
  }
  {
    size_t vertices[4] = {_nodes[nodeId].vertex(0), bottomVertex, centerVertex,
                          leftVertex};
    children.emplace_back(newNode(vertices, parentLevel + 1));
  }
  {
    size_t vertices[4] = {bottomVertex, _nodes[nodeId].vertex(1), rightVertex,
                          centerVertex};
    children.emplace_back(newNode(vertices, parentLevel + 1));
  }
  return children;
}

void CellGraph2::interconnectChildren(const std::vector<size_t> &children) {
  // interconnect children
  // 0  --  1
  // |      |
  // 2  --  3
  double middleX = _nodes[children[0]].region().upper().x();
  double middleY = _nodes[children[0]].region().lower().y();
  newEdge(Definitions::Side::LEFT, children[0], children[1],
          Point2d(middleX, _nodes[children[0]].region().center().y()));
  newEdge(Definitions::Side::LEFT, children[2], children[3],
          Point2d(middleX, _nodes[children[2]].region().center().y()));
  newEdge(Definitions::Side::TOP, children[0], children[2],
          Point2d(_nodes[children[0]].region().center().x(), middleY));
  newEdge(Definitions::Side::TOP, children[1], children[3],
          Point2d(_nodes[children[1]].region().center().x(), middleY));
}

void CellGraph2::connectChildrenToSmallNeighbor(
    const std::vector<size_t> &children, Definitions::Neighbor neighbor) {
  THROW(_nodes[neighbor.id].level() >= _nodes[children[0]].level(),
        "CellGraph2::connectChildrenToSmallNeighbor neighbor is not bigger "
        "than children.");
  // interconnect children
  // 0  --  1
  // |      |
  // 2  --  3
  double middleX = _nodes[children[0]].region().upper().x();
  double middleY = _nodes[children[0]].region().lower().y();
  auto neighborRegion = _nodes[neighbor.id].region();
  switch (neighbor.side) {
  case Definitions::Side::LEFT:
    newEdge(Definitions::oppositeSide(neighbor.side),
            (neighborRegion.center().y() > middleY) ? children[0] : children[2],
            neighbor.id,
            Point2d(neighborRegion.upper().x(), neighborRegion.center().y()));
    break;
  case Definitions::Side::RIGHT:
    newEdge(Definitions::oppositeSide(neighbor.side),
            (neighborRegion.center().y() > middleY) ? children[1] : children[3],
            neighbor.id,
            Point2d(neighborRegion.lower().x(), neighborRegion.center().y()));
    break;
  case Definitions::Side::BOTTOM:
    newEdge(Definitions::oppositeSide(neighbor.side),
            (neighborRegion.center().x() < middleX) ? children[2] : children[3],
            neighbor.id,
            Point2d(neighborRegion.center().x(), neighborRegion.upper().y()));
    break;
  case Definitions::Side::TOP:
    newEdge(Definitions::oppositeSide(neighbor.side),
            (neighborRegion.center().x() < middleX) ? children[0] : children[1],
            neighbor.id,
            Point2d(neighborRegion.center().x(), neighborRegion.lower().y()));
    break;
  default:
    THROW(false, "CellGraph2::connectChildrenToSmallNeighbor invalid side!");
  }
}

void CellGraph2::connectChildrenToBigNeighbor(
    const std::vector<size_t> &children, Definitions::Neighbor neighbor) {
  THROW(_nodes[neighbor.id].level() < _nodes[children[0]].level(),
        "CellGraph2::connectChildrenToBigNeighbor neighbor is not bigger than "
        "children.");
  switch (neighbor.side) {
  case Definitions::Side::RIGHT:
    newEdge(Definitions::Side::LEFT, children[1], neighbor.id,
            Point2d(_nodes[children[1]].region().upper().x(),
                    _nodes[children[1]].region().center().y()));
    newEdge(Definitions::Side::LEFT, children[3], neighbor.id,
            Point2d(_nodes[children[3]].region().upper().x(),
                    _nodes[children[3]].region().center().y()));
    break;
  case Definitions::Side::LEFT:
    newEdge(Definitions::Side::RIGHT, children[0], neighbor.id,
            Point2d(_nodes[children[0]].region().lower().x(),
                    _nodes[children[0]].region().center().y()));
    newEdge(Definitions::Side::RIGHT, children[2], neighbor.id,
            Point2d(_nodes[children[2]].region().lower().x(),
                    _nodes[children[2]].region().center().y()));
    break;
  case Definitions::Side::TOP:
    newEdge(Definitions::Side::BOTTOM, children[0], neighbor.id,
            Point2d(_nodes[children[0]].region().center().x(),
                    _nodes[children[0]].region().upper().y()));
    newEdge(Definitions::Side::BOTTOM, children[1], neighbor.id,
            Point2d(_nodes[children[1]].region().center().x(),
                    _nodes[children[1]].region().upper().y()));
    break;
  case Definitions::Side::BOTTOM:
    newEdge(Definitions::Side::TOP, children[2], neighbor.id,
            Point2d(_nodes[children[2]].region().center().x(),
                    _nodes[children[2]].region().lower().y()));
    newEdge(Definitions::Side::TOP, children[3], neighbor.id,
            Point2d(_nodes[children[3]].region().center().x(),
                    _nodes[children[3]].region().lower().y()));
    break;
  default:
    THROW(false, "CellGraph2::connectChildrenToBigNeighbor invalid side!");
  }
}

void CellGraph2::transferFaceFieldToVertexField(
    size_t faceFieldId, size_t vertexFieldId,
    Definitions::Orientation orientationFilter,
    const std::function<bool(size_t)> &isValid) {
  auto &vertexField =
      *dynamic_cast<Field<double> *>(this->_fields[vertexFieldId].get());
  for (auto &v : _vertices)
    if (v.isValid()) {
      vertexField[v.id()] = sampleFaceField(faceFieldId, v.position(),
                                            orientationFilter, isValid);
    }
}

double CellGraph2::maxScalarFieldValue(size_t fieldId) const {
  THROW(fieldId < this->_fields.size(),
        "CellGraph2::maxScalarFieldValue invalid field id");
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  double value = -INFINITY;
  switch (this->_fieldsLocations[fieldId]) {
  case Definitions::MeshLocation::CELL_CENTER:
    iterateCells([&](size_t i) { value = std::max(value, field[i]); });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
    iterateFaces([&](size_t i) { value = std::max(value, field[i]); });
    break;
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    iterateFaces([&](size_t i) {
      if (faceOrientation(i) == Definitions::Orientation::HORIZONTAL)
        value = std::max(value, field[i]);
    });
    break;
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    iterateFaces([&](size_t i) {
      if (faceOrientation(i) == Definitions::Orientation::VERTICAL)
        value = std::max(value, field[i]);
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    iterateVertices([&](size_t i) { value = std::max(value, field[i]); });
    break;
  default:
    THROW(false, "CellGraph2::maxScalarFieldValue invalid mesh location")
  }
  return value;
}

double CellGraph2::minScalarFieldValue(size_t fieldId) const {
  THROW(fieldId < this->_fields.size(),
        "CellGraph2::minScalarFieldValue invalid field id");
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  double value = INFINITY;
  switch (this->_fieldsLocations[fieldId]) {
  case Definitions::MeshLocation::CELL_CENTER:
    iterateCells([&](size_t i) { value = std::min(value, field[i]); });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
    iterateFaces([&](size_t i) { value = std::min(value, field[i]); });
    break;
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    iterateFaces([&](size_t i) {
      // TODO: check face orientation
      value = std::min(value, field[i]);
    });
    break;
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    iterateFaces([&](size_t i) {
      // TODO: check face orientation
      value = std::min(value, field[i]);
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    iterateVertices([&](size_t i) { value = std::min(value, field[i]); });
    break;
  default:
    THROW(false, "CellGraph2::minScalarFieldValue invalid mesh location")
  }
  return value;
}

void CellGraph2::transferCellFieldToFaceField(
    size_t cellFieldId, size_t faceFieldId,
    const std::function<bool(size_t)> &cellIsValid,
    const std::function<bool(size_t)> &faceIsValid) {
  auto &faceField =
      *dynamic_cast<Field<double> *>(this->_fields[faceFieldId].get());
  iterateFaces([&](size_t i) {
    switch (_fieldsLocations[faceFieldId]) {
    case Definitions::MeshLocation::FACE_CENTER:
      if (faceIsValid(i))
        faceField[i] =
            sampleCellField(faceFieldId, faceCenterPosition(i), cellIsValid);
      break;
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
      if (faceOrientation(i) == Definitions::Orientation::VERTICAL &&
          faceIsValid(i))
        faceField[i] =
            sampleCellField(faceFieldId, faceCenterPosition(i), cellIsValid);
      break;
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      if (faceOrientation(i) == Definitions::Orientation::HORIZONTAL &&
          faceIsValid(i))
        faceField[i] =
            sampleCellField(faceFieldId, faceCenterPosition(i), cellIsValid);
      break;
    default:
      THROW(false, "CellGraph2::transferCellFieldToFaceField face field not "
                   "located in faces");
    }
  });
}

double
CellGraph2::sampleField(size_t fieldId, Point2d target,
                        const std::function<bool(size_t)> &isValid) const {
  return 0;
}

void CellGraph2::iterateCellRing(size_t cellId,
                                 const std::function<void(size_t)> &f) const {
  THROW(cellId < _nodes.size() && _nodes[cellId].isValid(),
        "CellGraph2::iterateCellRing invalid center cell");
  auto _neighbors = neighbors(cellId);
  std::set<size_t> ring;
  for (auto neighbor : _neighbors) {
    if (neighbor.id <= 0 || !_nodes[neighbor.id].isValid())
      continue;
    ring.insert(neighbor.id);
    auto _nNeighbors = neighbors(neighbor.id);
    for (auto nNeighbor : _nNeighbors)
      if (nNeighbor.id > 0 && _nodes[nNeighbor.id].isValid() &&
          Definitions::sideOrientation(nNeighbor.side) !=
              Definitions::sideOrientation(neighbor.side)) {
        // TODO: since we don't store cells in vertices, we need to check
        // neighbors that share the same vertex
        bool found = false;
        for (size_t va = 0; va < 4u && !found; va++)
          for (size_t vb = 0; vb < 4u && !found; vb++)
            if (_nodes[cellId].vertex(va) == _nodes[nNeighbor.id].vertex(vb))
              found = true;
        if (found)
          ring.insert(nNeighbor.id);
      }
  }
  for (auto n : ring)
    f(n);
}

size_t CellGraph2::allCellDataCount() const { return _nodes.size(); }

size_t CellGraph2::allFaceDataCount() const { return _edges.size(); }

size_t CellGraph2::allVertexDataCount() const { return _vertices.size(); }

void CellGraph2::print() const {
  std::cout << "=======================  CellGraph Structure "
               "============================\n";
  std::cout << "All Nodes: " << _nodes.size() << std::endl;
  std::cout << "Valid Nodes: " << cellCount() << std::endl;
  std::cout << "Available ids: ";
  std::queue<size_t> availableNodeIds = _availableNodeIds;
  while (!availableNodeIds.empty()) {
    std::cout << availableNodeIds.front() << ",";
    availableNodeIds.pop();
  }
  std::cout << std::endl;
  std::cout << "List [id valid [vertices] [edges]]\n";
  for (auto n : _nodes) {
    std::cout << n.id() << "\t" << n.isValid() << "\t[";
    for (size_t i = 0; i < 4; i++)
      std::cout << n.vertex(i) << ",";
    std::cout << "]\t[";
    for (auto e : n.edgesCopy())
      std::cout << e << ",";
    std::cout << "]\n";
  }
  std::cout << "All Edges: " << _edges.size() << std::endl;
  std::cout << "Valid Edges: " << faceCount() << std::endl;
  std::cout << "Available ids: ";
  std::queue<size_t> availableEdgeIds = _availableEdgeIds;
  while (!availableEdgeIds.empty()) {
    std::cout << availableEdgeIds.front() << ",";
    availableEdgeIds.pop();
  }
  std::cout << std::endl;
  std::cout << "List [id valid orientation |  aside a -- b bside]\n";
  for (auto e : _edges) {
    std::cout << e.id() << "\t" << e.isValid() << "\t" << e.orientation()
              << " | ";
    std::cout << e.nodeNeighbor(e.nodes()[0]).side << " " << e.nodes()[0]
              << " -- ";
    std::cout << e.nodeNeighbor(e.nodes()[1]).side << " " << e.nodes()[1]
              << std::endl;
  }
  std::cout
      << "==============================================================\n";
}

} // namespace furoo
