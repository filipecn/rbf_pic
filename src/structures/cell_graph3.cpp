#include "cell_graph.h"

namespace furoo {

const CellGraph3::Node &CellGraph3::node(size_t nodeId) const {
  return _nodes[nodeId];
}

const CellGraph3::Edge &CellGraph3::edge(size_t edgeId) const {
  return _edges[edgeId];
}

const CellGraph3::Vertex &CellGraph3::vertex(size_t vertexId) const {
  return _vertices[vertexId];
}

size_t CellGraph3::cellCount() const {
  return _nodes.size() - 1 - _availableNodeIds.size();
}

size_t CellGraph3::faceCount() const {
  return _edges.size() - _availableEdgeIds.size();
}

size_t CellGraph3::vertexCount() const {
  return _vertices.size() - _availableVertexIds.size();
}

Point3d CellGraph3::cellCenterPosition(size_t cellId) const {
  return _nodes[cellId].region().center();
}

Point3d CellGraph3::faceCenterPosition(size_t faceId) const {
  return _edges[faceId].center();
}

Point3d CellGraph3::vertexPosition(size_t vertexId) const {
  return _vertices[vertexId].position();
}

BBox3d CellGraph3::cellRegion(size_t cellId) const {
  return _nodes[cellId].region();
}

BBox3d CellGraph3::domainRegion() const { return _domainRegion; }

void CellGraph3::iterateCells(const std::function<void(size_t)> &f) const {
  for (size_t i = 1; i < _nodes.size(); i++)
    if (_nodes[i].isValid())
      f(_nodes[i].id());
}

void CellGraph3::iterateFaces(const std::function<void(size_t)> &f) const {
  for (size_t i = 0; i < _edges.size(); i++)
    if (_edges[i].isValid())
      f((size_t)i);
}

void CellGraph3::iterateVertices(const std::function<void(size_t)> &f) const {
  for (size_t i = 0; i < _vertices.size(); i++)
    if (_vertices[i].isValid())
      f((size_t)i);
}

#ifdef _USE_OPENMP
void CellGraph3::iterateCells_par(const std::function<void(size_t)> &f) const {
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_nodes.size();

#pragma omp for
  for (index_type i = 1; i < n; i++)
    if (_nodes[i].isValid())
      f(_nodes[i].id());
}

void CellGraph3::iterateFaces_par(const std::function<void(size_t)> &f) const {
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_edges.size();

#pragma omp for
  for (index_type i = 0; i < n; i++)
    if (_edges[i].isValid())
      f(i);
}

void CellGraph3::iterateVertices_par(
    const std::function<void(size_t)> &f) const {
  using index_type = long long; // OpenMP requires signed integral type
  const auto n = (index_type)_vertices.size();

#pragma omp for
  for (index_type i = 0; i < n; i++)
    if (_vertices[i].isValid())
      f(i);
}
#endif // _USE_OPENMP

int CellGraph3::cellId(Point3d target) const {
  if (!_domainRegion.contains(target))
    return -1;
  // get any node connected to node 0
  auto randomNodes = _edges[*_nodes[0].edges().begin()].nodes();
  size_t curNode = (randomNodes[0]) ? randomNodes[0] : randomNodes[1];
  while (!_nodes[curNode].region().contains(target)) {
    Definitions::Side targetDirection = Definitions::Side::RIGHT;
    auto region = _nodes[curNode].region();
    if (target.x() < region.lower().x())
      targetDirection = Definitions::Side::LEFT;
    else if (target.x() > region.upper().x())
      targetDirection = Definitions::Side::RIGHT;
    else if (target.y() > region.upper().y())
      targetDirection = Definitions::Side::TOP;
    else if (target.y() < region.lower().y())
      targetDirection = Definitions::Side::BOTTOM;
    else if (target.z() > region.upper().z())
      targetDirection = Definitions::Side::FRONT;
    else if (target.z() < region.lower().z())
      targetDirection = Definitions::Side::BACK;
    else
      THROW(false, "CellGraph2::cellId: target not inside cell!");
    auto _neighbors = cellNeighbors(curNode);
    double minDistance = INFINITY;
    for (auto &neighbor : _neighbors)
      if (neighbor.side == targetDirection) {
        if (neighbor.id > 0) {
          auto center = _nodes[neighbor.id].region().center();
          Vector3d v = target - center;
          double distance = v[0] + v[1] + v[2];
          if (distance < minDistance) {
            curNode = neighbor.id;
            minDistance = distance;
          }
        }
      }
  }
  return static_cast<int>(curNode);
}

Definitions::Orientation CellGraph3::faceOrientation(size_t faceId) const {
  switch (_edges[faceId].orientation()) {
  case Definitions::Orientation::HORIZONTAL:
    return Definitions::Orientation::VERTICAL;
  case Definitions::Orientation::VERTICAL:
    return Definitions::Orientation::HORIZONTAL;
  case Definitions::Orientation::DEPTH:
    return Definitions::Orientation::DEPTH;
  default:
    THROW(false, "CellGraph3::faceOrienation invalid edge orientation.")
  }
  return Definitions::Orientation::CUSTOM;
}

std::vector<Definitions::Neighbor>
CellGraph3::cellNeighbors(size_t cellId) const {
  THROW(cellId < _nodes.size(), "CellGraph3::cellNeighbors invalid node id.");
  std::vector<Definitions::Neighbor> _neighbors;
  for (auto &n : neighbors(cellId)) {
    if (!n.id)
      n.id = -1;
    _neighbors.emplace_back(n);
  }
  return _neighbors;
}

std::vector<Definitions::Neighbor> CellGraph3::faceCells(size_t faceId) const {
  std::vector<Definitions::Neighbor> cells;
  auto nodes = _edges[faceId].nodes();
  for (auto n : nodes) {
    auto neighbor = _edges[faceId].nodeNeighbor(n);
    if (neighbor.id)
      cells.emplace_back(neighbor);
  }
  return cells;
}

std::vector<size_t> CellGraph3::cellFaces(size_t cellId) const {
  std::vector<size_t> faces;
  for (auto e : _nodes[cellId].edges())
    faces.emplace_back(e);
  return faces;
}

std::vector<size_t> CellGraph3::cellVertices(size_t cellId) const {
  std::vector<size_t> vertices;
  for (auto v : _nodes[cellId].vertices())
    vertices.emplace_back(v);
  return vertices;
}

void CellGraph3::iterateNeighborFaces(
    Point3d center, size_t star, Definitions::Orientation orientationFilter,
    const std::function<void(size_t)> &f) const {
  int centerCell = 0;
  centerCell = cellId(center);
  THROW(centerCell >= 0, "CellGraph3::neighborFaces: center outside domain!");
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
      if (faceOrientation(e) == orientationFilter ||
          orientationFilter == Definitions::Orientation::ANY)
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

void CellGraph3::iterateNeighborCells(
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

double CellGraph3::sampleFaceField(
    size_t fieldId, Point3d target, Definitions::Orientation orientationFilter,
    const std::function<bool(size_t)> &faceIsValid,
    const std::function<bool(size_t)> &cellIsValid,
    Definitions::BoundaryVelocity boundaryVelocity) const {
  std::vector<Point3d> points;
  std::vector<double> f;
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  iterateNeighborFaces(target, 1, orientationFilter, [&](size_t id) {
    if (faceIsValid && !faceIsValid(id))
      return;
    f.emplace_back(field[id]);
    points.emplace_back(faceCenterPosition(id));
  });
  if (boundaryVelocity != Definitions::BoundaryVelocity::NONE) {
    std::vector<bool> useVertex(8, false);
    int centerId = cellId(target);
    THROW(centerId > 0, "CellGraph<>::sampleFaceField target outside domain");
    std::vector<double> vertexValue(8, 0.);
    //   0    1     2     3    4    5
    // LEFT RIGHT BOTTOM TOP BACK FRONT
    double faceValues[6];
    if (_nodes[centerId].edges().size() == 6) {
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
        case Definitions::Side::BACK:
          faceValues[5] = field[edgeId];
          break;
        case Definitions::Side::FRONT:
          faceValues[4] = field[edgeId];
          break;
        }
      auto _neighbors = neighbors(cellId(target));
      for (auto neighbor : _neighbors) {
        if (neighbor.id < 0 || !cellIsValid(neighbor.id)) {
          switch (neighbor.side) {
          case Definitions::Side::LEFT:
            useVertex[0] = useVertex[2] = useVertex[4] = useVertex[6] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[0] = faceValues[2];
              vertexValue[2] = faceValues[3];
              vertexValue[4] = faceValues[2];
              vertexValue[6] = faceValues[3];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_DEPTH) {
              vertexValue[0] = faceValues[4];
              vertexValue[2] = faceValues[4];
              vertexValue[4] = faceValues[5];
              vertexValue[6] = faceValues[5];
            }
            break;
          case Definitions::Side::RIGHT:
            useVertex[1] = useVertex[3] = useVertex[5] = useVertex[7] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[1] = faceValues[2];
              vertexValue[3] = faceValues[3];
              vertexValue[5] = faceValues[2];
              vertexValue[7] = faceValues[3];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_DEPTH) {
              vertexValue[1] = faceValues[4];
              vertexValue[3] = faceValues[4];
              vertexValue[5] = faceValues[5];
              vertexValue[7] = faceValues[5];
            }
            break;
          case Definitions::Side::BOTTOM:
            useVertex[0] = useVertex[1] = useVertex[4] = useVertex[5] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[0] = faceValues[0];
              vertexValue[1] = faceValues[1];
              vertexValue[4] = faceValues[0];
              vertexValue[5] = faceValues[1];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_DEPTH) {
              vertexValue[0] = faceValues[4];
              vertexValue[1] = faceValues[4];
              vertexValue[4] = faceValues[5];
              vertexValue[5] = faceValues[5];
            }
            break;
          case Definitions::Side::TOP:
            useVertex[2] = useVertex[3] = useVertex[7] = useVertex[6] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[2] = faceValues[0];
              vertexValue[3] = faceValues[1];
              vertexValue[6] = faceValues[0];
              vertexValue[7] = faceValues[1];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_DEPTH) {
              vertexValue[2] = faceValues[4];
              vertexValue[3] = faceValues[4];
              vertexValue[6] = faceValues[5];
              vertexValue[7] = faceValues[5];
            }
            break;
          case Definitions::Side::BACK:
            useVertex[2] = useVertex[3] = useVertex[0] = useVertex[1] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[0] = faceValues[0];
              vertexValue[1] = faceValues[1];
              vertexValue[2] = faceValues[0];
              vertexValue[3] = faceValues[1];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[0] = faceValues[2];
              vertexValue[1] = faceValues[2];
              vertexValue[2] = faceValues[3];
              vertexValue[3] = faceValues[3];
            }
            break;
          case Definitions::Side::FRONT:
            useVertex[4] = useVertex[5] = useVertex[6] = useVertex[7] = true;
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_HORIZONTAL) {
              vertexValue[4] = faceValues[0];
              vertexValue[5] = faceValues[1];
              vertexValue[6] = faceValues[0];
              vertexValue[7] = faceValues[1];
            }
            if (boundaryVelocity ==
                Definitions::BoundaryVelocity::FREE_SLIP_VERTICAL) {
              vertexValue[4] = faceValues[2];
              vertexValue[5] = faceValues[2];
              vertexValue[6] = faceValues[3];
              vertexValue[7] = faceValues[3];
            }
            break;
          }
        }
      }
      for (size_t i = 0; i < 8; i++)
        if (useVertex[i]) {
          points.emplace_back(vertexPosition(_nodes[centerId].vertexAt(i)));
          f.emplace_back(vertexValue[i]);
        }
    }
  }
  THROW(points.size() != 0,
        "CellGraph3::sampleFaceField No points to inteporlate");

  return _interpolator->interpolateAt(target, points, f);
}

double
CellGraph3::sampleCellField(size_t fieldId, Point3d target,
                            const std::function<bool(size_t)> &isValid) const {
  int cell = cellId(target);
  THROW(cell > 0, "CellGraph3::sampleCellField target outside domain");
  std::vector<Point3d> points;
  std::vector<double> f;
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  iterateNeighborCells(cell, 1, [&](size_t id) {
    if (isValid && !isValid(id))
      return;
    f.emplace_back(field[id]);
    points.emplace_back(cellCenterPosition(id));
  });
  return _interpolator->interpolateAt(target, points, f);
}

double CellGraph3::sampleVertexField(size_t fieldId, Point3d target) const {
  NOT_IMPLEMENTED();
  return 0.;
}

double
CellGraph3::sampleField(size_t fieldId, Point3d target,
                        const std::function<bool(size_t)> &isValid) const {
  THROW(fieldId < _fields.size(), "CellGraph3::sampleField invalid field id.");
  switch (_fieldsLocations[fieldId]) {
  case Definitions::MeshLocation::CELL_CENTER:
    return sampleCellField(fieldId, target, isValid);
  case Definitions::MeshLocation::FACE_CENTER:
    return sampleFaceField(fieldId, target, Definitions::Orientation::ANY,
                           isValid);
  case Definitions::MeshLocation::VERTEX_CENTER:
    return sampleVertexField(fieldId, target);
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    return sampleFaceField(fieldId, target,
                           Definitions::Orientation::HORIZONTAL, isValid);
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    return sampleFaceField(fieldId, target, Definitions::Orientation::VERTICAL,
                           isValid);
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    return sampleFaceField(fieldId, target, Definitions::Orientation::DEPTH,
                           isValid);
  default:
    THROW(false, "CellGraph3::sampleField invalid location for field.");
  }
  return 0;
}

void CellGraph3::transferFaceFieldToVertexField(
    size_t faceFieldId, size_t vertexFieldId,
    Definitions::Orientation orientationFilter,
    const std::function<bool(size_t)> &isValid) {
  NOT_IMPLEMENTED();
}

void CellGraph3::transferCellFieldToFaceField(
    size_t cellFieldId, size_t faceFieldId,
    const std::function<bool(size_t)> &cellIsValid,
    const std::function<bool(size_t)> &faceIsValid) {
  NOT_IMPLEMENTED();
}

size_t CellGraph3::newNode(const size_t *vertices, size_t level) {
  size_t id = _nodes.size();
  BBox3d r(_vertices[vertices[0]].position(),
           _vertices[vertices[7]].position());
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
  for (size_t i = 0; i < 8; i++)
    _vertices[vertices[i]].addCell();
  return id;
}

size_t CellGraph3::newEdge(Definitions::Side aSide, size_t nodeA, size_t nodeB,
                           Point3d center) {
  size_t id = _edges.size();
  if (!_availableEdgeIds.empty()) {
    id = _availableEdgeIds.front();
    _availableEdgeIds.pop();
    _edges[id].set(aSide, nodeA, nodeB, center);
    for (auto fieldId : this->_faceFieldIds)
      this->_fields[fieldId]->reset(id);
  } else {
    _edges.emplace_back(aSide, nodeA, nodeB, center, id);
    for (auto fieldId : this->_faceFieldIds)
      this->_fields[fieldId]->increaseSize(1);
  }
  if (_nodes[nodeA].isValid())
    _nodes[nodeA].addEdge(id);
  if (_nodes[nodeB].isValid())
    _nodes[nodeB].addEdge(id);
  return id;
}

size_t CellGraph3::newVertex(const Point3d &position) {
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

void CellGraph3::destroyVertex(size_t vertexId) {
  THROW(vertexId < _vertices.size(),
        "CellGraph3::destroyVertex _vertexSize < vertex id")
  THROW(_vertices[vertexId].isValid(),
        "CellGraph3::destroyVertex invalid vertex.");
  _availableVertexIds.push(vertexId);
  _vertices[vertexId].destroy();
  // TODO: destroying a vertex should also destroy its nodes and edges...
}

void CellGraph3::destroyNode(size_t nodeId) {
  THROW(nodeId < _nodes.size(), "CellGraph3::destroyNode _nodesSize < node id")
  THROW(_nodes[nodeId].isValid(), "CellGraph3::destroyNode invalid node.");
  _availableNodeIds.push(nodeId);
  auto edges = _nodes[nodeId].edgesCopy();
  for (auto e : edges)
    destroyEdge(e);
  // decrement cell count of each vertex
  for (size_t i = 0; i < 8; i++) {
    _vertices[_nodes[nodeId].vertexAt(i)].removeCell();
    if (!_vertices[_nodes[nodeId].vertexAt(i)].cellCount())
      destroyVertex(_nodes[nodeId].vertexAt(i));
  }
  _nodes[nodeId].destroy();
}

void CellGraph3::destroyEdge(size_t edgeId) {
  THROW(edgeId < _edges.size() && _edges[edgeId].isValid(),
        "CellGraph3::destroyEdge invalid edge id or trying to destroy an "
        "invalid edge.");
  _availableEdgeIds.push(edgeId);
  auto nodes = _edges[edgeId].nodes();
  _nodes[nodes[0]].removeEdge(edgeId);
  _nodes[nodes[1]].removeEdge(edgeId);
  _edges[edgeId].destroy();
}

std::vector<size_t> CellGraph3::createChildren(size_t nodeId) {
  auto parentLevel = _nodes[nodeId].level();
  auto parentRegionCenter = _nodes[nodeId].region().center();
  // compute children regions
  //  auto childrenRegions = BBox3d::halfSplit(_nodes[nodeId].region());
  // children regions are in z order
  //                     topVertex
  //            p2 ------------------ p3
  //              |   c2    |    c3  |
  //              |         |        |
  // leftVertex    |---centerVertex---|  rightVertex
  //              |         |        |
  //              |   c0    |   c1   |
  //              --------------------
  //             p0   bottomVertex   p1
  std::vector<size_t> children;
  auto splitIndices = splitVertices(nodeId);
  for (size_t childZ = 0; childZ < 2; childZ++)
    for (size_t childY = 0; childY < 2; childY++)
      for (size_t childX = 0; childX < 2; childX++) {
        size_t childOriginVertexId = childZ * 9 + childY * 3 + childX;
        size_t childVertices[8];
        for (size_t z = 0; z < 2; z++)
          for (size_t y = 0; y < 2; y++)
            for (size_t x = 0; x < 2; x++)
              childVertices[z * 4 + y * 2 + x] =
                  splitIndices[childOriginVertexId + 9 * z + 3 * y + x];
        children.emplace_back(newNode(childVertices, parentLevel + 1));
      }
  return children;
}

std::vector<size_t> CellGraph3::splitVertices(size_t nodeId) {
  auto parentLevel = _nodes[nodeId].level();
  auto parentRegionCenter = _nodes[nodeId].region().center();
  auto addressCode = _nodes[nodeId].addressCode();
  size_t indexStep = 1u << (_maxLevel - (parentLevel + 1));
  // The indices of the split vertices
  std::vector<int> vertexIndices(27, -1);
  // The adress codes of the split vertices
  std::vector<Point3u> vertexAddressCodes;
  // Compute split vertices adress codes
  for (size_t z = 0; z <= 2; z++)
    for (size_t y = 0; y <= 2; y++)
      for (size_t x = 0; x <= 2; x++)
        vertexAddressCodes.emplace_back(addressCode +
                                        Vector3u(x, y, z) * indexStep);
  vertexIndices[0] = _nodes[nodeId].vertexAt(0);
  vertexIndices[2] = _nodes[nodeId].vertexAt(1);
  vertexIndices[6] = _nodes[nodeId].vertexAt(2);
  vertexIndices[8] = _nodes[nodeId].vertexAt(3);
  vertexIndices[18] = _nodes[nodeId].vertexAt(4);
  vertexIndices[20] = _nodes[nodeId].vertexAt(5);
  vertexIndices[24] = _nodes[nodeId].vertexAt(6);
  vertexIndices[26] = _nodes[nodeId].vertexAt(7);
  // Stores the address codes of the 8 vertices of a neighbor node
  Point3u neighborAdressCodes[8];
  // Iterate over all node edges in order to find existent vertex indices
  auto edgeIds = _nodes[nodeId].edges();
  for (auto &edgeId : edgeIds) {
    auto neighborId = _edges[edgeId].nodeNeighbor(nodeId).id;
    if (neighborId <= 0)
      continue;
    // Compute the adress codes of all vertices of the current neighbor node
    for (size_t neighborVertexId = 0; neighborVertexId < 8;
         neighborVertexId++) {
      neighborAdressCodes[neighborVertexId] =
          _nodes[neighborId].addressCode(neighborVertexId, _maxLevel);
    }
    // For each of our new vertices, check if there is a neighbor vertex with
    // the same address code
    for (size_t i = 0; i < vertexAddressCodes.size(); i++) {
      if (vertexIndices[i] >= 0)
        continue;
      for (size_t neighborVertexId = 0; neighborVertexId < 8;
           neighborVertexId++)
        if (neighborAdressCodes[neighborVertexId] == vertexAddressCodes[i])
          vertexIndices[i] = _nodes[neighborId].vertexAt(neighborVertexId);
    }
  }
  auto origin = _nodes[nodeId].region().lower();
  double side = _nodes[nodeId].region().size().x() / 2.;
  std::vector<size_t> splitIndices;
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    if (vertexIndices[i] < 0)
      splitIndices.emplace_back(
          newVertex(origin + Vector3d((i % 9) % 3, (i % 9) / 3, i / 9) * side));
    else
      splitIndices.emplace_back(vertexIndices[i]);
  }
  return splitIndices;
}

void CellGraph3::interconnectSiblings(const std::vector<size_t> &children) {
  newEdge(Definitions::Side::LEFT, children[0], children[1],
          Point3d(_nodes[children[0]].region().upper().x(),
                  _nodes[children[0]].region().center().y(),
                  _nodes[children[0]].region().center().z()));
  newEdge(Definitions::Side::LEFT, children[2], children[3],
          Point3d(_nodes[children[2]].region().upper().x(),
                  _nodes[children[2]].region().center().y(),
                  _nodes[children[2]].region().center().z()));
  newEdge(Definitions::Side::LEFT, children[6], children[7],
          Point3d(_nodes[children[6]].region().upper().x(),
                  _nodes[children[6]].region().center().y(),
                  _nodes[children[6]].region().center().z()));
  newEdge(Definitions::Side::LEFT, children[4], children[5],
          Point3d(_nodes[children[4]].region().upper().x(),
                  _nodes[children[4]].region().center().y(),
                  _nodes[children[4]].region().center().z()));
  newEdge(Definitions::Side::BOTTOM, children[0], children[2],
          Point3d(_nodes[children[0]].region().center().x(),
                  _nodes[children[0]].region().upper().y(),
                  _nodes[children[0]].region().center().z()));
  newEdge(Definitions::Side::BOTTOM, children[1], children[3],
          Point3d(_nodes[children[1]].region().center().x(),
                  _nodes[children[1]].region().upper().y(),
                  _nodes[children[1]].region().center().z()));
  newEdge(Definitions::Side::BOTTOM, children[4], children[6],
          Point3d(_nodes[children[4]].region().center().x(),
                  _nodes[children[4]].region().upper().y(),
                  _nodes[children[4]].region().center().z()));
  newEdge(Definitions::Side::BOTTOM, children[5], children[7],
          Point3d(_nodes[children[5]].region().center().x(),
                  _nodes[children[5]].region().upper().y(),
                  _nodes[children[5]].region().center().z()));
  newEdge(Definitions::Side::BACK, children[0], children[4],
          Point3d(_nodes[children[0]].region().center().x(),
                  _nodes[children[0]].region().center().y(),
                  _nodes[children[0]].region().upper().z()));
  newEdge(Definitions::Side::BACK, children[2], children[6],
          Point3d(_nodes[children[2]].region().center().x(),
                  _nodes[children[2]].region().center().y(),
                  _nodes[children[2]].region().upper().z()));
  newEdge(Definitions::Side::BACK, children[3], children[7],
          Point3d(_nodes[children[3]].region().center().x(),
                  _nodes[children[3]].region().center().y(),
                  _nodes[children[3]].region().upper().z()));
  newEdge(Definitions::Side::BACK, children[1], children[5],
          Point3d(_nodes[children[1]].region().center().x(),
                  _nodes[children[1]].region().center().y(),
                  _nodes[children[1]].region().upper().z()));
}

void CellGraph3::connectChildrenToBigNeighbor(
    const std::vector<size_t> &children, Definitions::Neighbor neighbor) {
  THROW(_nodes[neighbor.id].level() < _nodes[children[0]].level(),
        "CellGraph3::connectChildrenToBigNeighbor neighbor is not bigger than "
        "children.");
  switch (neighbor.side) {
  case Definitions::Side::LEFT: {
    size_t ids[4] = {0, 2, 4, 6};
    for (size_t i : ids)
      newEdge(Definitions::Side::RIGHT, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().lower().x(),
                      _nodes[children[i]].region().center().y(),
                      _nodes[children[i]].region().center().z()));
  } break;
  case Definitions::Side::RIGHT: {
    size_t ids[4] = {1, 3, 5, 7};
    for (size_t i : ids)
      newEdge(Definitions::Side::LEFT, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().upper().x(),
                      _nodes[children[i]].region().center().y(),
                      _nodes[children[i]].region().center().z()));
  } break;
  case Definitions::Side::BOTTOM: {
    size_t ids[4] = {0, 1, 4, 5};
    for (size_t i : ids)
      newEdge(Definitions::Side::TOP, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().center().x(),
                      _nodes[children[i]].region().lower().y(),
                      _nodes[children[i]].region().center().z()));
  } break;
  case Definitions::Side::TOP: {
    size_t ids[4] = {2, 3, 6, 7};
    for (size_t i : ids)
      newEdge(Definitions::Side::BOTTOM, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().center().x(),
                      _nodes[children[i]].region().upper().y(),
                      _nodes[children[i]].region().center().z()));
  } break;
  case Definitions::Side::BACK: {
    size_t ids[4] = {0, 1, 2, 3};
    for (size_t i : ids)
      newEdge(Definitions::Side::FRONT, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().center().x(),
                      _nodes[children[i]].region().center().y(),
                      _nodes[children[i]].region().lower().z()));
  } break;
  case Definitions::Side::FRONT: {
    size_t ids[4] = {4, 5, 6, 7};
    for (size_t i : ids)
      newEdge(Definitions::Side::BACK, children[i], neighbor.id,
              Point3d(_nodes[children[i]].region().center().x(),
                      _nodes[children[i]].region().center().y(),
                      _nodes[children[i]].region().upper().z()));
  } break;
  default:
    THROW(false, "CellGraph3::connectChildrenToBigNeighbor invalid side!");
  }
}

void CellGraph3::connectChildrenToSmallNeighbor(
    const std::vector<size_t> &children, Definitions::Neighbor neighbor) {
  auto neighborAddressCode = _nodes[neighbor.id].addressCode();
  switch (neighbor.side) {
  case Definitions::Side::LEFT: {
    auto middleAddressCode = _nodes[children[2]].addressCode(4, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.z() < middleAddressCode.z())
      childId = children[0];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.z() < middleAddressCode.z())
      childId = children[2];
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[4];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[6];
    newEdge(Definitions::Side::RIGHT, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().upper().x(),
                    _nodes[neighbor.id].region().center().y(),
                    _nodes[neighbor.id].region().center().z()));
  } break;
  case Definitions::Side::RIGHT: {
    auto middleAddressCode = _nodes[children[3]].addressCode(5, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.z() < middleAddressCode.z())
      childId = children[1];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.z() < middleAddressCode.z())
      childId = children[3];
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[5];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[7];
    newEdge(Definitions::Side::LEFT, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().lower().x(),
                    _nodes[neighbor.id].region().center().y(),
                    _nodes[neighbor.id].region().center().z()));
  } break;
  case Definitions::Side::BOTTOM: {
    auto middleAddressCode = _nodes[children[0]].addressCode(5, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.x() < middleAddressCode.x() &&
        neighborAddressCode.z() < middleAddressCode.z())
      childId = children[0];
    else if (neighborAddressCode.x() >= middleAddressCode.x() &&
             neighborAddressCode.z() < middleAddressCode.z())
      childId = children[1];
    if (neighborAddressCode.x() < middleAddressCode.x() &&
        neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[4];
    else if (neighborAddressCode.x() >= middleAddressCode.x() &&
             neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[5];
    newEdge(Definitions::Side::TOP, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().center().x(),
                    _nodes[neighbor.id].region().upper().y(),
                    _nodes[neighbor.id].region().center().z()));
  } break;
  case Definitions::Side::TOP: {
    auto middleAddressCode = _nodes[children[0]].addressCode(7, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.x() < middleAddressCode.x() &&
        neighborAddressCode.z() < middleAddressCode.z())
      childId = children[2];
    else if (neighborAddressCode.x() >= middleAddressCode.x() &&
             neighborAddressCode.z() < middleAddressCode.z())
      childId = children[3];
    if (neighborAddressCode.x() < middleAddressCode.x() &&
        neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[6];
    else if (neighborAddressCode.x() >= middleAddressCode.x() &&
             neighborAddressCode.z() >= middleAddressCode.z())
      childId = children[7];
    newEdge(Definitions::Side::BOTTOM, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().center().x(),
                    _nodes[neighbor.id].region().lower().y(),
                    _nodes[neighbor.id].region().center().z()));
  } break;
  case Definitions::Side::BACK: {
    auto middleAddressCode = _nodes[children[0]].addressCode(3, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.x() < middleAddressCode.x())
      childId = children[0];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.x() < middleAddressCode.x())
      childId = children[2];
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.x() >= middleAddressCode.x())
      childId = children[1];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.x() >= middleAddressCode.x())
      childId = children[3];
    newEdge(Definitions::Side::FRONT, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().center().x(),
                    _nodes[neighbor.id].region().center().y(),
                    _nodes[neighbor.id].region().upper().z()));
  } break;
  case Definitions::Side::FRONT: {
    auto middleAddressCode = _nodes[children[4]].addressCode(7, _maxLevel);
    size_t childId = 0;
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.x() < middleAddressCode.x())
      childId = children[4];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.x() < middleAddressCode.x())
      childId = children[6];
    if (neighborAddressCode.y() < middleAddressCode.y() &&
        neighborAddressCode.x() >= middleAddressCode.x())
      childId = children[5];
    else if (neighborAddressCode.y() >= middleAddressCode.y() &&
             neighborAddressCode.x() >= middleAddressCode.x())
      childId = children[7];
    newEdge(Definitions::Side::BACK, childId, neighbor.id,
            Point3d(_nodes[neighbor.id].region().center().x(),
                    _nodes[neighbor.id].region().center().y(),
                    _nodes[neighbor.id].region().lower().z()));
  } break;
  default:
    THROW(false, "CellGraph3::connectChildrenToSmallNeighbor invalid side!");
  }
}

std::vector<size_t> CellGraph3::refine(size_t nodeId) {
  THROW(nodeId < _nodes.size(), "CellGraph3::refine invalid node id.");
  // This code is quite delicate. We need to track update neighbors, edges and
  // vertices! Handling nodes:
  //    - Remove the parent node and replace it by 8 children as in a octree
  // Handling edges:
  //    - All edges are removed, however we need to save the neighboring nodes
  //    first
  //      8 new edges are created to connect the sibling nodes that were just
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
  Point3u parentAddressCode = _nodes[nodeId].addressCode();
  auto parentLevel = _nodes[nodeId].level();
  size_t indexStep = 1u << (_maxLevel - (parentLevel + 1));
  // Compute split vertices adress codes
  for (size_t z = 0; z < 2; z++)
    for (size_t y = 0; y < 2; y++)
      for (size_t x = 0; x < 2; x++)
        _nodes[children[z * 4 + y * 2 + x]].setAddressCode(
            parentAddressCode + Vector3u(x, y, z) * indexStep);
  // create edges between children
  interconnectSiblings(children);
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

void CellGraph3::refine(size_t nodeId,
                        std::function<bool(const CellGraph3::Node &)> f) {
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

std::vector<Definitions::Neighbor> CellGraph3::neighbors(size_t nodeId) const {
  THROW(nodeId < _nodes.size(), "CellGraph3::neighbors invalid node id.");
  auto edges = _nodes[nodeId].edges();
  std::vector<Definitions::Neighbor> _neighbors;
  for (auto e : edges)
    _neighbors.emplace_back(_edges[e].nodeNeighbor(nodeId));
  return _neighbors;
}

CellGraph3::~CellGraph3() = default;

CellGraph3::CellGraph3(BBox3d region, size_t maxLevel)
    : _maxLevel(maxLevel), _domainRegion(region),
      _interpolator(new MLS<3>(Definitions::PolynomialType::QUADRATIC)) {
  // create vertices
  newVertex(
      Point3d(region.lower().x(), region.lower().y(), region.lower().z()));
  newVertex(
      Point3d(region.upper().x(), region.lower().y(), region.lower().z()));
  newVertex(
      Point3d(region.lower().x(), region.upper().y(), region.lower().z()));
  newVertex(
      Point3d(region.upper().x(), region.upper().y(), region.lower().z()));
  newVertex(
      Point3d(region.lower().x(), region.lower().y(), region.upper().z()));
  newVertex(
      Point3d(region.upper().x(), region.lower().y(), region.upper().z()));
  newVertex(
      Point3d(region.lower().x(), region.upper().y(), region.upper().z()));
  newVertex(
      Point3d(region.upper().x(), region.upper().y(), region.upper().z()));
  size_t vertices[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  newNode(vertices, 0);
  _nodes[0].setAddressCode(
      Point3u(2 << maxLevel, 2 << maxLevel, 2 << maxLevel));
  newNode(vertices, 0);
  _nodes[1].setAddressCode(Point3u(0, 0, 0));
  auto regionCenter = region.center();
  // create edges and connect them to root node
  newEdge(Definitions::Side::LEFT, 0, 1,
          Point3d(region.lower().x(), regionCenter.y(), regionCenter.z()));
  newEdge(Definitions::Side::RIGHT, 0, 1,
          Point3d(region.upper().x(), regionCenter.y(), regionCenter.z()));
  newEdge(Definitions::Side::BOTTOM, 0, 1,
          Point3d(regionCenter.x(), region.lower().y(), regionCenter.z()));
  newEdge(Definitions::Side::TOP, 0, 1,
          Point3d(regionCenter.x(), region.upper().y(), regionCenter.z()));
  newEdge(Definitions::Side::BACK, 0, 1,
          Point3d(regionCenter.x(), regionCenter.y(), region.lower().z()));
  newEdge(Definitions::Side::FRONT, 0, 1,
          Point3d(regionCenter.x(), regionCenter.y(), region.upper().z()));
}

CellGraph3::CellGraph3(const BBox3d &region,
                       const std::function<bool(const CellGraph3::Node &)> &f,
                       size_t maxLevel) {
  refine(1, f);
}

bool CellGraph3::canCoarse(std::vector<size_t> nodesId) const {
  if (nodesId.size() != 8)
    return false;
  size_t level = _nodes[nodesId[0]].level();
  size_t indexStep = 1u << (_maxLevel - (level - 1u));
  Point3u nodeFamily = _nodes[nodesId[0]].addressCode() / indexStep;
  for (auto id : nodesId) {
    if (_nodes[id].level() != level)
      return false;
    if (nodeFamily != _nodes[id].addressCode() / indexStep)
      return false;
  }
  return true;
}

std::vector<size_t> CellGraph3::nodeSiblings(size_t nodeId) const {
  size_t level = _nodes[nodeId].level();
  size_t indexStep = 1u << (_maxLevel - (level - 1u));
  Point3u nodeFamily = _nodes[nodeId].addressCode() / indexStep;
  std::vector<size_t> siblings;
  std::set<size_t> visited;
  std::queue<size_t> q;
  q.push(nodeId);
  while (!q.empty()) {
    size_t curNode = q.front();
    if (visited.find(curNode) == visited.end())
      siblings.emplace_back(curNode);
    q.pop();
    visited.insert(curNode);
    auto neighbors = cellNeighbors(curNode);
    for (auto neighbor : neighbors)
      if (neighbor.id > 0 && visited.find(neighbor.id) == visited.end() &&
          level == _nodes[neighbor.id].level() &&
          nodeFamily == _nodes[neighbor.id].addressCode() / indexStep) {
        q.push(neighbor.id);
      }
  }
  return siblings;
}

size_t CellGraph3::coarse(size_t nodeId) {
  // To coarse a node we need to first find out its siblings.
  // A node can be coarsed only if we can find a group of 4 with same parent's
  // addres code
  auto siblings = nodeSiblings(nodeId);
  if (siblings.size() != 8)
    return 0;
  // Before making any removal, we need to find the vertices
  // for convenience, lets sort by x < y < z
  {
    auto compare = [&](size_t a, size_t b) -> bool {
      if (_nodes[a].addressCode().z() == _nodes[b].addressCode().z()) {
        if (_nodes[a].addressCode().y() == _nodes[b].addressCode().y())
          return _nodes[a].addressCode().x() < _nodes[b].addressCode().x();
        return _nodes[a].addressCode().y() < _nodes[b].addressCode().y();
      }
      return _nodes[a].addressCode().z() < _nodes[b].addressCode().z();
    };
    std::sort(siblings.begin(), siblings.end(), compare);
  }
  // now we may compute the vertices
  size_t vertices[8];
  for (size_t i = 0; i < 8; i++)
    vertices[i] = _nodes[siblings[i]].vertexAt(i);
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
    for (auto n : neighbors(s)) {
      bool good = true;
      for (size_t sibling : siblings)
        if (n.id == static_cast<int>(sibling))
          good = false;
      if (good)
        _neighbors.insert(n);
    }
  Point3u addressCode = _nodes[siblings[0]].addressCode();
  size_t level = _nodes[siblings[0]].level();
  // remove all siblings
  for (auto s : siblings)
    destroyNode(s);
  // create the new node and connect neighbors
  size_t parent = newNode(vertices, level - 1);
  _nodes[parent].setAddressCode(addressCode);
  std::vector<Definitions::Neighbor> neighborList;
  for (auto &neighbor : _neighbors)
    neighborList.emplace_back(neighbor);
  connectParentToNeighbors(parent, neighborList);
  return parent;
}

void CellGraph3::connectParentToNeighbors(
    size_t nodeId, const std::vector<Definitions::Neighbor> &_neighbors) {
  BBox3d parentRegion = _nodes[nodeId].region();
  for (auto n : _neighbors) {
    BBox3d region = parentRegion;
    Definitions::Side side = Definitions::oppositeSide(n.side);
    if (n.id && _nodes[n.id].level() > _nodes[nodeId].level()) {
      region = _nodes[n.id].region();
      side = n.side;
    }
    Point3d center = region.center();
    switch (side) {
    case Definitions::Side::LEFT:
      center = Point3d(region.upper().x(), center.y(), center.z());
      break;
    case Definitions::Side::RIGHT:
      center = Point3d(region.lower().x(), center.y(), center.z());
      break;
    case Definitions::Side::BOTTOM:
      center = Point3d(center.x(), region.upper().y(), center.z());
      break;
    case Definitions::Side::TOP:
      center = Point3d(center.x(), region.lower().y(), center.z());
      break;
    case Definitions::Side::BACK:
      center = Point3d(center.x(), center.y(), region.upper().y());
      break;
    case Definitions::Side::FRONT:
      center = Point3d(center.x(), center.y(), region.lower().y());
      break;
    default:
      THROW(false, "invalid neighbor side");
    }
    newEdge(n.side, n.id, nodeId, center);
  }
}

void CellGraph3::iterateCellRing(size_t cellId,
                                 const std::function<void(size_t)> &f) const {
  THROW(cellId < _nodes.size() && _nodes[cellId].isValid(),
        "CellGraph3::iterateCellRing invalid center cell");
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
        for (size_t va = 0; va < 8u && !found; va++)
          for (size_t vb = 0; vb < 8u && !found; vb++)
            if (_nodes[cellId].vertexAt(va) ==
                _nodes[nNeighbor.id].vertexAt(vb))
              found = true;
        if (found)
          ring.insert(nNeighbor.id);
      }
  }
  for (auto n : ring)
    f(n);
}

size_t CellGraph3::allCellDataCount() const { return _nodes.size(); }

size_t CellGraph3::allFaceDataCount() const { return _edges.size(); }

size_t CellGraph3::allVertexDataCount() const { return _vertices.size(); }

void CellGraph3::setGraph(const std::vector<CellGraph3::Vertex> &vertices,
                          const std::vector<CellGraph3::Node> &nodes,
                          const std::vector<CellGraph3::Edge> &edges) {
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
      for (size_t j = 0; j < 8; j++)
        _vertices[_nodes[i].vertexAt(j)].addCell();
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
  // update field sizes
  for (size_t i = 0; i < _fields.size(); i++)
    switch (_fieldsLocations[i]) {
    case Definitions::MeshLocation::CELL_CENTER:
      _fields[i]->increaseSize(_nodes.size());
      break;
    case Definitions::MeshLocation::FACE_CENTER:
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      _fields[i]->increaseSize(_edges.size());
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      _fields[i]->increaseSize(_vertices.size());
      break;
    }
} // namespace furoo

void CellGraph3::print() const {
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
    for (size_t i = 0; i < 8; i++)
      std::cout << n.vertexAt(i) << ",";
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