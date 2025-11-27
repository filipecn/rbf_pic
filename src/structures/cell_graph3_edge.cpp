#include "cell_graph.h"

namespace furoo {

CellGraph3::Edge::Edge(Definitions::Side aSide, size_t vertexA, size_t vertexB,
                       Point3d center, size_t id) : _id(id) {
  set(aSide, vertexA, vertexB, center);
}

void CellGraph3::Edge::set(Definitions::Side aSide,
                           size_t vertexA,
                           size_t vertexB,
                           Point3d center) {
  _nodeIds[0] = vertexA;
  _nodeSides[0] = aSide;
  _nodeIds[1] = vertexB;
  _nodeSides[1] = Definitions::oppositeSide(aSide);
  _center = center;
  _orientation = Definitions::sideOrientation(aSide);
  _isValid = true;
}

Definitions::Orientation CellGraph3::Edge::orientation() const { return _orientation; }

Definitions::Neighbor CellGraph3::Edge::nodeNeighbor(size_t nodeId) const {
  Definitions::Neighbor neighbor{0, Definitions::Side::CUSTOM};
  if (nodeId == _nodeIds[0]) {
    neighbor.id = static_cast<int>(_nodeIds[1]);
    neighbor.side = _nodeSides[1];
  } else if (nodeId == _nodeIds[1]) {
    neighbor.id = static_cast<int>(_nodeIds[0]);
    neighbor.side = _nodeSides[0];
  }
  return neighbor;
}

bool CellGraph3::Edge::isValid() const { return _isValid; }

void CellGraph3::Edge::destroy() {
  _isValid = false;
}

std::vector<size_t> CellGraph3::Edge::nodes() const {
  return {_nodeIds[0], _nodeIds[1]};
}
Point3d CellGraph3::Edge::center() const { return _center; }

size_t CellGraph3::Edge::id() const {
  return _id;
}

} // furoo namespace
