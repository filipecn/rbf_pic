#include "cell_graph.h"

namespace furoo {

CellGraph3::Node::Node(const size_t *vertices, const BBox3d &r, size_t id,
                       size_t level)
    : _level(level), _id(id), _isValid(true) {
  set(vertices, r, level);
}

void CellGraph3::Node::set(const size_t *vertices, const BBox3d &r,
                           size_t level) {
  for (size_t i = 0; i < 8u; i++)
    _cornerVertices[i] = vertices[i];
  _region = r;
  _level = level;
  _isValid = true;
  _edges.clear();
}

BBox3d CellGraph3::Node::region() const { return _region; }

size_t CellGraph3::Node::level() const { return _level; }

size_t CellGraph3::Node::vertexAt(size_t i) const { return _cornerVertices[i]; }

const std::set<size_t> &CellGraph3::Node::edges() const { return _edges; }

std::vector<size_t> CellGraph3::Node::vertices() const {
  std::vector<size_t> cornerVertices;
  for (int i = 0; i < 8; i++)
    cornerVertices.emplace_back(_cornerVertices[i]);
  return cornerVertices;
}

std::set<size_t> CellGraph3::Node::edgesCopy() const { return _edges; }

void CellGraph3::Node::addEdge(size_t e) { _edges.insert(e); }

void CellGraph3::Node::removeEdge(size_t e) {
  THROW(_edges.find(e) != _edges.end(),
        "CellGraph3::Node::removeEdge() edge not found!");
  _edges.erase(e);
}

bool CellGraph3::Node::isValid() const { return _isValid; }

size_t CellGraph3::Node::id() const { return _id; }

void CellGraph3::Node::destroy() {
  _isValid = false;
  _edges.clear();
}

Point3u CellGraph3::Node::addressCode() const { return _addressCode; }

void CellGraph3::Node::setAddressCode(Point3u addressCode) {
  _addressCode = addressCode;
}

Point3u CellGraph3::Node::addressCode(size_t vertexNodeId,
                                      size_t maxLevel) const {
  size_t indexStep = 1u << (maxLevel - _level);
  Vector3u delta;
  if (vertexNodeId % 2)
    delta += Vector3u(indexStep, 0u, 0u);
  if (vertexNodeId > 3)
    delta += Vector3u(0u, 0u, indexStep);
  if (vertexNodeId == 2 || vertexNodeId == 3 || vertexNodeId == 6 ||
      vertexNodeId == 7)
    delta += Vector3u(0u, indexStep, 0u);
  return _addressCode + delta;
}

} // namespace furoo