#include "cell_graph.h"

namespace furoo {

CellGraph2::Node::Node(const size_t *vertices, const BBox2d &r, size_t id,
                       size_t level)
    : _level(level), _id(id), _isValid(true) {
  set(vertices, r, level);
}

void CellGraph2::Node::set(const size_t *vertices, const BBox2d &r,
                           size_t level) {
  for (size_t i = 0; i < 4u; i++)
    _cornerVertices[i] = vertices[i];
  _region = r;
  _level = level;
  _isValid = true;
  _edges.clear();
}

BBox2d CellGraph2::Node::region() const { return _region; }

size_t CellGraph2::Node::level() const { return _level; }

size_t CellGraph2::Node::vertex(size_t i) const { return _cornerVertices[i]; }

const std::set<size_t> &CellGraph2::Node::edges() const { return _edges; }

std::vector<size_t> CellGraph2::Node::vertices() const {
  std::vector<size_t> cornerVertices;
  for (int i = 0; i < 4; i++)
    cornerVertices.emplace_back(_cornerVertices[i]);
  return cornerVertices;
}

std::set<size_t> CellGraph2::Node::edgesCopy() const { return _edges; }

void CellGraph2::Node::addEdge(size_t e) { _edges.insert(e); }

void CellGraph2::Node::removeEdge(size_t e) {
  THROW(_edges.find(e) != _edges.end(),
        "CellGraph2::Node::removeEdge() edge not found!");
  _edges.erase(e);
}

bool CellGraph2::Node::isValid() const { return _isValid; }

size_t CellGraph2::Node::id() const { return _id; }

void CellGraph2::Node::destroy() {
  _isValid = false;
  _edges.clear();
}

Point2u CellGraph2::Node::addressCode() const { return _addressCode; }

void CellGraph2::Node::setAddressCode(Point2u addressCode) {
  _addressCode = addressCode;
}

} // namespace furoo
