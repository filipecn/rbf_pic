#include "cell_graph.h"

namespace furoo {
CellGraph2::Vertex::Vertex(Point2d position, size_t id)
    : _position(position), _isValid(true), _id(id) {}

void CellGraph2::Vertex::set(Point2d position, size_t id) {
  _position = position;
  _isValid = true;
  _id = id;
}

void CellGraph2::Vertex::destroy() {
  _cellCount = 0;
  _isValid = false;
}

void CellGraph2::Vertex::removeCell() {
  THROW(_cellCount > 0,
        "CellGraph2::Vertex::decrementCell() no cell to decrement");
  _cellCount--;
};

void CellGraph2::Vertex::addCell() { _cellCount++; }

size_t CellGraph2::Vertex::cellCount() const { return _cellCount; }

Point2d CellGraph2::Vertex::position() const { return _position; }

bool CellGraph2::Vertex::isValid() const { return _isValid; }

size_t CellGraph2::Vertex::id() const {
  return _id;
}

} // furoo namespace

