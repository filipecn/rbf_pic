#include "cell_graph.h"

namespace furoo {
CellGraph3::Vertex::Vertex(Point3d position, size_t id)
    : _position(position), _isValid(true), _id(id) {}

void CellGraph3::Vertex::set(Point3d position, size_t id) {
  _position = position;
  _isValid = true;
  _id = id;
}

void CellGraph3::Vertex::destroy() {
  _cellCount = 0;
  _isValid = false;
}

void CellGraph3::Vertex::removeCell() {
  THROW(_cellCount > 0,
        "CellGraph3::Vertex::decrementCell() no cell to decrement");
  _cellCount--;
};

void CellGraph3::Vertex::addCell() { _cellCount++; }

size_t CellGraph3::Vertex::cellCount() const { return _cellCount; }

Point3d CellGraph3::Vertex::position() const { return _position; }

bool CellGraph3::Vertex::isValid() const { return _isValid; }

size_t CellGraph3::Vertex::id() const {
  return _id;
}

} // furoo namespace

