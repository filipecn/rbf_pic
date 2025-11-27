template <int D> size_t Mesh<D>::elementCount() const {
  return _elements.size();
}

template <int D> size_t Mesh<D>::vertexCount() const {
  return _vertices.size();
}

template <int D> Point<double, D> Mesh<D>::vertex(size_t i) const {
  THROW(i < _vertices.size(), "Mesh<>::vertex invalid index");
  return _vertices[i];
}

template <int D> Vector<int, D> Mesh<D>::element(size_t i) const {
  THROW(i < _elements.size(), "Mesh<>::element invalid index");
  return _elements[i];
}

template <int D> BBox<double, D> Mesh<D>::elementBBox(size_t i) const {
  BBox<double, D> box(_vertices[_elements[i][0]], _vertices[_elements[i][1]]);
  if (D == 2)
    return box;
  return BBox<double, D>::combine(box, _vertices[_elements[i][2]]);
}

template <>
inline bool Mesh<2>::intersectElement(size_t elementId, Point<double, 2> origin,
                                      Vector<double, 2> direction,
                                      double *hit) const {
  return GeometricQueries::segmentRayIntersection(
      _vertices[_elements[elementId][0]], _vertices[_elements[elementId][1]],
      origin, direction, hit);
}

template <>
inline bool Mesh<3>::intersectElement(size_t elementId, Point<double, 3> origin,
                                      Vector<double, 3> direction,
                                      double *hit) const {
  return GeometricQueries::triangleRayIntersection(
      _vertices[_elements[elementId][0]], _vertices[_elements[elementId][1]],
      _vertices[_elements[elementId][2]], origin, direction, hit);
}

template <int D> void Mesh<D>::addVertex(Point<double, D> p) {
  _vertices.emplace_back(p);
}

template <int D> void Mesh<D>::addElement(Vector<int, D> e) {
  _elements.emplace_back(e);
}
