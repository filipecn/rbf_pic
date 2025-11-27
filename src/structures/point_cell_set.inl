template<typename T, int D>
PointCellSet<T, D>::PointCellSet() {
  _end = 0;
  _lastId = 0;
  _needsUpdate = true;
}

template<typename T, int D>
PointCellSet<T, D>::~PointCellSet() = default;

template<typename T, int D>
void PointCellSet<T, D>::update() {
  if (_points.empty())
    return;
  if (!_needsUpdate)
    return;
  std::sort(&_points[0], &_points[0] + _end,
            [](const PointElement &a, const PointElement &b) {
              if (!a.active)
                return false;
              if (!b.active)
                return true;
              return a.cellId < b.cellId;
            });
  if (_indices.size() < _end)
    _indices.resize(_end);
  for (size_t i = 0; i < _end; i++)
    _indices[_points[i].id] = i;
  // compute new end in case some points have been deleted
  if (_end > 0) {
    while (_end > 0 && !_points[_end - 1].active)
      _end--;
  }
  _needsUpdate = false;
}

template<typename T, int D>
size_t PointCellSet<T, D>::size() {
  return _end;
}

template<typename T, int D>
size_t PointCellSet<T, D>::add(Point <T, D> p) {
  if (_end == _points.size()) {
    _points.emplace_back();
    _points[_end].id = _lastId++;
  }
  if (_points[_end].id == _positions.size())
    _positions.emplace_back();
  _points[_end].active = true;
  _positions[_points[_end].id] = p;
  _end++;
  _needsUpdate = true;
  return _points[_end - 1].id;
}

template<typename T, int D>
void PointCellSet<T, D>::setPosition(size_t i, Point <T, D> p) {
  _positions[i] = p;
}

template<typename T, int D>
void PointCellSet<T, D>::remove(size_t i) {
  THROW(i < _indices.size(), "");
  THROW(_indices[i] < _points.size(), "");
  _points[_indices[i]].active = false;
  _needsUpdate = true;
}

template<typename T, int D>
Point <T, D> PointCellSet<T, D>::operator[](size_t i) const {
  THROW(i < _positions.size(), "");
  return _positions[i];
}

template<typename T, int D>
void PointCellSet<T, D>::search(const BBox <T, D> &b,
                                const std::function<void(size_t)> &f) {
  UNUSED_VARIABLE(b);
  UNUSED_VARIABLE(f);
  NOT_IMPLEMENTED();
}

template<typename T, int D>
void PointCellSet<T, D>::iteratePoints(const std::function<void(size_t,
                                                                Point <T,
                                                                D
                                                                >)> &f) const {
  for (size_t i = 0; i < _end; i++)
    if (_points[i].active)
      f(_points[i].id, _positions[_points[i].id]);
}

template<typename T, int D>
void PointCellSet<T, D>::iterateClosestPoints(Point <T, D> center,
                                              size_t n,
                                              const std::function<void(size_t,
                                                                       Point
                                                                       <T,
                                                                       D
                                                                       >)> &f,
                                              std::function<bool(size_t)> isValid) {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(center);
  UNUSED_VARIABLE(n);
  UNUSED_VARIABLE(f);
  UNUSED_VARIABLE(isValid);
}

template<typename T, int D>
void PointCellSet<T, D>::setDomainRegion(const BBox <T, D> &region) {
  UNUSED_VARIABLE(region);
  NOT_IMPLEMENTED();
}

template<typename T, int D>
void PointCellSet<T, D>::setCell(size_t pointId, size_t cellId) {
  if (_points[_indices[pointId]].cellId != cellId)
    _needsUpdate = true;
  _points[_indices[pointId]].cellId = cellId;
}

template<typename T, int D>
void PointCellSet<T, D>::getCell(size_t pointId) const {
  return _points[_indices[pointId]].cellId;
}

template<typename T, int D>
void PointCellSet<T, D>::iterateCellPoints(size_t cellId,
                                           const std::function<void(size_t,
                                                                    Point <T, D>)> &f) {

}
