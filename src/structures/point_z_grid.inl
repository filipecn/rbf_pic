template <typename T, int D> PointZGrid<T, D>::PointZGrid() {
  _end = 0;
  _lastId = 0;
  _nbits = 0;
  _maxDepth = 0;
  _maxZCode = 0;
  _needUpdate = true;
  _needZUpdate = true;
  _sizeChanged = true;
}

template <typename T, int D> PointZGrid<T, D>::~PointZGrid() = default;

template <typename T, int D>
PointZGrid<T, D>::PointZGrid(size_t maxCoordinates) : PointZGrid() {
  THROW(isPowerOf2(maxCoordinates), "maxCoordinates not power of 2!");
  _resolution = Point<T, D>(maxCoordinates);
  // TODO: solve this
  // auto h = 1.f / maxCoordinates;
  // toSet_ = scale(h, h, h);
  _maxZCode = computeIndex(_resolution);
  int n = maxCoordinates - 1;
  for (_nbits = 0; n; n >>= 1)
    _nbits++;
  _maxDepth = _nbits;
}

template <typename T, int D> void PointZGrid<T, D>::update() {
  if (_points.empty())
    return;
  if (_needZUpdate) {
    for (size_t i = 0; i < _end; i++) {
      Point<T, D> gp = _toSet(_positions[_points[i].id]);
      _points[i].zcode = computeIndex(gp);
      THROW(_points[i].zcode <= _maxZCode, "");
    }
    _needZUpdate = false;
    _needUpdate = true;
  }
  if (!_needUpdate)
    return;
  std::sort(&_points[0], &_points[0] + _end,
            [](const PointElement &a, const PointElement &b) {
              if (!a.active)
                return false;
              if (!b.active)
                return true;
              return a.zcode < b.zcode;
            });
  if (_indices.size() < _end)
    _indices.resize(_end);
  for (size_t i = 0; i < _end; i++) {
    _indices[_points[i].id] = i;
  }
  // compute new end in case some points have been deleted
  if (_end > 0) {
    while (_end > 0 && !_points[_end - 1].active)
      _end--;
  }
  _needUpdate = false;
  _sizeChanged = false;
}

template <typename T, int D> void PointZGrid<T, D>::clear() {
  _end = 0;
  _lastId = 0;
  _needUpdate = true;
  _needZUpdate = true;
  _sizeChanged = true;
  _points.clear();
  //  _sizeChanged = true;
  //  for (auto &p:_points)
  //    p.active = false;
  update();
}

template <typename T, int D> size_t PointZGrid<T, D>::size() {
  if (_sizeChanged)
    update();
  return _end;
}

template <typename T, int D> size_t PointZGrid<T, D>::add(Point<T, D> p) {
  if (_end == _points.size()) {
    _points.emplace_back();
    _points[_end].id = _lastId++;
    _sizeChanged = true;
  }
  if (_points[_end].id == _positions.size())
    _positions.emplace_back();
  _points[_end].active = true;
  _positions[_points[_end].id] = p;
  Point<T, D> sp = _toSet(p);
  THROW((sp >= Point<T, D>() && sp <= _resolution), "");
  _points[_end].zcode = computeIndex(sp);
  THROW(_points[_end].zcode <= _maxZCode, "");
  _end++;
  _needUpdate = true;
  return _points[_end - 1].id;
}

template <typename T, int D>
void PointZGrid<T, D>::setPosition(size_t i, Point<T, D> p) {
  _positions[i] = p;
  _needZUpdate = true;
  _needUpdate = true;
}

template <typename T, int D> void PointZGrid<T, D>::remove(size_t i) {
  THROW(i < _indices.size(), "");
  THROW(_indices[i] < _points.size(), "");
  _points[_indices[i]].active = false;
  _needUpdate = true;
  _sizeChanged = true;
}

template <typename T, int D>
Point<T, D> PointZGrid<T, D>::operator[](size_t i) const {
  THROW(i < _positions.size(), "");
  return _positions[i];
}

//#ifdef _USE_OPENMP
template <typename T, int D>
void PointZGrid<T, D>::search(const BBox<T, D> &b,
                              const std::function<void(size_t)> &f) {
  if (_needUpdate || _needZUpdate)
    update();
  // perform an implicit search
  std::function<void(size_t, size_t, const BBox<T, D> &, const BBox<T, D> &,
                     const std::function<void(size_t)> &)>
      implicitTraverse = [&](size_t level, size_t zcode,
                             const BBox<T, D> &region, const BBox<T, D> &bbox,
                             const std::function<void(size_t)> &callBack) {
        if (level > _maxDepth)
          return;
        auto gridBBox = _toSet(bbox);
        // check if node is fully contained by bbox or node is leaf, refine
        // otherwise
        if (gridBBox.contains(region) || level == _maxDepth) {
          if (gridBBox.intersect(region))
            for (Iterator it(*this, zcode, level); it.next(); ++it)
              if (bbox.contains(it.getWorldPosition()) &&
                  _points[_indices[it.getId()]].active)
                f(it.getId());
          return;
        }
        size_t d = (_nbits - level - 1) * D;
        auto regions = BBox<T, D>::halfSplit(region);
        for (size_t i = 0; i < regions.size(); i++)
          implicitTraverse(level + 1, zcode | (i << d), regions[i], bbox,
                           callBack);
      };
  implicitTraverse(0, 0, BBox<T, D>(Point<T, D>(), _resolution), b, f);
}

template <typename T, int D>
void PointZGrid<T, D>::searchZ(const BBox<T, D> &b, size_t level,
                               const std::function<void(size_t)> &f) {
  if (_needUpdate || _needZUpdate)
    update();
  Point<T, D> sp = _toSet(b.lower());
  level = std::min(level, _maxDepth);
  size_t zcode = computeIndex(sp);
  // size_t count = 0, countTrue = 0;
  for (Iterator it(*this, zcode, level); it.next(); ++it) {
    // count++;
    if (b.contains(it.getWorldPosition()) &&
        _points[_indices[it.getId()]].active) {
      f(it.getId());
      // countTrue++;
    }
  }
  // std::cerr << "count(" << count << ") true(" << countTrue << ")" <<
  // std::endl;
}

template <typename T, int D>
void PointZGrid<T, D>::iteratePoints(
  const std::function<void(size_t, Point<T, D>)> &f) const {
  for (size_t i = 0; i < _end; i++)
    if (_points[i].active)
      f(_points[i].id, _positions[_points[i].id]);
}

#ifdef _USE_OPENMP
template <typename T, int D>
void PointZGrid<T, D>::iteratePoints_par(
  const std::function<void(size_t, Point<T, D>)> &f) const
{
  using index_type = long long; // OpenMP requires signed integral type

#pragma omp parallel for
  for (index_type i = 0; i < (index_type)_end; ++i)
    if (_points[i].active)
      f(_points[i].id, _positions[_points[i].id]);
}
#endif // _USE_OPENMP

template <typename T, int D>
size_t PointZGrid<T, D>::computeIndex(Point<T, D> p) {
  if (D == 2)
    return mortonCode(p[0], p[1]);
  return mortonCode(p[0], p[1], p[2]);
}

template <>
inline void PointZGrid<double, 2>::iterateClosestPoints(
    Point<double, 2> center, size_t n,
    const std::function<void(size_t, Point<double, 2>)> &f,
    std::function<bool(size_t)> isValid) {
  auto gp = _toSet(center);
  THROW((gp >= Point2d(0.) && gp <= _resolution), "query outside domain!");
  if (_needUpdate || _needZUpdate)
    update();
  Point2i centerCell;
  for (int d = 0; d < 2; d++)
    centerCell[d] = gp[d];
  // define distance between two points
  std::function<bool(size_t, size_t)> cmp = [&](size_t left, size_t right) {
    return distance2(center, _positions[left]) <
           distance2(center, _positions[right]);
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> maxHeap(cmp);
  size_t maxRing = _resolution[0];
  for (int d = 1; d < 2; d++)
    maxRing = std::max(maxRing, static_cast<size_t>(_resolution[d]));
  // comparing function between two points
  std::function<int(const PointElement &, const size_t &)> comp =
      [](const PointElement &p, const size_t &v) {
        if (p.zcode < v)
          return -1;
        if (p.zcode > v)
          return 1;
        return 0;
      };
  /// process points from a single cell
  auto process = [&](Point2i cell) {
    Point2d cellPoint;
    for (int i = 0; i < 2; i++)
      if (cell[i] < 0 || cell[i] > static_cast<int>(_resolution[i]))
        return;
      else
        cellPoint[i] = cell[i];
    auto zcode = computeIndex(cellPoint);
    size_t cur = static_cast<size_t>(
        lower_bound<PointElement, size_t>(&_points[0], _end, zcode, comp) + 1);
    for (; cur < _end && _points[cur].zcode == zcode; cur++) {
      if (!_points[cur].active)
        continue;
      if (isValid && !isValid(_points[cur].id))
        continue;
      if (maxHeap.size() < n) {
        maxHeap.push(_points[cur].id);
      } else {
        size_t top = maxHeap.top();
        if (distance2(center, _positions[_points[cur].id]) <
            distance2(center, _positions[top])) {
          maxHeap.pop();
          maxHeap.push(_points[cur].id);
        }
      }
    }
  };
  process(centerCell);
  for (size_t ring = 1; (ring < maxRing && maxHeap.size() < n) || ring < 2;
       ring++) {
    for (int i = centerCell.x() - ring;
         i <= centerCell.x() + static_cast<int>(ring); i++) {
      process(Point2i(i, centerCell.y() - ring));
      process(Point2i(i, centerCell.y() + ring));
    }
    for (int j = centerCell.y() - ring + 1;
         j < centerCell.y() + static_cast<int>(ring); j++) {
      process(Point2i(centerCell.x() - ring, j));
      process(Point2i(centerCell.x() + ring, j));
    }
  }
  while (!maxHeap.empty()) {
    f(maxHeap.top(), _positions[maxHeap.top()]);
    maxHeap.pop();
  }
}

template <>
inline void PointZGrid<double, 3>::iterateClosestPoints(
    Point<double, 3> center, size_t n,
    const std::function<void(size_t, Point<double, 3>)> &f,
    std::function<bool(size_t)> isValid) {
  auto gp = _toSet(center);
  THROW((gp >= Point3d(0.) && gp <= _resolution), "query outside domain!");
  if (_needUpdate || _needZUpdate)
    update();
  Point3i centerCell;
  for (int d = 0; d < 3; d++)
    centerCell[d] = gp[d];
  // define distance between two points
  std::function<bool(size_t, size_t)> cmp = [&](size_t left, size_t right) {
    return distance2(center, _positions[left]) <
           distance2(center, _positions[right]);
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> maxHeap(cmp);
  size_t maxRing = _resolution[0];
  for (int d = 1; d < 3; d++)
    maxRing = std::max(maxRing, static_cast<size_t>(_resolution[d]));
  // comparing function between two points
  std::function<int(const PointElement &, const size_t &)> comp =
      [](const PointElement &p, const size_t &v) {
        if (p.zcode < v)
          return -1;
        if (p.zcode > v)
          return 1;
        return 0;
      };
  /// process points from a single cell
  auto process = [&](Point3i cell) {
    Point3d cellPoint;
    for (int i = 0; i < 3; i++)
      if (cell[i] < 0 || cell[i] > static_cast<int>(_resolution[i]))
        return;
      else
        cellPoint[i] = cell[i];
    auto zcode = computeIndex(cellPoint);
    size_t cur = static_cast<size_t>(
        lower_bound<PointElement, size_t>(&_points[0], _end, zcode, comp) + 1);
    for (; cur < _end && _points[cur].zcode == zcode; cur++) {
      if (!_points[cur].active)
        continue;
      if (isValid && !isValid(_points[cur].id))
        continue;
      if (maxHeap.size() < n) {
        maxHeap.push(_points[cur].id);
      } else {
        size_t top = maxHeap.top();
        if (distance2(center, _positions[_points[cur].id]) <
            distance2(center, _positions[top])) {
          maxHeap.pop();
          maxHeap.push(_points[cur].id);
        }
      }
    }
  };
  process(centerCell);
  for (size_t ring = 1; (ring < maxRing && maxHeap.size() < n) || ring < 2;
       ring++) {
    for (int i = centerCell.x() - ring;
         i <= centerCell.x() + static_cast<int>(ring); i++) {
      process(Point3i(i, centerCell.y() - ring, centerCell.z() - ring));
      process(Point3i(i, centerCell.y() + ring, centerCell.z() - ring));
      process(Point3i(i, centerCell.y() - ring, centerCell.z() + ring));
      process(Point3i(i, centerCell.y() + ring, centerCell.z() + ring));
    }
    for (int j = centerCell.y() - ring + 1;
         j < centerCell.y() + static_cast<int>(ring); j++) {
      process(Point3i(centerCell.x() - ring, j, centerCell.z() - ring));
      process(Point3i(centerCell.x() + ring, j, centerCell.z() - ring));
      process(Point3i(centerCell.x() - ring, j, centerCell.z() + ring));
      process(Point3i(centerCell.x() + ring, j, centerCell.z() + ring));
    }
    for (int k = centerCell.z() - ring + 1;
         k < centerCell.z() + static_cast<int>(ring); k++) {
      process(Point3i(centerCell.x() - ring, centerCell.y() - ring, k));
      process(Point3i(centerCell.x() + ring, centerCell.y() - ring, k));
      process(Point3i(centerCell.x() - ring, centerCell.y() + ring, k));
      process(Point3i(centerCell.x() + ring, centerCell.y() + ring, k));
    }
  }
  while (!maxHeap.empty()) {
    f(maxHeap.top(), _positions[maxHeap.top()]);
    maxHeap.pop();
  }
}

template <typename T, int D>
void PointZGrid<T, D>::setDomainRegion(const BBox<T, D> &region) {
  Vector<double, D> scaleVector, regionSize = region.size();
  for (int i = 0; i < D; i++)
    scaleVector[i] = _resolution[i] / regionSize[i];
  _toSet = Transform<D>::scale(scaleVector);
}

template <typename T, int D>
bool PointZGrid<T, D>::containsPoint(const BBox<T, D> &bbox) {
  if (_needUpdate || _needZUpdate)
    update();
  // perform an implicit search
  bool found = false;
  auto gridBBox = _toSet(bbox);
  std::function<void(size_t, size_t, const BBox<T, D> &)> implictTraverse =
      [&](size_t level, size_t zcode, const BBox<T, D> &region) {
        if (found || level > _maxDepth)
          return;
        // check if node is fully contained by bbox or node is leaf, refine
        // otherwise
        if (gridBBox.contains(region) || level == _maxDepth) {
          if (gridBBox.intersect(region))
            for (Iterator it(*this, zcode, level); it.next() && !found; ++it)
              if (bbox.contains(it.getWorldPosition()) &&
                  _points[_indices[it.getId()]].active)
                found = true;
          return;
        }
        size_t d = (_nbits - level - 1) * D;
        auto regions = BBox<T, D>::halfSplit(region);
        for (size_t i = 0; i < regions.size() && !found; i++)
          implictTraverse(level + 1, zcode | (i << d), regions[i]);
      };
  implictTraverse(0, 0, BBox<T, D>(Point<T, D>(), _resolution));
  return found;
}

template <typename T, int D>
bool PointZGrid<T, D>::containsPointZ(const BBox<T, D> &b, size_t level) {
  if (_needUpdate || _needZUpdate)
    update();
  Point<T, D> sp = _toSet(b.lower());
  level = std::min(level, _maxDepth);
  size_t zcode = computeIndex(sp);
  for (Iterator it(*this, zcode, level); it.next(); ++it) {
    if (b.contains(it.getWorldPosition()) &&
        _points[_indices[it.getId()]].active)
      return true;
  }
  return false;
}
