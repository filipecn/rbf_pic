#include <structures/point_z_grid.h>

namespace furoo {

PointZGrid2::PointZGrid2() {
  _end = 0;
  _lastId = 0;
  _tree = nullptr;
  _needUpdate = true;
  _needZUpdate = false;
  _sizeChanged = false;
}

PointZGrid2::~PointZGrid2() {}

size_t PointZGrid2::add(Point2d p) {
  if (_end == _points.size()) {
    _points.emplace_back();
    _points[_end].id = _lastId++;
  }
  if (_points[_end].id == _positions.size())
    _positions.emplace_back();
  _points[_end].active = true;
  _positions[_points[_end].id] = p;
  Point2d gp = _grid.transformToGrid(p);
  ASSERT_FATAL(gp >= Point2d(0., 0.) &&
      gp <= Point2d(_grid.gridSize().x(), _grid.gridSize().y()));
  _points[_end].zcode = _computeIndex(gp);
  ASSERT_FATAL(_points[_end].zcode <= _maxZCode);
  _end++;
  _needUpdate = true;
  _sizeChanged = true;
  return _points[_end - 1].id;
}

void PointZGrid2::search(const BBox2D &b,
                         const std::function<void(size_t)> &f) {
  if (/*!_tree ||*/ _needUpdate || _needZUpdate)
    update();
  if (!_tree)
    return;
  _tree->iteratePoints(b, [&](PointElement *p) { f(p->id); });
  return;
  // perform an implicit search
  std::function<void(size_t, size_t, const BBox2D &, const BBox2D &,
                     const std::function<void(size_t)> &)>
      implicitTraverse = [&](size_t level, size_t zcode, const BBox2D &region,
                             const BBox2D &bbox,
                             const std::function<void(size_t)> &callBack) {
    if (level >= _maxDepth)
      return;
    // check if node is fully contained by bbox or node is leaf, refine
    // otherwise
    if (bbox.inside(region) || level == _maxDepth - 1) {
      if (bbox_bbox_intersection(bbox, region))
        for (PointIterator it(*this, zcode, level); it.next(); ++it)
          if (bbox.inside(it.getWorldPosition()))
            f(it.getId());
      return;
    }
    size_t d = (_nbits - level - 1) * 3;
    auto regions = region.splitBy4();
    for (size_t i = 0; i < 4; i++)
      implicitTraverse(level + 1, zcode | (i << d), regions[i], bbox,
                       callBack);

  };
  implicitTraverse(
      0, 0,
      BBox2D(Point2d(), Point2d(_grid.gridSize()[0], _grid.gridSize()[1])), b,
      f);
}

void PointZGrid2::update() {
  if (!_points.size())
    return;
  if (_needZUpdate) {
    for (size_t i = 0; i < _end; i++) {
      Point2d gp = _grid.transformToGrid(_positions[_points[i].id]);
      _points[i].zcode = mortonCode(gp.x(), gp.y());
      ASSERT_FATAL(_points[i].zcode <= _maxZCode);
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
              return a.zcode < b.zcode ? true : false;
            });
  if (_indices.size() < _end)
    _indices.resize(_end);
  for (size_t i = 0; i < _end; i++) {
    _indices[_points[i].id] = i;
  }
  // compute new end in case some particles have been deleted
  if (_end > 0) {
    while (_end > 0 && !_points[_end - 1].active)
      _end--;
  }
  delete _tree;
  _tree = new SearchTree(*this);
  _needUpdate = false;
  _sizeChanged = false;
}

void PointZGrid2::set() {
  _needUpdate = true;
  _sizeChanged = true;
  size_t w = _grid.gridSize().x();
  size_t h = _grid.gridSize().y();
  ASSERT_FATAL(isPowerOf2(w) && isPowerOf2(h));
  ASSERT(w == h); // TODO TEMPORARY
  int n = std::max(w, h) - 1;
  for (_nbits = 0; n; n >>= 1) {
    _nbits++;
  }
  _maxDepth = _nbits;
  _maxZCode = mortonCode(w, h);
}

size_t PointZGrid2::size() {
  if (_sizeChanged)
    update();
  return _end;
}

Point2d PointZGrid2::operator[](size_t
                                i) const {
  ASSERT_FATAL(i < _positions.size());
  return _positions[i];
}

void PointZGrid2::setPosition(unsigned int i, Point2d p) {
  _positions[i] = p;
  _needZUpdate = true;
  _needUpdate = true;
}

void PointZGrid2::remove(unsigned int i) {
  ASSERT_FATAL(i < _indices.size());
  ASSERT_FATAL(_indices[i] < _points.size());
  _points[_indices[i]].active = false;
  _needUpdate = true;
  _sizeChanged = true;
}

void PointZGrid2::iteratePoints(
    const std::function<void(size_t, Point2d)> &f) const {
  // TODO we should update point order before iterating!!
  for (size_t i = 0; i < _end; i++)
    if (_points[i].active)
      f(_points[i].id, _positions[_points[i].id]);
}

void PointZGrid2::iterateClosestPoints(Point2d center,
                                       size_t n,
                                       const std::function<void(size_t,
                                                                Point2d)> &f,
                                       std::function<bool(size_t)> isValid) {
  Point2d gp = _grid.transformToGrid(center);
  ASSERT_FATAL(gp >= Point2d(0., 0.) &&
      gp <= Point2d(_grid.gridSize().x(), _grid.gridSize().y()));
  if (_needUpdate || _needZUpdate)
    update();
  Point2i centerCell(static_cast<int>(gp.x()), static_cast<int>(gp.y()));
  std::function<bool(size_t, size_t)> cmp = [&](size_t left, size_t right) {
    return distance2(center, _positions[left])
        < distance2(center, _positions[right]);
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> maxHeap(cmp);
  size_t maxRing = std::max(_grid.gridSize().x(), _grid.gridSize().y());
  std::function<int(const PointElement &, const size_t &)> comp =
      [](const PointElement &p, const size_t &v) {
        if (p.zcode < v)
          return -1;
        if (p.zcode > v)
          return 1;
        return 0;
      };
  auto process = [&](int i, int j) {
    if (i < 0 || j < 0 || i >= static_cast<int>(_grid.gridSize().x()) ||
        j >= static_cast<int>(_grid.gridSize().y()))
      return;
    auto zcode = _computeIndex(Point2d(i, j));
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
        if (distance2(center, _positions[_points[cur].id])
            < distance2(center, _positions[top])) {
          maxHeap.pop();
          maxHeap.push(_points[cur].id);
        }
      }
    }
  };
  process(centerCell.x(), centerCell.y());
  for (size_t ring = 1; (ring < maxRing && maxHeap.size() < n) || ring < 2;
       ring++) {
    for (int i = centerCell.x() - ring;
         i <= centerCell.x() + static_cast<int>(ring); i++) {
      process(i, centerCell.y() - ring);
      process(i, centerCell.y() + ring);
    }
    for (int j = centerCell.y() - ring + 1;
         j < centerCell.y() + static_cast<int>(ring); j++) {
      process(centerCell.x() - ring, j);
      process(centerCell.x() + ring, j);
    }
  }
  while (!maxHeap.empty()) {
    f(maxHeap.top(), _positions[maxHeap.top()]);
    maxHeap.pop();
  }
}

} // namespace furoo
