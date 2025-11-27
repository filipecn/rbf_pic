template <typename T, int D> BBox<T, D>::BBox() = default;

template <typename T, int D>
BBox<T, D>::BBox(const Point<T, D> &p1, const Point<T, D> &p2) {
  for (int i = 0; i < D; i++) {
    _lower[i] = std::min(p1[i], p2[i]);
    _upper[i] = std::max(p1[i], p2[i]);
  }
}

template<typename T, int D>
BBox<T,D>::BBox(const Point<T, D> &p1, double radius){
  Vector<T,D> vRadius(radius);
  _lower = p1 - vRadius;
  _upper = p1 + vRadius;
}


template <typename T, int D> BBox<T, D> BBox<T, D>::squareBox(T length) {
  return BBox<T, D>(Point<T, D>(T(0)), Point<T, D>(T(length)));
}

template <typename T, int D>
std::vector<BBox<T, D>> BBox<T, D>::halfSplit(const BBox<T, D> &bbox) {
  auto low = bbox.lower();
  auto mid = bbox.center();
  auto upp = bbox.upper();
  std::vector<BBox<T, D>> children;
  /// The children are enumerated in Z-order, like
  ///           y
  ///  2  3     |
  ///  0  1      --x
  for (int childIndex = 0; childIndex < (1 << D); childIndex++) {
    Point<T, D> lower, upper;
    /// For each new child, we decide in which portion of the new cell
    /// it belongs to: lower until mid, or mid until top
    for (int dimension = 0; dimension < D; dimension++) {
      bool firstHalf = childIndex & (1 << dimension);
      if (!firstHalf) {
        lower[dimension] = low[dimension];
        upper[dimension] = mid[dimension];
      } else {
        lower[dimension] = mid[dimension];
        upper[dimension] = upp[dimension];
      }
    }
    children.emplace_back(lower, upper);
  }
  return children;
}

template <typename T, int D> Point<T, D> BBox<T, D>::lower() const {
  return _lower;
}

template <typename T, int D> Point<T, D> BBox<T, D>::upper() const {
  return _upper;
}

template <typename T, int D> Point<T, D> BBox<T, D>::center() const {
  Point<T, D> center;
  for (int i = 0; i < D; i++)
    center[i] = (_lower[i] + _upper[i]) * 0.5;
  return center;
}
template <typename T, int D>

bool BBox<T, D>::contains(Point<T, D> p) const {
  for (int i = 0; i < D; i++)
    if (p[i] < _lower[i] || p[i] > _upper[i])
      return false;
  return true;
}

template <typename T, int D>
BBox<T, D> BBox<T, D>::combine(const BBox<T, D> &a, const BBox<T, D> &b) {
  auto box = combine(a, b.lower());
  box = combine(box, b.upper());
  return box;
}

template <typename T, int D>
BBox<T, D> BBox<T, D>::combine(const BBox<T, D> &a, const Point<T, D> &b) {
  Point<T, D> lower = a.lower(), upper = a.upper();
  for (int i = 0; i < D; i++) {
    lower[i] = std::min(lower[i], b[i]);
    upper[i] = std::max(upper[i], b[i]);
  }
  return BBox<T, D>(lower, upper);
}

template <typename T, int D> bool BBox<T, D>::contains(BBox<T, D> b) const {
  return contains(b.lower()) && contains(b.upper());
}

template <typename T, int D> bool BBox<T, D>::intersect(BBox<T, D> b) const {
  for (int i = 0; i < D; i++)
    if (lower()[i] > b.upper()[i] || upper()[i] < b.lower()[i])
      return false;
  return true;
}

template <typename T, int D>
BBox<T, D> BBox<T, D>::expanded(Vector<double, D> expansionRatio) const {
  auto s = expansionRatio * size();
  auto c = center();
  return BBox<T, D>(c - s, c + s);
}

template <typename T, int D> Vector<T, D> BBox<T, D>::size() const {
  return _upper - _lower;
}

template <typename T, int D>
bool BBox<T, D>::intersect(Point<T, D> initialPoint, Vector<T, D> direction,
                           double *hit1, double *hit2) const {
  Vector<double, D> d = (direction).normalized();
  Point<double, D> p = initialPoint;
  double tmin = 0.;
  double tmax = INFINITY;
  for (int i = 0; i < D; i++) {
    if (std::fabs(d[i]) < 1e-8) {
      if (p[i] < lower()[i] || p[i] > upper()[i])
        return false;
    } else {
      double ood = 1. / d[i];
      double t1 = (lower()[i] - p[i]) * ood;
      double t2 = (upper()[i] - p[i]) * ood;
      if (t1 > t2)
        std::swap(t1, t2);
      if (t1 > tmin)
        tmin = t1;
      if (t2 > tmax)
        tmax = t2;
      if (tmin > tmax)
        return false;
    }
  }
  if (hit1)
    (*hit1) = tmin;
  if (hit2)
    (*hit2) = tmax;
  return true;
}

template <typename T, int D> Point<T,D> BBox<T, D>::clamp(Point<T, D> point) {
  Point<T, D> newPoint;
  for (int d = 0; d < D; d++) {
    if (point[d] < _lower[d])
      newPoint[d] = _lower[d];
    else if (point[d] > _upper[d])
      newPoint[d] = _upper[d];
    else
      newPoint[d] = point[d];
  }
  return newPoint;
}