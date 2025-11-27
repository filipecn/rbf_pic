template <int D>
std::vector<Point<double, D>>
Distribution<D>::linear(Point<double, D> a, Point<double, D> b, size_t n) {
  Vector<double, D> dir = b - a;
  double dist = dir.length();
  dir = dir.normalized();
  double step = dist / (n + 1);
  std::vector<Point<double, D>> points(1, a);
  for (size_t i = 1; i <= n; i++)
    points.emplace_back(a + i * step * dir);
  points.emplace_back(b);
  return points;
}

template <>
inline std::vector<Point<double, 2>> Distribution<2>::linear(BBox<double, 2> bbox,
                                                      Definitions::Side side,
                                                      size_t n) {
  Point2d a, b;
  switch (side) {
  case Definitions::Side::LEFT:
    a = bbox.lower();
    b[0] = bbox.lower().x();
    b[1] = bbox.upper().y();
    break;
  case Definitions::Side::RIGHT:
    b = bbox.upper();
    a[0] = bbox.upper().x();
    a[1] = bbox.lower().y();
    break;
  case Definitions::Side::BOTTOM:
    a = bbox.lower();
    b[0] = bbox.upper().x();
    b[1] = bbox.lower().y();
    break;
  case Definitions::Side::TOP:
    b = bbox.upper();
    a[0] = bbox.lower().x();
    a[1] = bbox.upper().y();
    break;
  default:
    THROW(false, "Distribution<>::linear invalid box side");
  }
  return linear(a, b, n);
}

template <>
inline std::vector<Point<double, 3>> Distribution<3>::linear(BBox<double, 3> bbox,
                                                      Definitions::Side side,
                                                      size_t n) {
  NOT_IMPLEMENTED();
  Point3d a, b;
  switch (side) {
  case Definitions::Side::LEFT:
    a = bbox.lower();
    b[0] = bbox.lower().x();
    b[1] = bbox.upper().y();
    break;
  case Definitions::Side::RIGHT:
    b = bbox.upper();
    a[0] = bbox.upper().x();
    a[1] = bbox.lower().y();
    break;
  case Definitions::Side::BOTTOM:
    a = bbox.lower();
    b[0] = bbox.upper().x();
    b[1] = bbox.lower().y();
    break;
  case Definitions::Side::TOP:
    b = bbox.upper();
    a[0] = bbox.lower().x();
    a[1] = bbox.upper().y();
    break;
  default:
    THROW(false, "Distribution<>::linear invalid box side");
  }
  return std::vector<Point3d>();
}
