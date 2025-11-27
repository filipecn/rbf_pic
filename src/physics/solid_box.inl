template <>
inline SolidBox<2>::SolidBox(const BBox<double, 2> &box) : _box(box) {
  _toUnitBox = Transform<2>::toUnitBox(_box);
  //    2
  //  3----2
  // 3 |    | 1
  //  0___ 1
  //    0
  _mesh.addVertex(Point2d(_box.lower().x(), _box.lower().y()));
  _mesh.addVertex(Point2d(_box.upper().x(), _box.lower().y()));
  _mesh.addVertex(Point2d(_box.upper().x(), _box.upper().y()));
  _mesh.addVertex(Point2d(_box.lower().x(), _box.upper().y()));
  for (int i = 0; i < 4; i++)
    _mesh.addElement(Vector2i(i, (i + 1) % 4));
}

template <>
inline SolidBox<3>::SolidBox(const BBox<double, 3> &box) : _box(box) {
  _toUnitBox = Transform<3>::toUnitBox(_box);
  _mesh.addVertex( // 0
      Point3d(_box.lower().x(), _box.lower().y(), _box.lower().z()));
  _mesh.addVertex( // 1
      Point3d(_box.upper().x(), _box.lower().y(), _box.lower().z()));
  _mesh.addVertex( // 2
      Point3d(_box.upper().x(), _box.upper().y(), _box.lower().z()));
  _mesh.addVertex( // 3
      Point3d(_box.lower().x(), _box.upper().y(), _box.lower().z()));
  _mesh.addVertex( // 4
      Point3d(_box.lower().x(), _box.lower().y(), _box.upper().z()));
  _mesh.addVertex( // 5
      Point3d(_box.upper().x(), _box.lower().y(), _box.upper().z()));
  _mesh.addVertex( // 6
      Point3d(_box.upper().x(), _box.upper().y(), _box.upper().z()));
  _mesh.addVertex( // 7
      Point3d(_box.lower().x(), _box.upper().y(), _box.upper().z()));
  _mesh.addElement(Vector3i(0, 1, 2));
  _mesh.addElement(Vector3i(0, 2, 3));
  _mesh.addElement(Vector3i(4, 5, 6));
  _mesh.addElement(Vector3i(4, 6, 7));
  _mesh.addElement(Vector3i(7, 6, 3));
  _mesh.addElement(Vector3i(6, 2, 3));
  _mesh.addElement(Vector3i(4, 0, 5));
  _mesh.addElement(Vector3i(0, 5, 1));
  _mesh.addElement(Vector3i(0, 3, 4));
  _mesh.addElement(Vector3i(7, 3, 4));
  _mesh.addElement(Vector3i(5, 1, 6));
  _mesh.addElement(Vector3i(6, 1, 2));
}

template <int D>
bool SolidBox<D>::contains(const Point<double, D> &point) const {
  return _box.contains(point);
}

template <int D>
Definitions::Intersection<D>
SolidBox<D>::intersectionPoint(const Point<double, D> &origin,
                               const Point<double, D> &destination) const {
  Definitions::Intersection<D> intersection{};
  Vector<double, D> d = (destination - origin);
  if (d.length2() < 1e-8) {
    intersection.isValid = false;
    return intersection;
  }
  d = d.normalized();
  intersection.isValid = true;
  Point<double, D> p = origin;
  double tmin = 0.;
  double tmax = INFINITY;
  for (int i = 0; i < D; i++) {
    if (std::fabs(d[i]) < 1e-8) {
      if (p[i] < _box.lower()[i] || p[i] > _box.upper()[i])
        THROW(false, "SolidBox<>::intersectionPoint no intersection");
    } else {
      double ood = 1. / d[i];
      double t1 = (_box.lower()[i] - p[i]) * ood;
      double t2 = (_box.upper()[i] - p[i]) * ood;
      if (t1 > t2)
        std::swap(t1, t2);
      if (t1 > tmin)
        tmin = t1;
      if (t2 > tmax)
        tmax = t2;
      if (tmin > tmax)
        THROW(false, "SolidBox<>::intersectionPoint no intersection");
    }
  }
  intersection.point = p + d * tmin;
  auto v = _toUnitBox(intersection.point) - Point<double, D>(0.5);
  double m = 0.;
  int md = 0;
  int signal = 1;
  for (int d = 0; d < D; d++)
    if (std::fabs(v[d]) > m) {
      m = std::fabs(v[d]);
      md = d;
      if (v[d] < 0)
        signal = -1;
      else
        signal = 1;
    }
  intersection.normal[md] = signal;
  return intersection;
}

template <int D> inline Mesh<D> *SolidBox<D>::mesh() { return &_mesh; }

template <int D> BBox<double, D> SolidBox<D>::box() const { return _box; }
