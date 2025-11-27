#include <geometry/bbox.h>

namespace furoo {

BBox2D::BBox2D() {
  _pMin = Point2d({INFINITY, INFINITY});
  _pMax = Point2d({-INFINITY, -INFINITY});
}

BBox2D::BBox2D(const Point2d &p1, const Point2d &p2) {
  _pMin = Point2d(std::min(p1[0], p2[0]), std::min(p1[1], p2[1]));
  _pMax = Point2d(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]));
}

BBox2D::BBox2D(const Point2d &c, double x, double y) {
  _pMin = c - Vector2d(x, (y <= 0) ? x : y);
  _pMax = c + Vector2d(x, (y <= 0) ? x : y);
}

BBox2D::BBox2D(const std::vector<Point2d> &points) : BBox2D() {
  for (auto p : points)
    (*this) = make_union((*this), p);
}

bool BBox2D::inside(const Point2d &p) const {
  return (p.x() >= _pMin.x() && p.x() <= _pMax.x() && p.y() >= _pMin.y() &&
          p.y() <= _pMax.y());
}

double BBox2D::size(int d) const {
  d = std::max(0, std::min(1, d));
  return _pMax[d] - _pMin[d];
}

Vector2d BBox2D::extends() const { return _pMax - _pMin; }

Point2d BBox2D::center() const { return _pMin + (_pMax - _pMin) * .5; }

int BBox2D::maxExtent() const {
  Vector2d diag = _pMax - _pMin;
  if (diag.x() > diag.y())
    return 0;
  return 1;
}

const Point2d &BBox2D::operator[](int i) const {
  return (i == 0) ? _pMin : _pMax;
}

Point2d &BBox2D::operator[](int i) { return (i == 0) ? _pMin : _pMax; }

Point2d BBox2D::lower() const { return _pMin; }

Point2d BBox2D::upper() const { return _pMax; }

BBox2D BBox2D::make_union(const BBox2D &b, const Point2d &p) {
  BBox2D ret = b;
  ret[0][0] = std::min(b.lower().x(), p.x());
  ret[0][1] = std::min(b.lower().y(), p.y());
  ret[1][0] = std::max(b.upper().x(), p.x());
  ret[1][1] = std::max(b.upper().y(), p.y());
  return ret;
}

BBox2D BBox2D::make_union(const BBox2D &a, const BBox2D &b) {
  BBox2D ret = make_union(a, b.lower());
  return make_union(ret, b.upper());
}

BBox2D BBox2D::make_unit_bbox() { return BBox2D(Point2d(), Point2d({1, 1})); }

bool BBox2D::intersect(Point2d initialPoint, Vector2d direction,
                       Point2d &intersectionPoint, Vector2d &normal) const {
  // TODO Genealize for any BBox
  // Vector2d direction = finalPoint - initialPoint;
  double x, y;
  if (direction[1] == 0) {
    x = direction[0] > 0 ? 1. + _pMax.x() : _pMin.x() - 1.;
  } else if (direction[1] > 0) { // TOP SIDE
    x = initialPoint[0] +
        direction[0] * (_pMax.y() - initialPoint[1]) / direction[1];
  } else { // BOTTOM SIDE
    x = initialPoint[0] +
        direction[0] * (_pMin.y() - initialPoint[1]) / direction[1];
  }
  if (x < _pMin.x()) {
    intersectionPoint[0] = _pMin.x();
    normal[0] = 1;
  } else if (x > _pMax.x()) {
    intersectionPoint[0] = _pMax.x();
    normal[0] = -1;
  } else {
    intersectionPoint[0] = x;
    normal[0] = 0;
  }

  if (direction[0] == 0) {
    y = direction[1] > 0 ? 1. + _pMax.y() : _pMin.y() - .1;
  } else if (direction[0] > 0) { // RIGHT SIDE
    y = (initialPoint[1]) +
        direction[1] * (_pMax.x() - initialPoint[0]) / direction[0];
  } else { // LEFT SIDE
    y = (initialPoint[1]) +
        direction[1] / direction[0] * (_pMin.x() - initialPoint[0]);
  }
  if (y < _pMin.y()) {
    intersectionPoint[1] = _pMin.y();
    normal[1] = 1;
  } else if (y > _pMax.y()) {
    intersectionPoint[1] = _pMax.y();
    normal[1] = -1;
  } else {
    intersectionPoint[1] = y;
    normal[1] = 0;
  }

  if (normal == Vector2d(0., 0.))
    normal -= direction;
  return true;
}

Point2d BBox2D::pointCoordinates(const Point2d& p) const {
  auto d = p - _pMin;
  auto e = _pMax - _pMin;
  return {d.x() / e.x(), d.y() / e.y()};
}

bool BBox2D::inside(const BBox2D& b) const {
  return inside(b._pMin) && inside(b._pMax);
}

BBox3D::BBox3D() {
  _pMin = Point3d(INFINITY);
  _pMax = Point3d(-INFINITY);
}

BBox3D::BBox3D(const Point3d &p1, const Point3d &p2) {
  _pMin = Point3d(std::min(p1[0], p2[0]), std::min(p1[1], p2[1]),
                  std::min(p1[2], p2[2]));
  _pMax = Point3d(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]),
                  std::max(p1[2], p2[2]));
}

BBox3D::BBox3D(const Point3d &c, double x, double y, double z) {
  _pMin = c - Vector3d(x, (y >= 0) ? y : x, (z >= 0) ? z : x);
  _pMax = c + Vector3d(x, (y >= 0) ? y : x, (z >= 0) ? z : x);
}

BBox3D::BBox3D(const std::vector<Point3d> &points) : BBox3D() {
  for (auto p : points)
    (*this) = make_union((*this), p);
}

bool BBox3D::inside(const Point3d &p) const {
  return (p.x() >= _pMin.x() && p.x() <= _pMax.x() && p.y() >= _pMin.y() &&
          p.y() <= _pMax.y() && p.z() >= _pMin.z() && p.z() <= _pMax.z());
}

double BBox3D::size(int d) const {
  d = std::max(0, std::min(2, d));
  return _pMax[d] - _pMin[d];
}

Vector3d BBox3D::extends() const { return _pMax - _pMin; }

Point3d BBox3D::center() const { return _pMin + (_pMax - _pMin) * .5; }

int BBox3D::maxExtent() const {
  Vector3d diag = _pMax - _pMin;
  if (diag.x() > diag.y() && diag.x() > diag.z())
    return 0;
  else if (diag.y() > diag.x() && diag.y() > diag.z())
    return 1;
  return 2;
}

const Point3d &BBox3D::operator[](int i) const {
  return (i == 0) ? _pMin : _pMax;
}

Point3d &BBox3D::operator[](int i) { return (i == 0) ? _pMin : _pMax; }

Point3d BBox3D::lower() const { return _pMin; }

Point3d BBox3D::upper() const { return _pMax; }

BBox3D BBox3D::make_union(const BBox3D &b, const Point3d &p) {
  BBox3D ret = b;
  ret[0][0] = std::min(b.lower().x(), p.x());
  ret[0][1] = std::min(b.lower().y(), p.y());
  ret[0][2] = std::min(b.lower().z(), p.z());
  ret[1][0] = std::max(b.upper().x(), p.x());
  ret[1][1] = std::max(b.upper().y(), p.y());
  ret[1][2] = std::max(b.upper().z(), p.z());
  return ret;
}

BBox3D BBox3D::make_union(const BBox3D &a, const BBox3D &b) {
  BBox3D ret = make_union(a, b.lower());
  return make_union(ret, b.upper());
}

BBox3D BBox3D::make_unit_bbox() { return BBox3D(Point3d(), Point3d(1.)); }

bool BBox3D::intersect(Point3d initialPoint, Vector3d direction,
                       Point3d &intersectionPoint, Vector3d &normal) const {
  NOT_IMPLEMENTED();
  UNUSED_VARIABLE(initialPoint);
  UNUSED_VARIABLE(direction);
  UNUSED_VARIABLE(intersectionPoint);
  UNUSED_VARIABLE(normal);
  return false;
}

} // namespace furoo
