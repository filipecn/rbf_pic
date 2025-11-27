#include <geometry/numeric.h>
#include <geometry/transform.h>

namespace furoo {
Transform2D::Transform2D() {
  _m.setIdentity();
  _m_inv.setIdentity();
}

Transform2D::Transform2D(const Transform2D &other) {
  _m = other._m;
  _m_inv = other._m_inv;
}

Transform2D::Transform2D(const Matrix3d &mat, const Matrix3d inv_mat)
    : _m(mat), _m_inv(inv_mat) {}

Transform2D::Transform2D(const BBox2D &bbox) {
  _m(2, 0) = _m(2, 1) = 0.0;
  _m(2, 2) = 1.0;
  _m(0, 0) = 1. / bbox.size(0);
  _m(1, 1) = 1. / bbox.size(1);
  _m(0, 2) = -(1. / bbox.size(0)) * bbox.lower().x();
  _m(1, 2) = -(1. / bbox.size(1)) * bbox.lower().y();
  _m_inv = furoo::inverse(_m);
}

void Transform2D::reset() { _m.setIdentity(); }

Transform2D Transform2D::applyTranslate(const Vector2d &d) {
  (*this) *= translate(d);
  return *this;
}

Transform2D Transform2D::applyScale(double x, double y) {
  (*this) *= scale(x, y);
  return *this;
}

Transform2D Transform2D::applyRotate(double angle) {
  (*this) *= rotate(angle);
  return *this;
}

void Transform2D::operator()(const Point2d &p, Point2d *r) const {
  double x = p[0], y = p[1];

  (*r)[0] = _m(0, 0) * x + _m(0, 1) * y + _m(0, 2);
  (*r)[1] = _m(1, 0) * x + _m(1, 1) * y + _m(1, 2);
  double wp = _m(2, 0) * x + _m(2, 1) * y + _m(2, 2);
  if (wp != 1.) {
    *r /= wp;
  }
}

void Transform2D::operator()(const Vector2d &v, Vector2d *r) const {
  double x = v.x(), y = v.y();

  (*r)[0] = _m(0, 0) * x + _m(0, 1) * y;
  (*r)[1] = _m(1, 0) * x + _m(1, 1) * y;
}

Vector2d Transform2D::operator()(const Vector2d &v) const {
  double x = v.x(), y = v.y();

  return Vector2d(_m(0, 0) * x + _m(0, 1) * y, _m(1, 0) * x + _m(1, 1) * y);
}

Point2d Transform2D::operator()(const Point2d &p) const {
  double x = p.x(), y = p.y();
  double xp = _m(0, 0) * x + _m(0, 1) * y + _m(0, 2);
  double yp = _m(1, 0) * x + _m(1, 1) * y + _m(1, 2);
  double wp = _m(2, 0) * x + _m(2, 1) * y + _m(2, 2);

  if (wp == 1.f) {
    return Point2d(xp, yp);
  }
  return Point2d(xp / wp, yp / wp);
}

Point2d Transform2D::operator()(const Point2i &p) const {
  return (*this)(Point2d(p.x(), p.y()));
}

BBox2D Transform2D::operator()(const BBox2D &b) const {
  const Transform2D &M = *this;
  BBox2D ret;

  ret = BBox2D::make_union(ret, M(Point2d(b.lower().x(), b.lower().y())));
  ret = BBox2D::make_union(ret, M(Point2d(b.upper().x(), b.lower().y())));
  ret = BBox2D::make_union(ret, M(Point2d(b.upper().x(), b.upper().y())));
  ret = BBox2D::make_union(ret, M(Point2d(b.lower().x(), b.upper().y())));
  return ret;
}

Transform2D &Transform2D::operator=(const Transform2D &t) {
  _m = t._m;
  _m_inv = t._m_inv;
  return *this;
}

Transform2D &Transform2D::operator*=(const Transform2D &t) {
  _m = _m * t._m;
  _m_inv = t._m_inv * _m_inv;
  return *this;
}

Transform2D Transform2D::operator*(const Transform2D &t) const {
  Matrix3d m1 = _m * t._m;
  Matrix3d m1_inv = t._m_inv * _m_inv;
  return Transform2D(m1, m1_inv);
}

Matrix3d Transform2D::getMatrix() const { return _m; }

Transform2D Transform2D::scale(double x, double y) {
  Matrix3d m({x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 1.0});
  return Transform2D(m, furoo::inverse(m));
}

Transform2D Transform2D::rotate(double angle) {
  double sin_a = sin(TO_RADIANS(angle));
  double cos_a = cos(TO_RADIANS(angle));
  Matrix3d m({cos_a, -sin_a, 0., sin_a, cos_a, 0., 0., 0., 1.});

  return Transform2D(m, m.transpose());
}

Transform2D Transform2D::translate(const Vector2d &v) {
  Matrix3d m({1., 0., v.x(), 0., 1., v.y(), 0., 0., 1.});
  Matrix3d m_inv({1., 0., -v.x(), 0., 1., -v.y(), 0., 0., 1.});

  return Transform2D(m, m_inv);
}

Transform2D Transform2D::inverse() const { return Transform2D(_m_inv, _m); }

Vector2d Transform2D::applyReflection(const Vector2d &v, const Vector2d &n) {
  return v - n * (2. * v.dot(n));
}

Transform3D::Transform3D() {
  _m.setIdentity();
  _m_inv.setIdentity();
}

Transform3D::Transform3D(const Transform3D &other) {
  _m = other._m;
  _m_inv = other._m_inv;
}

Transform3D::Transform3D(const Matrix4d &mat, const Matrix4d inv_mat)
    : _m(mat), _m_inv(inv_mat) {}

Transform3D::Transform3D(const BBox3D &bbox) {
  _m(3, 0) = _m(3, 1) = _m(3, 2) = 0.0;
  _m(3, 3) = 1.0;
  _m(0, 0) = 1. / bbox.size(0);
  _m(1, 1) = 1. / bbox.size(1);
  _m(2, 2) = 1. / bbox.size(2);
  _m(0, 3) = -(1. / bbox.size(0)) * bbox.lower().x();
  _m(1, 3) = -(1. / bbox.size(1)) * bbox.lower().y();
  _m(2, 3) = -(1. / bbox.size(2)) * bbox.lower().z();
  _m_inv = furoo::inverse(_m);
}

void Transform3D::reset() { _m.setIdentity(); }

Transform3D Transform3D::applyTranslate(const Vector3d &d) {
  (*this) *= translate(d);
  return *this;
}

Transform3D Transform3D::applyScale(double x, double y, double z) {
  (*this) *= scale(x, y, z);
  return *this;
}

Transform3D Transform3D::applyRotate(double angle, Vector3d axis) {
  (*this) *= rotate(angle, axis);
  return *this;
}

void Transform3D::operator()(const Point3d &p, Point3d *r) const {
  double x = p[0], y = p[1], z = p[2];

  (*r)[0] = _m(0, 0) * x + _m(0, 1) * y + _m(0, 2) * z + _m(0, 3);
  (*r)[1] = _m(1, 0) * x + _m(1, 1) * y + _m(1, 2) * z + _m(1, 3);
  (*r)[2] = _m(2, 0) * x + _m(2, 1) * y + _m(2, 2) * z + _m(2, 3);
  double wp = _m(3, 0) * x + _m(3, 1) * y + _m(3, 2) * z + _m(3, 3);
  if (wp != 1.) {
    *r /= wp;
  }
}

void Transform3D::operator()(const Vector3d &v, Vector3d *r) const {
  double x = v.x(), y = v.y(), z = v.z();

  (*r)[0] = _m(0, 0) * x + _m(0, 1) * y + _m(0, 2) * z;
  (*r)[1] = _m(1, 0) * x + _m(1, 1) * y + _m(1, 2) * z;
  (*r)[2] = _m(2, 0) * x + _m(2, 1) * y + _m(2, 2) * z;
}

Vector3d Transform3D::operator()(const Vector3d &v) const {
  double x = v.x(), y = v.y(), z = v.z();

  return Vector3d(_m(0, 0) * x + _m(0, 1) * y + _m(0, 2) * z,
                  _m(1, 0) * x + _m(1, 1) * y + _m(1, 2) * z,
                  _m(2, 0) * x + _m(2, 1) * y + _m(2, 2) * z);
}

Point3d Transform3D::operator()(const Point3d &p) const {
  double x = p.x(), y = p.y(), z = p.z();
  double xp = _m(0, 0) * x + _m(0, 1) * y + _m(0, 2) * z + _m(0, 3);
  double yp = _m(1, 0) * x + _m(1, 1) * y + _m(1, 2) * z + _m(1, 3);
  double zp = _m(2, 0) * x + _m(2, 1) * y + _m(2, 2) * z + _m(2, 3);
  double wp = _m(3, 0) * x + _m(3, 1) * y + _m(3, 2) * z + _m(3, 3);

  if (wp == 1.f) {
    return Point3d(xp, yp, zp);
  }
  return Point3d(xp / wp, yp / wp, zp / wp);
}

Point3d Transform3D::operator()(const Point3i &p) const {
  return (*this)(Point3d(p.x(), p.y(), p.z()));
}

BBox3D Transform3D::operator()(const BBox3D &b) const {
  const Transform3D &M = *this;
  BBox3D ret;
  ret = BBox3D::make_union(
      ret, M(Point3d(b.lower().x(), b.lower().y(), b.lower().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.upper().x(), b.lower().y(), b.lower().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.upper().x(), b.upper().y(), b.lower().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.lower().x(), b.upper().y(), b.lower().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.lower().x(), b.lower().y(), b.upper().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.upper().x(), b.lower().y(), b.upper().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.upper().x(), b.upper().y(), b.upper().z())));
  ret = BBox3D::make_union(
      ret, M(Point3d(b.lower().x(), b.upper().y(), b.upper().z())));
  return ret;
}

Transform3D &Transform3D::operator=(const Transform3D &t) {
  _m = t._m;
  _m_inv = t._m_inv;
  return *this;
}

Transform3D &Transform3D::operator*=(const Transform3D &t) {
  _m = _m * t._m;
  _m_inv = t._m_inv * _m_inv;
  return *this;
}

Transform3D Transform3D::operator*(const Transform3D &t) const {
  Matrix4d m1 = _m * t._m;
  Matrix4d m1_inv = t._m_inv * _m_inv;
  return Transform3D(m1, m1_inv);
}

Matrix4d Transform3D::getMatrix() const { return _m; }

Transform3D Transform3D::scale(double x, double y, double z) {
  Matrix4d m({x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0,
              0.0, 1.0});
  return Transform3D(m, furoo::inverse(m));
}

Transform3D Transform3D::rotate(double angle, Vector3d axis) {
  Vector3d a = axis.normalized();
  float s = sinf(TO_RADIANS(angle));
  float c = cosf(TO_RADIANS(angle));
  Matrix4d m;

  m(0, 0) = a.x() * a.x() + (1.f - a.x() * a.x()) * c;
  m(0, 1) = a.x() * a.y() * (1.f - c) - a.z() * s;
  m(0, 2) = a.x() * a.z() * (1.f - c) + a.y() * s;
  m(0, 3) = 0;

  m(1, 0) = a.x() * a.y() * (1.f - c) + a.z() * s;
  m(1, 1) = a.y() * a.y() + (1.f - a.y() * a.y()) * c;
  m(1, 2) = a.y() * a.z() * (1.f - c) - a.x() * s;
  m(1, 3) = 0;

  m(2, 0) = a.x() * a.z() * (1.f - c) - a.y() * s;
  m(2, 1) = a.y() * a.z() * (1.f - c) + a.x() * s;
  m(2, 2) = a.z() * a.z() + (1.f - a.z() * a.z()) * c;
  m(2, 3) = 0;

  m(3, 0) = 0;
  m(3, 1) = 0;
  m(3, 2) = 0;
  m(3, 3) = 1;

  return Transform3D(m, m.transpose());
}

Transform3D Transform3D::translate(const Vector3d &v) {
  Matrix4d m({1., 0., 0., v.x(), 0., 1., 0., v.y(), 0., 0., 1., v.z(), 0., 0.,
              0., 1.});
  Matrix4d m_inv({1., 0., 0., -v.x(), 0., 1., 0., -v.y(), 0., 0., 1., -v.z(),
                  0., 0., 0., 1.});

  return Transform3D(m, m_inv);
}

Transform3D Transform3D::inverse() const { return Transform3D(_m_inv, _m); }

Vector3d Transform3D::applyReflection(const Vector3d &v, const Vector3d &n) {
  return v - n * (2. * v.dot(n));
}

} // furoo namespace
