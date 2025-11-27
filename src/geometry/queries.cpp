#include <geometry/matrix.h>
#include <geometry/queries.h>

namespace furoo {

bool bbox_bbox_intersection(const BBox2D &a, const BBox2D &b) {
  for (int i = 0; i < 2; i++)
    if (!((a.lower()[i] <= b.lower()[i] && a.upper()[i] >= b.lower()[i]) ||
          (b.lower()[i] <= a.lower()[i] && b.upper()[i] >= a.lower()[i])))
      return false;
  return true;
}

bool GeometricQueries::triangleRayIntersection(
    const Point3d &p1, const Point3d &p2, const Point3d &p3,
    const Point3d &origin, const Vector3d &direction, double *tHit, double *b1,
    double *b2) {
  // Compute s1
  Vector3d e1 = p2 - p1;
  Vector3d e2 = p3 - p1;
  Vector3d s1 = cross(direction, e2);
  double divisor = s1.dot(e1);
  if (divisor == 0.f)
    return false;
  double invDivisor = 1.f / divisor;
  // Compute first barycentric coordinate
  Vector3d d = origin - p1;
  double b1_ = d.dot(s1) * invDivisor;
  if (b1_ < 0.f || b1_ > 1.f)
    return false;
  // Compute second barycentric coordinate
  Vector3d s2 = cross(d, e1);
  double b2_ = direction.dot(s2) * invDivisor;
  if (b2_ < 0. || b1_ + b2_ > 1.)
    return false;
  // Compute t to intersection point
  double t = e2.dot(s2) * invDivisor;
  if (t < 0.)
    return false;
  if (tHit != nullptr)
    *tHit = t;
  if (b1)
    *b1 = b1_;
  if (b2)
    *b2 = b2_;
  return true;
}

bool GeometricQueries::segmentRayIntersection(const Point2d &a,
                                              const Point2d &b,
                                              const Point2d &origin,
                                              const Vector2d &direction,
                                              double *t) {
  Vector2d dd = b - a;
  auto ro = Point3d(origin.x(), origin.y(), 0.);
  auto rd = Vector3d(direction.x(), direction.y(), 0.);
  auto r2o = Point3d(a.x(), a.y(), 0.f);
  auto r2d = Vector3d(dd.x(), dd.y(), 0.f);
  double d = cross(rd, r2d).length2();
  if (d == 0.f)
    // TODO handle this case
    return false;
  double t1 = Matrix3d({r2o - ro, r2d, cross(rd, r2d)}).determinant() / d;
  double t2 = Matrix3d({r2o - ro, rd, cross(rd, r2d)}).determinant() / d;
  // TODO add user threshold
  double dist = ((ro + t1 * rd) - (r2o + t2 * r2d)).length();
  if (t1 >= 0.f && t2 >= 0.f && t2 <= dd.length() && dist < 0.1f) {
    if (t)
      *t = t1;
    return true;
  }
  return false;
}

} // namespace furoo
