#ifndef FUROO_GEOMETRY_QUERIES_H
#define FUROO_GEOMETRY_QUERIES_H

#include <geometry/bbox.h>

namespace furoo {

template <typename T, int D>
bool bbox_bbox_intersection(const BBox<T, D> &a, const BBox<T, D> &b) {
  for (int i = 0; i < D; i++)
    if (!((a.lower()[i] <= b.lower()[i] && a.upper()[i] >= b.lower()[i]) ||
          (b.lower()[i] <= a.lower()[i] && b.upper()[i] >= a.lower()[i])))
      return false;
  return true;
}

bool bbox_bbox_intersection(const BBox2D &a, const BBox2D &b);

class GeometricQueries {
public:
  static bool triangleRayIntersection(const Point3d &p1, const Point3d &p2,
                                      const Point3d &p3, const Point3d &origin,
                                      const Vector3d &direction,
                                      double *tHit = nullptr,
                                      double *b1 = nullptr,
                                      double *b2 = nullptr);

  static bool segmentRayIntersection(const Point2d &a, const Point2d &b,
                                     const Point2d &origin,
                                     const Vector2d &direction,
                                     double *t = nullptr);
};

} // namespace furoo

#endif // FUROO_GEOMETRY_QUERIES_H
