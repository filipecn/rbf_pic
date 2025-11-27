#ifndef FUROO_GEOMETRY_BBOX_H
#define FUROO_GEOMETRY_BBOX_H

#include <algorithm>
#include <geometry/point.h>
#include <iostream>

namespace furoo {

/// Represents a bounding box
/// \tparam T coordinate type
/// \tparam D number of dimensions
template <typename T, int D> class BBox {
public:
  BBox();

  /// \param p1 lower point
  /// \param p2 upper point
  explicit BBox(const Point<T, D> &p1, const Point<T, D> &p2);

  /// \param p1 lower point
  /// \param p2 upper point
  explicit BBox(const Point<T, D> &p1, double radius);

  /// \return point with smallest coordinates
  Point<T, D> lower() const;

  /// \return point with biggest coordinates
  Point<T, D> upper() const;

  /// \return box center point
  Point<T, D> center() const;

  /// \return side sizes
  Vector<T, D> size() const;

  /// \param p query point
  /// \return true if p is inside bbox
  bool contains(Point<T, D> p) const;

  /// \param b bounding box
  /// \return true if b is fully inside bbox
  bool contains(BBox<T, D> b) const;

  /// Clamp the point to the bounding box. Each coordinate is set to limits to 
  /// the bounding box, if hte point is outside
  /// \param point point to be clamped
  /// \return Point<T,D> New point thats inside bbox limits
  Point<T,D> clamp(Point<T,D> point);

  /// \param b bounding box
  /// \return true if intersection exists
  bool intersect(BBox<T, D> b) const;

  /// Computes the expanded box
  /// \param expansionRatio percentage of expansion in each dimension
  /// \return expanded box
  BBox<T, D> expanded(Vector<double, D> expansionRatio) const;

  /// Intersects ray with the bbox
  /// \param initialPoint ray origin
  /// \param direction ray direction
  /// \param hit1 [optional] first hit parametric point
  /// \param hit2 [optional] second hit parametric point
  /// \return true if intersection exists
  bool intersect(Point<T, D> initialPoint, Vector<T, D> direction,
                 double *hit1 = nullptr, double *hit2 = nullptr) const;

  /// Splits the bounding box into 2^D bounding boxes
  /// with half of box size
  /// \param bbox
  /// \return list of bounding boxes
  static std::vector<BBox<T, D>> halfSplit(const BBox<T, D> &bbox);

  /// \param length box side size
  /// \return a square box of dimensions [0]x[length]
  static BBox<T, D> squareBox(T length = 1);

  /// Performs an union operation between boxes a and b
  /// \param a box 1
  /// \param b box 2
  /// \return the union of a and b
  static BBox<T, D> combine(const BBox<T, D> &a, const BBox<T, D> &b);

  /// Performs an union operation between box a and point b
  /// \param a box
  /// \param b point
  /// \return the union of a and b
  static BBox<T, D> combine(const BBox<T, D> &a, const Point<T, D> &b);

  friend std::ostream &operator<<(std::ostream &o, const BBox<T, D> &b) {
    o << "BBOX [";
    for (size_t d = 0; d < D; d++)
      o << b.lower()[d] << ", ";
    o << "]x[";
    for (size_t d = 0; d < D; d++)
      o << b.upper()[d] << ", ";
    o << "]";
    return o;
  }

private:
  Point<T, D> _lower;
  Point<T, D> _upper;
};

#include "bbox.inl"

typedef BBox<double, 2> BBox2d;
typedef BBox<double, 3> BBox3d;

class BBox2D {
public:
  BBox2D();
  BBox2D(const Point2d &p1, const Point2d &p2);
  BBox2D(const Point2d &c, double x, double y = -1.);
  BBox2D(const std::vector<Point2d> &);
  /// \param p target point
  /// \return coordinates of p in bbox space [0,1]x[0,1]
  Point2d pointCoordinates(const Point2d &p) const;
  bool inside(const Point2d &p) const;
  bool inside(const BBox2D &b) const;
  double size(int d) const;
  Vector2d extends() const;
  Point2d center() const;
  int maxExtent() const;
  const Point2d &operator[](int i) const;
  Point2d &operator[](int i);
  Point2d lower() const;
  Point2d upper() const;
  static BBox2D make_union(const BBox2D &b, const Point2d &p);
  static BBox2D make_union(const BBox2D &a, const BBox2D &b);
  static BBox2D make_unit_bbox();
  friend std::ostream &operator<<(std::ostream &o, const BBox2D &b) {
    o << "BBOX [" << b.lower().x() << ", " << b.lower().y() << "]x["
      << b.upper().x() << ", " << b.upper().y() << "]\n";
    return o;
  }

  bool intersect(Point2d, Vector2d, Point2d &, Vector2d &) const;
  std::vector<BBox2D> splitBy4() const {
    auto mid = center();
    std::vector<BBox2D> children;
    children.emplace_back(Point2d(_pMin.x(), mid.y()),
                          Point2d(mid.x(), _pMax.y()));
    children.emplace_back(mid, _pMax);
    children.emplace_back(_pMin, mid);
    children.emplace_back(Point2d(mid.x(), _pMin.y()),
                          Point2d(_pMax.x(), mid.y()));
    return children;
  }

private:
  Point2d _pMin, _pMax;
};

class BBox3D {
public:
  BBox3D();
  BBox3D(const Point3d &p1, const Point3d &p3);
  BBox3D(const Point3d &c, double x, double y = -1., double z = -1.);
  BBox3D(const std::vector<Point3d> &);

  bool inside(const Point3d &p) const;
  double size(int d) const;
  Vector3d extends() const;
  Point3d center() const;
  int maxExtent() const;
  const Point3d &operator[](int i) const;
  Point3d &operator[](int i);
  Point3d lower() const;
  Point3d upper() const;
  static BBox3D make_union(const BBox3D &b, const Point3d &p);
  static BBox3D make_union(const BBox3D &a, const BBox3D &b);
  static BBox3D make_unit_bbox();
  friend std::ostream &operator<<(std::ostream &o, const BBox3D &b) {
    o << "BBOX [" << b.lower().x() << ", " << b.lower().y() << "]x["
      << b.upper().x() << ", " << b.upper().y() << "]\n";
    return o;
  }

  bool intersect(Point3d, Vector3d, Point3d &, Vector3d &) const;

private:
  Point3d _pMin, _pMax;
};

} // namespace furoo

#endif
