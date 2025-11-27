// DEPRECATED
#ifndef FUROO_STRUCTURES_REGULAR_GRID_H
#define FUROO_STRUCTURES_REGULAR_GRID_H

#include <geometry/transform.h>

namespace furoo {
class RegularGrid {
 public:
  RegularGrid();
  RegularGrid(size_t n);
  RegularGrid(size_t n, const Transform2D &t);
  RegularGrid(size_t n, double x);
  RegularGrid(size_t n, size_t m, const Transform2D &t);
  RegularGrid(size_t n, size_t m, const BBox2D &r);
  RegularGrid(size_t n, size_t m, double x, double y);
  virtual ~RegularGrid();
  Point2d fieldPosition(size_t, size_t);
  Vector<size_t, 2> gridSize() const;
  size_t cellCount() const;
  Point2d transformToWorld(Point2d p) const;
  Point2d transformToGrid(Point2d p) const;
  Transform2D getGridTransform() const { return _toGrid; }
  Vector2d spacing() const;

 private:
  double _hx, _hy;
  size_t _n, _m;
  Transform2D _toGrid;
  Transform2D _toWorld;
};

class RegularGrid3 {
 public:
  RegularGrid3();
  RegularGrid3(unsigned n);
  RegularGrid3(unsigned n, const Transform3D &t);
  RegularGrid3(unsigned n, double x);
  RegularGrid3(unsigned n, unsigned m, unsigned p, const Transform3D &t);
  RegularGrid3(unsigned n, unsigned m, unsigned p, const BBox3D &r);
  RegularGrid3(unsigned n, unsigned m, unsigned p, double x, double y,
               double z);
  virtual ~RegularGrid3();
  Point3d fieldPosition(unsigned, unsigned, unsigned);
  Vector<unsigned, 3> gridSize() const;
  unsigned cellCount() const;
  Point3d transformToWorld(Point3d p) const;
  Point3d transformToGrid(Point3d p) const;
  Transform3D getGridTransform() const { return _toGrid; }
  Vector3d spacing() const;

 private:
  double _hx, _hy, _hz;
  unsigned _n, _m, _p;
  Transform3D _toGrid;
  Transform3D _toWorld;
};

}  // namespace furoo

#endif  // FUROO_STRUCTURES_DOMAIN_H
