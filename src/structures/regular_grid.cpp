#include <common/debug.h>
#include <structures/regular_grid.h>

namespace furoo {
RegularGrid::RegularGrid() : _n(0), _m(0) {}

RegularGrid::~RegularGrid() {}

RegularGrid::RegularGrid(size_t n) : _n(n), _m(n) { _hx = _hy = 1; }

RegularGrid::RegularGrid(size_t n, const Transform2D &t)
    : _n(n), _m(n), _toWorld(t) {
  _hx = (t(Point2d(1, 0)) - t(Point2d(0, 0)))[0];
  _hy = (t(Point2d(0, 1)) - t(Point2d(0, 0)))[1];
  _toGrid = _toWorld.inverse();
}

RegularGrid::RegularGrid(size_t n, double x) : _n(n), _m(n) {
  _hx = _hy = x / n;
  _toWorld =
      Transform2D::translate(Vector2d(_hx / 2, _hy / 2)).applyScale(_hx, _hy);
  _toGrid = _toWorld.inverse();
}

RegularGrid::RegularGrid(size_t n, size_t m, const Transform2D &t)
    : _n(n), _m(m), _toWorld(t) {
  _toGrid = _toWorld.inverse();
}

RegularGrid::RegularGrid(size_t n, size_t m, const BBox2D &r) : _n(n), _m(m) {
  _toGrid = Transform2D::scale(n, m) * Transform2D(r);
  _toWorld = _toGrid.inverse();
}

RegularGrid::RegularGrid(size_t n, size_t m, double x, double y)
    : _n(n), _m(m) {
  _hx = x / n;
  _hy = y / m;
  _toWorld =
      Transform2D::translate(Vector2d(_hx / 2, _hy / 2)).applyScale(_hx, _hy);
  _toGrid = _toWorld.inverse();
}

Vector<size_t, 2> RegularGrid::gridSize() const {
  return Vector<size_t, 2>({_n, _m});
}

size_t RegularGrid::cellCount() const { return _n * _m; }

Point2d RegularGrid::transformToWorld(Point2d p) const { return _toWorld(p); }

Point2d RegularGrid::transformToGrid(Point2d p) const { return _toGrid(p); }

Point2d RegularGrid::fieldPosition(size_t i, size_t j) {
  return _toWorld(Point2d(i, j));
}

Vector2d RegularGrid::spacing() const { return Vector2d(_hx, _hy); }

RegularGrid3::RegularGrid3() : _n(0), _m(0), _p(0) {}

RegularGrid3::~RegularGrid3() {}

RegularGrid3::RegularGrid3(unsigned n) : _n(n), _m(n), _p(n) {
  _hx = _hy = _hz = 1;
}

RegularGrid3::RegularGrid3(unsigned n, const Transform3D &t)
    : _n(n), _m(n), _p(n), _toWorld(t) {
  _hx = (t(Point3d(1, 0, 0)) - t(Point3d(0.)))[0];
  _hy = (t(Point3d(0, 1, 0)) - t(Point3d(0.)))[1];
  _hz = (t(Point3d(0, 0, 1)) - t(Point3d(0.)))[2];
  _toGrid = _toWorld.inverse();
}

RegularGrid3::RegularGrid3(unsigned n, double x) : _n(n), _m(n), _p(n) {
  _hx = _hy = _hz = x / n;
  _toWorld = Transform3D::translate(Vector3d(_hx / 2, _hy / 2, _hz / 2))
                 .applyScale(_hx, _hy, _hz);
  _toGrid = _toWorld.inverse();
}

RegularGrid3::RegularGrid3(unsigned n, unsigned m, unsigned p,
                           const Transform3D &t)
    : _n(n), _m(m), _p(p), _toWorld(t) {
  _toGrid = _toWorld.inverse();
}

RegularGrid3::RegularGrid3(unsigned n, unsigned m, unsigned p, const BBox3D &r)
    : _n(n), _m(m), _p(p) {
  _toGrid = Transform3D::scale(n, m, p) * Transform3D(r);
  _toWorld = _toGrid.inverse();
}

RegularGrid3::RegularGrid3(unsigned n, unsigned m, unsigned p, double x,
                           double y, double z)
    : _n(n), _m(m), _p(p) {
  _hx = x / n;
  _hy = y / m;
  _hz = z / p;
  _toWorld = Transform3D::translate(Vector3d(_hx / 2, _hy / 2, _hz / 2))
                 .applyScale(_hx, _hy, _hz);
  _toGrid = _toWorld.inverse();
}

Vector<unsigned, 3> RegularGrid3::gridSize() const {
  return Vector<unsigned, 3>({_n, _m, _p});
}

unsigned RegularGrid3::cellCount() const { return _n * _m * _p; }

Point3d RegularGrid3::transformToWorld(Point3d p) const { return _toWorld(p); }

Point3d RegularGrid3::transformToGrid(Point3d p) const { return _toGrid(p); }

Point3d RegularGrid3::fieldPosition(unsigned i, unsigned j, unsigned k) {
  return _toWorld(Point3d(i, j, k));
}

Vector3d RegularGrid3::spacing() const { return Vector3d(_hx, _hy, _hz); }

} // namespace furoo
