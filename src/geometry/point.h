#ifndef FUROO_GEOMETRY_POINT_H
#define FUROO_GEOMETRY_POINT_H

#include <initializer_list>

#include "common/debug.h"
#include "geometry/vector.h"

namespace furoo {
template <class T, int D> class Point {
public:
  Point();
  explicit Point(const std::vector<T> &v);
  explicit Point(const T *v);
  Point(T V);
  Point(std::initializer_list<T>);
  Point(T x_, T y_);
  Point(T x_, T y_, T z_);

  operator float() const { return 0.; }

  T operator[](int i) const {
    ASSERT_FATAL(i >= 0 && i <= static_cast<int>(D));
    return _v[i];
  }

  T &operator[](int i) {
    THROW(i >= 0 && i < static_cast<int>(D),
          "Point<>::operator[] invalid index");
    return _v[i];
  }

  Point<T, D> &operator+=(const Vector<T, D> &V) {
    for (int i = 0; i < D; i++) {
      _v[i] += V[i];
    }
    return *this;
  }

  Point<T, D> &operator-=(const Vector<T, D> &V) {
    for (int i = 0; i < D; i++) {
      _v[i] -= V[i];
    }
    return *this;
  }

  bool operator==(const Point<T, D> &) const;

  bool operator!=(const Point<T, D> &rhs) { return !(*this == rhs); }

  bool operator>=(const Point<T, D> &p) const;

  bool operator>(const Point<T, D> &p) const {
    for (int i = 0; i < D; i++) {
      if (_v[i] <= p[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator<=(const Point<T, D> &p) const {
    for (int i = 0; i < D; i++) {
      if (_v[i] > p[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator<(const Point<T, D> &p) const {
    for (int i = 0; i < D; i++) {
      if (_v[i] >= p[i]) {
        return false;
      }
    }
    return true;
  }

  Point<T, D> operator+(const Point<T, D> &V) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] + V[i];
    }
    return P;
  }

  Point<T, D> operator+(const Vector<T, D> &V) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] + V[i];
    }
    return P;
  }

  Point<T, D> &operator+=(const Point<T, D> &p) {
    for (int i = 0; i < D; i++) {
      _v[i] += p[i];
    }
    return *this;
  }

  Point<T, D> &operator+=(const double d) {
    for (int i = 0; i < D; i++) {
      _v[i] += d;
    }
    return *this;
  }

  Point<T, D> operator+(const double &f) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] + f;
    }
    return P;
  }

  Vector<T, D> operator-(const Point<T, D> &p) const {
    Vector<T, D> V;
    for (int i = 0; i < D; i++) {
      V[i] = _v[i] - p[i];
    }
    return V;
  }

  Point<T, D> operator-(const Vector<T, D> &V) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] - V[i];
    }
    return P;
  }

  Point<T, D> operator-(const double &f) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] - f;
    }
    return P;
  }

  Point<T, D> operator*(const T d) const {
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = static_cast<T>(_v[i] * d);
    }
    return P;
  }

  Point<T, D> operator/(const T d) const {
    ASSERT_FATAL(IS_NOT_ZERO(d));
    Point<T, D> P;
    for (int i = 0; i < D; i++) {
      P[i] = _v[i] / d;
    }
    return P;
  }

  Point<T, D> &operator/=(const double d) {
    ASSERT_FATAL(IS_NOT_ZERO(d));
    for (int i = 0; i < D; i++) {
      _v[i] /= d;
    }
    return *this;
  }

  void setX(T value) { _v[0] = value; }
  void setY(T value) { _v[1] = value; }
  void setZ(T value) { _v[2] = value; }

  T x() const { return _v[0]; }

  T y() const { return _v[1]; }

  T z() const { return _v[2]; }

  Point<T, 2> doubleXY(size_t x = 0, size_t y = 1) const {
    return Point<T, 2>(static_cast<double>(_v[x]), static_cast<double>(_v[y]));
  }

  double distance(const Point<T, D> &b) { return (*this - b).length(); }

  double distance2(const Point<T, D> &b) { return (*this - b).length2(); }

  size_t dimension() const { return D; }

  friend std::ostream &operator<<(std::ostream &os, const Point &p) {
    os << "[Point<" << D << ">]";
    for (size_t i = 0; i < D; i++) {
      os << p[i] << " ";
    }
    return os;
  }

  const T *asConstArray() const { return _v; }

private:
  T _v[D];
};

#include "geometry/point.inl"

typedef Point<size_t, 2> Point2u;
typedef Point<size_t, 3> Point3u;
typedef Point<int, 2> Point2i;
typedef Point<double, 2> Point2d;
typedef Point<int, 3> Point3i;
typedef Point<double, 3> Point3d;
} // namespace furoo
#endif
