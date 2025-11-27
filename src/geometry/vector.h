#ifndef FUROO_GEOMETRY_VECTOR_H
#define FUROO_GEOMETRY_VECTOR_H

#include "common/debug.h"
#include "geometry/numeric.h"

#include <cstring>
#include <initializer_list>
#include <vector>

namespace furoo {
template <typename T = double, size_t D = 3> class Vector {
public:
  /// Constructors
  Vector();
  Vector(std::initializer_list<T> values);
  Vector(size_t n, const T *t);
  Vector(const T &t);
  Vector(const T &x, const T &y);
  Vector(const T &x, const T &y, const T &z);
  Vector(const Vector<T, D> &other);
  /// Operators
  T operator[](int i) const;
  T &operator[](int i);
  bool operator==(const Vector<T, D> &_v) const;
  bool operator!=(const Vector<T, D> &_v) const;
  bool operator<=(const Vector<T, D> &_v) const;
  bool operator<(const Vector<T, D> &_v) const;
  bool operator>=(const Vector<T, D> &_v) const;
  bool operator>(const Vector<T, D> &_v) const;
  Vector<T, D> operator-(const Vector<T, D> &_v) const;
  Vector<T, D> operator+(const Vector<T, D> &_v) const;
  Vector<T, D> operator*(const Vector<T, D> &_v) const;
  Vector<T, D> operator*(const T &f) const;
  Vector<T, D> operator/(const Vector<T, D> &_v) const;
  Vector<T, D> operator/=(T f);
  Vector<T, D> operator+=(const Vector<T, D> &_v);
  Vector<T, D> operator-=(const Vector<T, D> &_v);
  Vector<T, D> operator/(T f) const;
  Vector<T, 2> xy(size_t x = 0, size_t y = 1) const;
  Vector<double, 2> doubleXY(size_t x = 0, size_t y = 1) const;
  Vector<double, 3> doubleXYZ(size_t x = 0, size_t y = 1, size_t z = 2) const;
  T x() const { return v[0]; }
  T y() const { return v[1]; }
  T z() const { return v[2]; }
  size_t dimension() const { return D; }

  // Methods
  T max() const;
  double length2() const;
  double length() const;
  Vector<double, D> normalized() const;
  Vector<T, 2> right() const;
  Vector<T, 2> left() const;
  double dot(const Vector<T, D> &) const;

  ///
  /// Project self vector into given target vecot
  /// \param target vector in which projection occur
  /// \return Vector<T, D> a new projected vector
  Vector<T, D> projectInto(Vector<T, D> target);

  friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
    os << "Vector[<" << D << ">]";
    for (size_t i = 0; i < D; i++) {
      os << v[i] << " ";
    }
    return os;
  }

private:
  T v[D];
};

template <typename T, size_t D>
Vector<T, D> operator*(T f, const Vector<T, D> &v);
template <typename T> Vector<int, 3> round(const Vector<T, 3> &v);
template <typename T> Vector<int, 2> round(const Vector<T, 2> &v);
template <typename T> Vector<int, 3> ceil(const Vector<T, 3> &v);
template <typename T> Vector<int, 2> ceil(const Vector<T, 2> &v);
template <typename T> Vector<int, 3> floor(const Vector<T, 3> &v);
template <typename T> Vector<int, 2> floor(const Vector<T, 2> &v);
template <typename T> Vector<T, 3> min(Vector<T, 3> a, Vector<T, 3> b);
template <typename T> Vector<T, 3> max(Vector<T, 3> a, Vector<T, 3> b);
template <typename T> Vector<T, 2> min(Vector<T, 2> a, Vector<T, 2> b);
template <typename T> Vector<T, 2> max(Vector<T, 2> a, Vector<T, 2> b);
template <typename T, size_t D> void normalize(Vector<T, D> &v);
template <typename T>
inline Vector<T, 3> cross(const Vector<T, 3> &a, const Vector<T, 3> &b) {
  return Vector<T, 3>((a.y() * b.z()) - (a.z() * b.y()),
                      (a.z() * b.x()) - (a.x() * b.z()),
                      (a.x() * b.y()) - (a.y() * b.x()));
}

#include "geometry/vector.inl"

typedef Vector<size_t, 2> Vector2u;
typedef Vector<size_t, 3> Vector3u;
typedef Vector<double, 2> Vector2d;
typedef Vector<double, 3> Vector3d;
typedef Vector<int, 2> Vector2i;
typedef Vector<int, 3> Vector3i;
} // namespace furoo

#endif
