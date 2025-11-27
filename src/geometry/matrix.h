#ifndef FUROO_GEOMETRY_MATRIX_H
#define FUROO_GEOMETRY_MATRIX_H

#include "geometry/vector.h"
#include <blas/matrix_interface.h>
#include <vector>

namespace furoo {
template <typename T = double, size_t N = 3> class Matrix {
public:
  Matrix();
  Matrix(std::initializer_list<T> vs, bool cm = false);
  Matrix(std::initializer_list<Vector<T, N>> vs, bool columns = false);
  Vector<T, N> operator*(const Vector<T, N> &v) const;
  Matrix<T, N> operator*(const Matrix<T, N> &b) const;
  Matrix<T, N> transpose() const;
  T operator()(size_t i, size_t j) const;
  T &operator()(size_t i, size_t j);
  void setIdentity();
  bool isIdentity() const;
  T determinant() const;

  friend std::ostream &operator<<(std::ostream &os, const Matrix<T, N> &ma) {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++)
        os << ma.m[i][j] << "\t";
      os << std::endl;
    }
    return os;
  }
  T m[N][N];
};

template <typename T, size_t N> Matrix<T, N> inverse(const Matrix<T, N> &m);

#include "geometry/matrix.inl"

typedef Matrix<double, 2> Matrix2d;
typedef Matrix<double, 3> Matrix3d;
typedef Matrix<double, 4> Matrix4d;
} // namespace furoo

#endif
