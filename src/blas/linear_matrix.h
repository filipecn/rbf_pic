#ifndef FUROO_BLAS_LINEAR_MATRIX_H
#define FUROO_BLAS_LINEAR_MATRIX_H

#include <blas/linear_vector.h>
#include <blas/matrix_interface.h>
#include <common/debug.h>
#include <vector>

namespace furoo {

class LinearMatrix : public MatrixInterface {
public:
  LinearMatrix() : _N(0), _M(0), _data(nullptr) {}
  LinearMatrix(size_t n) : _N(n), _M(n), _data(nullptr) { resize(n, n); }
  LinearMatrix(size_t n, size_t m) : _N(n), _M(m), _data(nullptr) { resize(n, m); }
  LinearMatrix(size_t n, size_t m, std::vector<double> v) : LinearMatrix(n, m) {
    for (size_t i = 0; i < _N; i++)
      for (size_t j = 0; j < _M; j++)
        _data[i][j] = v[i * _N + j];
  }
  LinearMatrix(const LinearMatrix &A) : LinearMatrix() {
    resize(A._N, A._M);
    // TODO make this copy faster
    for (size_t r = 0; r < A._N; r++)
      for (size_t c = 0; c < A._M; c++)
        _data[r][c] = A._data[r][c];
  }
  virtual ~LinearMatrix() {
    if (_data)
      clearData();
  }
  LinearMatrix &operator=(const LinearMatrix& A) {
    if (this == &A)
      return *this;

    resize(A._N, A._M);
    for (size_t r = 0; r < A._N; r++)
      std::copy(&A._data[r][0], &A._data[r][0] + A._M, &_data[r][0]);
    return *this;
  }
  LinearMatrix operator*(const LinearMatrix &A) const {
    ASSERT_FATAL(_M == A._N && _N == A._M);
    LinearMatrix R(_N, A._M);
    for (size_t r = 0; r < R._N; r++)
      for (size_t c = 0; c < R._M; c++) {
        R(r, c) = 0;
        for (size_t i = 0; i < _N; i++)
          R(r, c) += _data[r][i] * A._data[i][c];
      }
    return R;
  }
  LinearVector operator*(const LinearVector &v) const {
    ASSERT_FATAL(_M == v.size());
    LinearVector r(_N);
    for (size_t i = 0; i < _N; ++i) {
      r[i] = 0.;
      for (size_t j = 0; j < _M; j++)
        r[i] += v[j] * _data[i][j];
    }
    return r;
  }
  void resize(size_t rows, size_t columns) {
    _N = rows, _M = columns;
    if (_data) {
      clearData();
    }
    _data = new double *[_N];
    for (size_t i = 0; i < _N; i++) {
      _data[i] = new double[_M];
      memset(_data[i], 0, sizeof(double) * _M);
    }
  }
  double operator()(size_t i, size_t j) const override { return _data[i][j]; }
  double &operator()(size_t i, size_t j) override { return _data[i][j]; }
  void
  iterateRow(size_t i,
             std::function<void(const double &, size_t)> f) const override {
    for (size_t j = 0; j < _M; j++)
      f(_data[i][j], j);
  }
  void iterateRow(size_t i, std::function<void(double &, size_t)> f) override {
    for (size_t j = 0; j < _M; j++)
      f(_data[i][j], j);
  }
  void
  iterateColumn(size_t j,
                std::function<void(const double &, size_t)> f) const override {
    for (size_t i = 0; i < _N; i++)
      f(_data[i][j], i);
  }
  void iterateColumn(size_t j,
                     std::function<void(double &, size_t)> f) override {
    for (size_t i = 0; i < _N; i++)
      f(_data[i][j], i);
  }
  size_t rowCount() const override { return _N; }
  size_t columnCount() const override { return _M; }

  friend std::ostream &operator<<(std::ostream &os, const LinearMatrix &m) {
    for (size_t i = 0; i < m._N; i++) {
      for (size_t j = 0; j < m._M; j++)
        os << m(i, j) << " ";
      os << std::endl;
    }
    return os;
  }

private:
  size_t _N, _M;
  double **_data;

  // methods:
  void clearData()
  {
    for (size_t i = 0; i < _N; i++)
      delete[] _data[i];
    delete[] _data;
  }
};

} // namespace furoo

#endif
