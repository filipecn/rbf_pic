#ifndef FUROO_BLAS_CUSTOM_MATRIX_H
#define FUROO_BLAS_CUSTOM_MATRIX_H

#include <assert.h>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>

namespace furoo {

class VectorInterface {
public:
  VectorInterface() {}
  virtual ~VectorInterface() {}
  virtual double operator[](size_t i) const = 0;
  virtual double &operator[](size_t i) = 0;
  virtual size_t size() const = 0;
};

class MatrixInterface {
public:
  MatrixInterface() {}
  virtual ~MatrixInterface() {}
  virtual double operator()(size_t i, size_t j) const = 0;
  virtual double &operator()(size_t i, size_t j) = 0;
  virtual void
  iterateRow(size_t i, std::function<void(const double &, size_t)> f) const = 0;
  virtual void iterateRow(size_t i,
                          std::function<void(double &, size_t)> f) = 0;
  virtual void
  iterateColumn(size_t j,
                std::function<void(const double &, size_t)> f) const = 0;
  virtual void iterateColumn(size_t j,
                             std::function<void(double &, size_t)> f) = 0;
  virtual size_t rowCount() const = 0;
  virtual size_t columnCount() const = 0;
};

class MatrixConstAccessor {
public:
  MatrixConstAccessor(const MatrixInterface *m, bool t = false)
      : _m(m), transposeAccess(t) {}
  virtual ~MatrixConstAccessor() {}
  double operator()(size_t i, size_t j) const { return (*_m)(i, j); }
  void iterateRow(size_t i,
                  std::function<void(const double &, size_t)> f) const {
    if (transposeAccess) {
      _m->iterateColumn(i, f);
    } else {
      _m->iterateRow(i, f);
    }
  }

  size_t rowCount() const { return _m->rowCount(); }
  size_t columnCount() const { return _m->columnCount(); }

private:
  const MatrixInterface *_m;
  bool transposeAccess;
};
} // furoo namespace
#endif
