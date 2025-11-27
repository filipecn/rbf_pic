#ifndef FUROO_BLAS_LINEAR_VECTOR_H
#define FUROO_BLAS_LINEAR_VECTOR_H

#include <blas/matrix_interface.h>
#include <vector>

namespace furoo {
class LinearVector : public VectorInterface {
public:
  LinearVector() {}
  LinearVector(std::vector<double> v) : _v(v) {}
  LinearVector(size_t n, double v = 0.) {
    if (n > 0u)
      _v = std::vector<double>(n, v);
  }
  LinearVector &operator=(const std::vector<double> &v) {
    _v = v;
    return *this;
  }
  void emplace_back(double v) { _v.emplace_back(v); }
  void resize(size_t i, double v = 0.0) { _v.resize(i, v); }
  double operator[](size_t i) const override { return _v[i]; }
  double &operator[](size_t i) override { return _v[i]; }
  LinearVector operator+(LinearVector);
  LinearVector operator-(LinearVector);
  size_t size() const override { return _v.size(); }
  void setZero() { memset(&_v[0], 0, _v.size() * sizeof(double)); }

  friend std::ostream &operator<<(std::ostream &os, const LinearVector &v) {
    os << "LinearVector[" << v.size() << "]: ";
    for (size_t i = 0; i < v.size(); i++) {
      os << v[i] << " ";
    }
    os << std::endl;
    return os;
  }

private:
  std::vector<double> _v;
};
} // furoo namespace

#endif
