#include <blas/linear_vector.h>
#include <cmath>
#include <algorithm>

namespace furoo {
LinearVector LinearVector::operator+(LinearVector v) {
  LinearVector w(_v.size());

  for (size_t i = 0; i < _v.size(); i++) {
    w[i] = _v[i] + v[i];
  }
  return w;
}
LinearVector LinearVector::operator-(LinearVector v) {
  size_t size = std::min(_v.size(), v.size());
  LinearVector w(size);

  for (size_t i = 0; i < size; i++) {
    w[i] = _v[i] - v[i];
  }
  return w;
}
} // namespace furoo
