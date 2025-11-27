#include <algorithm>
#include <blas/blas.h>
#include <common/debug.h>
#include <geometry/numeric.h>

namespace furoo {

void Blas::set(const VectorInterface *b, VectorInterface *a) {
  for (unsigned i = 0; i < b->size(); i++)
    (*a)[i] = (*b)[i];
}

void Blas::mvm(const MatrixConstAccessor &m, const VectorInterface *v,
               VectorInterface *r) {
  for (size_t i = 0; i < m.rowCount(); i++) {
    double s = 0.0;
    m.iterateRow(i, [&s, &v](const double &a, size_t j) { s += a * (*v)[j]; });
    (*r)[i] = s;
  }
}
void Blas::axpy(double a, const VectorInterface *x, VectorInterface *y,
                VectorInterface *r) {
  size_t size = x->size();

  for (size_t i = 0; i < size; i++) {
    (*r)[i] = a * (*x)[i] + (*y)[i];
  }
}

void Blas::residual(const MatrixConstAccessor &A, const VectorInterface *x,
                    const VectorInterface *b, VectorInterface *r) {
  mvm(A, x, r);
  for (size_t i = 0; i < b->size(); i++) {
    (*r)[i] = (*b)[i] - (*r)[i];
  }
}

double Blas::norm(const VectorInterface *v) {
  size_t size = v->size();
  double sum = 0.0;
  for (size_t i = 0; i < size; i++) {
    sum += (*v)[i] * (*v)[i];
  }
  return sqrt(sum);
}

double Blas::infnorm(const VectorInterface *v) {
  size_t size = v->size();
  double m = 0.0;
  for (size_t i = 0; i < size; i++) {
    m = std::max(m, fabs((*v)[i]));
  }
  return m;
}

void Blas::dSolve(const MatrixConstAccessor &A, const VectorInterface *b,
                  VectorInterface *x) {
  for (size_t i = 0; i < A.rowCount(); i++) {
    if (A(i, i) != 0.0) {
      (*x)[i] = (*b)[i] / A(i, i);
    } else {
      (*x)[i] = (*b)[i];
    }
  }
}

double Blas::dot(const VectorInterface *a, const VectorInterface *b) {
  double sum = 0.0;

  for (size_t i = 0; i < a->size(); i++) {
    sum += (*a)[i] * (*b)[i];
  }
  return sum;
}

} // furoo namespace
