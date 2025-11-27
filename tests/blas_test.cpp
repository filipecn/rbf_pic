#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(BLAS, mvm) {
  LinearMatrix m(4, 4, linspace(1, 16, 16));
  LinearVector v(linspace(1, 4, 4));
  LinearVector r(linspace(1, 4, 4));
  MatrixConstAccessor acc(&m);
  Blas::mvm(acc, &v, &r);
  for (size_t i = 0; i < 4; i++) {
    double s = 0;
    for (size_t j = 0; j < 4; j++)
      s += m(i, j) * v[j];
    ASSERT_NEAR(s, r[i], 1e-8);
  }
  MatrixConstAccessor acct(&m, true);
  Blas::mvm(acct, &v, &r);
  for (size_t i = 0; i < 4; i++) {
    double s = 0;
    for (size_t j = 0; j < 4; j++)
      s += m(j, i) * v[j];
    ASSERT_NEAR(s, r[i], 1e-8);
  }
}

TEST(Blas, axpy) {
  LinearVector x(linspace(1, 16, 16));
  LinearVector y(linspace(1, 16, 16));
  LinearVector r(linspace(1, 16, 16));
  Blas::axpy(-2, &x, &y, &r);
  for (size_t i = 0; i < 16; i++)
    ASSERT_NEAR(r[i], -2.0 * x[i] + y[i], 1e-8);
}

TEST(Blas, residual) {
  LinearMatrix m(4, 4, linspace(1, 16, 16));
  LinearVector x(linspace(1, 4, 4));
  LinearVector b(linspace(1, 4, 4));
  LinearVector r(linspace(1, 4, 4));
  MatrixConstAccessor acc(&m);
  Blas::residual(acc, &x, &b, &r);
  for (size_t i = 0; i < 4; i++) {
    double s = 0;
    for (size_t j = 0; j < 4; j++)
      s += m(i, j) * x[j];
    ASSERT_NEAR(r[i], b[i] - s, 1e-8);
  }
  MatrixConstAccessor acct(&m, true);
  Blas::residual(acct, &x, &b, &r);
  for (size_t i = 0; i < 4; i++) {
    double s = 0;
    for (size_t j = 0; j < 4; j++)
      s += m(j, i) * x[j];
    ASSERT_NEAR(r[i], b[i] - s, 1e-8);
  }
}

TEST(LinearMatrix, norm) {
  LinearVector x(linspace(1, 16, 16));
  double s = Blas::norm(&x);
  double sum = 0;
  for (size_t i = 0; i < x.size(); i++)
    sum += x[i] * x[i];
  ASSERT_NEAR(s, sqrt(sum), 1e-8);
}

TEST(LinearMatrix, infnorm) {
  LinearVector x(linspace(1, 16, 16));
  x[10] = -10.0;
  double s = Blas::infnorm(&x);
  double sum = 0;
  for (size_t i = 0; i < x.size(); i++)
    sum = std::max(sum, fabs(x[i]));
  ASSERT_NEAR(sum, s, 1e-8);
}

TEST(LinearMatrix, dSolve) {
  LinearMatrix m(4, 4, linspace(1, 16, 16));
  LinearVector b(linspace(2, 5, 4));
  LinearVector x(linspace(1, 4, 4));
  MatrixConstAccessor acc(&m);
  Blas::dSolve(acc, &b, &x);
  for (size_t i = 0; i < x.size(); i++)
    ASSERT_NEAR(x[i], ((m(i, i) != 0) ? b[i] / m(i, i) : b[i]), 1e-8);
}

TEST(LinearMatrix, dot) {
  LinearVector b(linspace(2, 5, 4));
  LinearVector x(linspace(1, 4, 4));
  double d = Blas::dot(&b, &x);
  double sum = 0;
  for (size_t i = 0; i < x.size(); i++)
    sum += x[i] * b[i];
  ASSERT_NEAR(d, sum, 1e-8);
}
