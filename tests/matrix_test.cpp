#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Matrix, Constructors) {
  Matrix2d v({0.1, 0.2, 1.1, 1.2});
  ASSERT_NEAR(0.1, v(0, 0), 1e-8);
  ASSERT_NEAR(0.2, v(0, 1), 1e-8);
  ASSERT_NEAR(1.1, v(1, 0), 1e-8);
  ASSERT_NEAR(1.2, v(1, 1), 1e-8);
  Matrix2d b({0.1, 1.1, 0.2, 1.2}, true);
  ASSERT_NEAR(0.1, b(0, 0), 1e-8);
  ASSERT_NEAR(1.1, b(1, 0), 1e-8);
  ASSERT_NEAR(0.2, b(0, 1), 1e-8);
  ASSERT_NEAR(1.2, b(1, 1), 1e-8);
  Matrix<double, 100> a;
  for (size_t i = 0; i < 100; i++)
    for (size_t j = 0; j < 100; j++)
      ASSERT_NEAR(a(i, j), 0.0, 1e-8);
}

TEST(Matrix, Access) {
  Matrix3d b;
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++)
      b(i, j) = i * 10.0 + j;
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++)
      ASSERT_NEAR(b(i, j), i * 10.0 + j, 1e-8);
}

TEST(Matrix, Multiplication) {
  Vector2d v(1., 2.);
  Matrix2d m({2.0, 1.0, 4.0, 2.0});
  Vector2d M = m * v;
  ASSERT_NEAR(M[0], 4.0, 1e-8);
  ASSERT_NEAR(M[1], 8.0, 1e-8);
  Matrix2d m2({-2.0, -1.0, -4.0, -2.0});
  Matrix2d r = m2 * m;
  ASSERT_NEAR(r(0, 0), -8.0, 1e-8);
  ASSERT_NEAR(r(0, 1), -4.0, 1e-8);
  ASSERT_NEAR(r(1, 0), -16.0, 1e-8);
  ASSERT_NEAR(r(1, 1), -8.0, 1e-8);
}

TEST(Matrix, Idendity) {
  Matrix<double, 100> a;
  a.setIdentity();
  for (size_t i = 0; i < 100; i++)
    for (size_t j = 0; j < 100; j++)
      if (i == j)
        ASSERT_NEAR(a(i, j), 1.0, 1e-8);
      else
        ASSERT_NEAR(a(i, j), 0.0, 1e-8);
  EXPECT_EQ(a.isIdentity(), true);
}

TEST(Matrix, Transpose) {
  Matrix3d a({1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3});
  Matrix3d t = a.transpose();
  ASSERT_NEAR(t(0, 0), 1.1, 1e-8);
  ASSERT_NEAR(t(0, 1), 2.1, 1e-8);
  ASSERT_NEAR(t(0, 2), 3.1, 1e-8);
  ASSERT_NEAR(t(1, 0), 1.2, 1e-8);
  ASSERT_NEAR(t(1, 1), 2.2, 1e-8);
  ASSERT_NEAR(t(1, 2), 3.2, 1e-8);
  ASSERT_NEAR(t(2, 0), 1.3, 1e-8);
  ASSERT_NEAR(t(2, 1), 2.3, 1e-8);
  ASSERT_NEAR(t(2, 2), 3.3, 1e-8);
}
