#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(LinearMatrix, Constructors) {
  LinearMatrix m(10);
  for (size_t i = 0; i < m.rowCount(); i++)
    for (size_t j = 0; j < m.columnCount(); j++)
      m(i, j) = i * 10 + j;
  LinearMatrix m2(m);
  for (size_t i = 0; i < m.rowCount(); i++)
    for (size_t j = 0; j < m.columnCount(); j++)
      ASSERT_NEAR(m2(i, j), i * 10. + j, 1e-8);
}

TEST(LinearMatrix, Iterators) {
  LinearMatrix lm;
  EXPECT_EQ(lm.rowCount(), 0u);
  EXPECT_EQ(lm.columnCount(), 0u);
  lm.resize(10, 100);
  EXPECT_EQ(lm.rowCount(), 10u);
  EXPECT_EQ(lm.columnCount(), 100u);
  for (size_t i = 0; i < 10; i++)
    for (size_t j = 0; j < 100; j++) {
      ASSERT_NEAR(lm(i, j), 0.0, 1e-8);
      lm(i, j) = 10.0 * i + 1.0 * j;
      ASSERT_NEAR(lm(i, j), 10.0 * i + 1.0 * j, 1e-8);
    }
  LinearMatrix lm2 = lm;
  for (size_t i = 0; i < 10; i++)
    lm2.iterateRow(i, [i](const double &v, size_t j) {
      ASSERT_NEAR(i * 10.0 + 1.0 * j, v, 1e-8);
    });
  for (size_t j = 0; j < 100; j++)
    lm2.iterateColumn(j, [j](const double &v, size_t i) {
      ASSERT_NEAR(10.0 * i + 1.0 * j, v, 1e-8);
    });
  for (size_t i = 0; i < 10; i++)
    lm2.iterateRow(i, [i](double &v, size_t j) {
      ASSERT_NEAR(i * 10.0 + 1.0 * j, v, 1e-8);
    });
  for (size_t j = 0; j < 100; j++)
    lm2.iterateColumn(j, [j](double &v, size_t i) {
      ASSERT_NEAR(10.0 * i + 1.0 * j, v, 1e-8);
    });
}

TEST(LinearMatrix, ConstAccessor) {
  LinearMatrix m(4, 4, linspace(1.0, 16.0, 16));
  MatrixConstAccessor acc(&m);
  EXPECT_EQ(acc.rowCount(), m.rowCount());
  EXPECT_EQ(acc.columnCount(), m.columnCount());
  for (size_t i = 0; i < 4; i++)
    acc.iterateRow(i, [](const double &v, size_t j) {
      UNUSED_VARIABLE(j);
      static int k = 1;
      ASSERT_NEAR(1.0 * k, v, 1e-8);
      k++;
    });
  MatrixConstAccessor acct(&m, true);
  for (size_t i = 0; i < 4; i++) {
    int k = i + 1;
    acct.iterateRow(i, [&k](const double &v, size_t j) {
      UNUSED_VARIABLE(j);
      ASSERT_NEAR(1.0 * k, v, 1e-8);
      k += 4;
    });
  }
}

TEST(LinearMatrix, Operators) {
  {
    LinearMatrix a(4, 4, linspace(1., 16., 16));
    LinearMatrix b(4, 4, linspace(1., 16., 16));
    LinearMatrix m = a * b;
    for (size_t r = 0; r < 4; r++)
      for (size_t c = 0; c < 4; c++) {
        double e = 0.;
        for (size_t i = 0; i < 4; i++)
          e += a(r, i) * b(i, c);
        ASSERT_NEAR(e, m(r, c), 1e-8);
      }
  }
  {
    LinearMatrix m(4, 5, linspace(1., 16., 20));
    LinearVector v(linspace(1., 5., 5));
    LinearVector r = m * v;
    for (size_t i = 0; i < 4; i++) {
      double e = 0;
      for (size_t j = 0; j < 5; j++)
        e += m(i, j) * v[j];
      ASSERT_NEAR(e, r[i], 1e-8);
    }
  }
}
