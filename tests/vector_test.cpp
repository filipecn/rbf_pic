#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Vector, Constructors) {
  {
    Vector3d v;
    ASSERT_NEAR(0.0, v[0], 1e-8);
    ASSERT_NEAR(0.0, v[1], 1e-8);
    ASSERT_NEAR(0.0, v[2], 1e-8);
    std::cout << v;
  }
  {
    Vector2d v;
    ASSERT_NEAR(0.0, v[0], 1e-8);
    ASSERT_NEAR(0.0, v[1], 1e-8);
  }
  {
    Vector3d v({0.1, 0.2, 0.3});
    ASSERT_NEAR(0.1, v[0], 1e-8);
    ASSERT_NEAR(0.2, v[1], 1e-8);
    ASSERT_NEAR(0.3, v[2], 1e-8);
  }
  {
    Vector2d v({1.1, 1.2});
    ASSERT_NEAR(1.1, v[0], 1e-8);
    ASSERT_NEAR(1.2, v[1], 1e-8);
  }
  {
    Vector3d v(0.1, 0.2, 0.3);
    ASSERT_NEAR(0.1, v[0], 1e-8);
    ASSERT_NEAR(0.2, v[1], 1e-8);
    ASSERT_NEAR(0.3, v[2], 1e-8);
  }
  {
    Vector2d v(1.1, 1.2);
    ASSERT_NEAR(1.1, v[0], 1e-8);
    ASSERT_NEAR(1.2, v[1], 1e-8);
  }
  {
    Vector3i v;
    EXPECT_EQ(0, v[0]);
    EXPECT_EQ(0, v[1]);
    EXPECT_EQ(0, v[2]);
  }
  {
    Vector2i v;
    EXPECT_EQ(0, v[0]);
    EXPECT_EQ(0, v[1]);
  }
  {
    Vector3i v({1, 2, 3});
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
    EXPECT_EQ(3, v[2]);
  }
  {
    Vector2i v({1, 2});
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
  }
  {
    Vector3i v(1, 2, 3);
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
    EXPECT_EQ(3, v[2]);
  }
  {
    Vector2i v(1, 2);
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
  }

  Vector<int, 10> a({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  Vector<int, 10> b = a;
  for (int i = 0; i < 10; i++)
    EXPECT_EQ(a[i], b[i]);
}

TEST(Vector, Boolean) {
  Vector3i a(1, 2, 3);
  Vector3i b(1, 2, 3);
  Vector3i c(1, 0, 0);
  Vector3i d(10, 10, 10);
  EXPECT_EQ(a == b, true);
  EXPECT_EQ(a != b, false);
  EXPECT_EQ(a == c, false);
  EXPECT_EQ(a != c, true);
  EXPECT_EQ(c <= a, true);
  EXPECT_EQ(a >= c, true);
  EXPECT_EQ(d > a, true);
  EXPECT_EQ(b > a, false);
  EXPECT_EQ(a < d, true);
  EXPECT_EQ(c < a, false);
}

TEST(Vector, Arithmetic) {
  EXPECT_EQ(Vector2i(1, 1) - Vector2i(2, -3), Vector2i(-1, 4));
  EXPECT_EQ(Vector2i(0, 0) + Vector2i(2, -3), Vector2i(2, -3));
  EXPECT_EQ(Vector2i(1, 2) * Vector2i(2, -3), Vector2i(2, -6));
  EXPECT_EQ(Vector2i(1, 2) * 2, Vector2i(2, 4));
  EXPECT_EQ(Vector2i(10, 20) / Vector2i(2, -2), Vector2i(5, -10));
  Vector2d a(5.0, 10.0);
  EXPECT_EQ(a / 2.0, Vector2d(2.5, 5.0));
  a /= 2.0;
  EXPECT_EQ(a, Vector2d(2.5, 5.0));
  a -= 2.0;
  EXPECT_EQ(a, Vector2d(0.5, 3.0));
}

TEST(Vector, Gets) {
  EXPECT_EQ(Vector3i(1, 2, 3).xy(), Vector2i(1, 2));
  EXPECT_EQ(Vector3i(1, 2, 3).xy(1, 0), Vector2i(2, 1));
  Vector<double, 10> a({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
  EXPECT_EQ(a.xy(5, 9), Vector2d(6.0, 10.0));
  a[5] = 60.0;
  EXPECT_EQ(a.xy(0, 5), Vector2d(1.0, 60.0));
  EXPECT_EQ(Vector3i(1, 2, 3).doubleXY(), Vector2d(1.0, 2.0));
  EXPECT_EQ(Vector3i(1, 2, 3).doubleXY(1, 0), Vector2d(2.0, 1.0));
  Vector<char, 4> b({'a', 'b', 'c', 'd'});
  EXPECT_EQ(b.doubleXYZ(), Vector3d({97.0, 98.0, 99.0}));
  EXPECT_EQ(b.doubleXYZ(3, 0, 2), Vector3d({100.0, 97.0, 99.0}));
  EXPECT_EQ(b.x(), 'a');
  EXPECT_EQ(b.y(), 'b');
  EXPECT_EQ(b.z(), 'c');
}

TEST(Vector, Analytic) {
  Vector<double, 10> a(
      {1.0, 254.0, 3.0, 4.0, 50.0, 6.0, 7.0, -823.0, 9.0, 10.0});
  ASSERT_NEAR(a.max(), 254.0, 1e-8);
  ASSERT_NEAR(Vector2i(3, 4).length2(), 25.0, 1e-8);
  ASSERT_NEAR(Vector2i(3, 4).length(), 5.0, 1e-8);
  EXPECT_EQ(Vector2i(3, 4).normalized(), Vector2d(3.0, 4.0) / 5.0);
  EXPECT_EQ(Vector2i(3, 4).left(), Vector2i(-4, 3));
  EXPECT_EQ(Vector2i(3, 4).right(), Vector2i(4, -3));
}

TEST(Vector, Functions) {
  EXPECT_EQ(2 * Vector2i(1, 5), Vector2i(2, 10));
  EXPECT_EQ(round(Vector2d(-10.9, 2.01)), Vector2i(-11, 2));
  EXPECT_EQ(round(Vector3d(-10.499, 2.51, 1.0)), Vector3i(-10, 3, 1));
  EXPECT_EQ(ceil(Vector2d(-10.9, 2.01)), Vector2i(-9, 3));
  EXPECT_EQ(ceil(Vector3d(-10.499, 2.51, 1.0)), Vector3i(-9, 3, 2));
  EXPECT_EQ(floor(Vector2d(-10.9, 2.01)), Vector2i(-10, 2));
  EXPECT_EQ(floor(Vector3d(-10.499, 2.51, 1.0)), Vector3i(-10, 2, 1));
  EXPECT_EQ(min(Vector3i(0, 4, -1), Vector3i(-1, 8, 0)), Vector3i(-1, 4, -1));
  EXPECT_EQ(max(Vector3i(0, 4, -1), Vector3i(-1, 8, 0)), Vector3i(0, 8, 0));
  EXPECT_EQ(min(Vector2i(0, 4), Vector2i(-1, 8)), Vector2i(-1, 4));
  EXPECT_EQ(max(Vector2i(0, 4), Vector2i(-1, 8)), Vector2i(0, 8));
  Vector2d a(3.0, 4.0);
  normalize(a);
  EXPECT_EQ(a, Vector2d(3.0, 4.0) / 5.0);
}
