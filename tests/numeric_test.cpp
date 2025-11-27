#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Numeric, Macros) {
  EXPECT_EQ(SQR(8), 64);
  EXPECT_EQ(CUBE(3), 27);
  ASSERT_NEAR(TO_DEGREES(PI_2), 360.0, 1e-8);
  ASSERT_NEAR(TO_RADIANS(90.0), PI / 2.0, 1e-8);
  EXPECT_TRUE(IS_ZERO(1e-9));
  EXPECT_TRUE(IS_NOT_ZERO(2e-8));
  EXPECT_TRUE(IS_EQUAL(PI, PI + 1e-8));
  EXPECT_TRUE(IS_EQUAL_ERROR(1e-5, 1e-6, 1e-5));
  EXPECT_TRUE(IS_BETWEEN(2, 1, 3));
  EXPECT_FALSE(IS_BETWEEN(1, 1, 3));
  EXPECT_TRUE(IS_BETWEEN_CLOSE(2, 1, 3));
  EXPECT_FALSE(IS_BETWEEN_CLOSE(4, 1, 3));
  EXPECT_TRUE(IS_BETWEEN_CLOSE(1, 1, 3));
}

TEST(Numeric, sign) {
  EXPECT_EQ(sign(-1), -1);
  EXPECT_EQ(sign(0), 1);
  EXPECT_EQ(sign(1), 1);
  ASSERT_NEAR(sign(-1.), -1., 1e-8);
  ASSERT_NEAR(sign(0.), 1., 1e-8);
  ASSERT_NEAR(sign(1.), 1., 1e-8);
}

TEST(Numeric, IsPowerOf2) {
  EXPECT_FALSE(isPowerOf2(15));
  EXPECT_FALSE(isPowerOf2(100));
  for (int i = 0; i < 30; i++)
    EXPECT_TRUE(isPowerOf2(2 << i));
}

TEST(Numeric, linspace) {
  auto v = linspace(1., 16., 16);
  for (size_t i = 1; i <= 16; i++)
    EXPECT_EQ(static_cast<size_t>(v[i - 1]), i);
}
TEST(Numeric, Clamp) {
  EXPECT_EQ(clamp(7, 1, 5), 5);
  EXPECT_EQ(clamp(-5, 1, 5), 1);
  ASSERT_NEAR(clamp(0.3, 0.5, 0.8), 0.5, 1e-8);
  ASSERT_NEAR(clamp(10., 0.5, 1.8), 1.8, 1e-8);
}

TEST(Numeric, SeparateBy1) { EXPECT_EQ(separateBy1(5), 17u); }

TEST(Numeric, MortonCode) {
  EXPECT_EQ(mortonCode((unsigned)1.9, (unsigned)2.1), 9u);
}

TEST(Numeric, lerp) {
  ASSERT_NEAR(lerp(.5, 0., 10.), 5., 1e-8);
  ASSERT_NEAR(lerp(0., -5., 5.), -5., 1e-8);
  ASSERT_NEAR(lerp(1., -5., 5.), 5., 1e-8);
  ASSERT_NEAR(lerp(.5, -5., 5.), 0., 1e-8);
}

TEST(Numeric, rounding) {
  EXPECT_EQ(floor2Int(1e-7), 0);
  EXPECT_EQ(floor2Int(1.887), 1);
  EXPECT_EQ(floor2Int(3.59494), 3);
  EXPECT_EQ(floor2Int(77 - 1e-8), 76);
  EXPECT_EQ(floor2Int(-0.5), -1);

  EXPECT_EQ(ceil2Int(1e-7), 1);
  EXPECT_EQ(ceil2Int(1.887), 2);
  EXPECT_EQ(ceil2Int(3.59494), 4);
  EXPECT_EQ(ceil2Int(77 - 1e-8), 77);
  EXPECT_EQ(ceil2Int(-0.5), 0);

  EXPECT_EQ(round2Int(1e-7), 0);
  EXPECT_EQ(round2Int(1.887), 2);
  EXPECT_EQ(round2Int(3.5000000001), 4);
  EXPECT_EQ(round2Int(3.4999999999), 3);
  EXPECT_EQ(round2Int(77 - 1e-8), 77);
}

TEST(Numeric, maxMin) {
  EXPECT_EQ(max(1e-7, 2e-7), 2e-7);
  EXPECT_EQ(max(1.887, 1.88), 1.887);
  EXPECT_EQ(max(-1, -5), -1);
  EXPECT_EQ(max(-997123, -997122), -997122);

  EXPECT_EQ(min(1e-7, 2e-7), 1e-7);
  EXPECT_EQ(min(1.887, 1.88), 1.88);
  EXPECT_EQ(min(-1, -5), -5);
  EXPECT_EQ(min(-997123, -997122), -997123);
}

TEST(Numeric, smoothStep) {
  ASSERT_NEAR(smoothStep(0.5, 0., 1.), 0.5, 1e-8);
  ASSERT_NEAR(smoothStep(5., 0., 10.), 0.5, 1e-8);
}

TEST(Numeric, linearStep) {
  ASSERT_NEAR(linearStep(0.5, 0., 1.), 0.5, 1e-8);
  ASSERT_NEAR(linearStep(5., 0., 10.), 0.5, 1e-8);
  ASSERT_NEAR(linearStep(2.5, 0., 10.), 0.25, 1e-8);
  ASSERT_NEAR(linearStep(75., 0., 100.), 0.75, 1e-8);
}

TEST(Numeric, mod) {
  EXPECT_EQ(mod(5., 2), 1);
  EXPECT_EQ(mod(5., 1), 0);
  EXPECT_EQ(mod(100, 7), 2);
  EXPECT_EQ(mod(-7, 7), 0);
  EXPECT_EQ(mod(4, -7), 4 % (-7));
  // EXPECT_EQ(mod(4, .75), 0.25);
}

TEST(Numeric, roundUpPow2) {
  EXPECT_EQ(roundUpPow2(5), 8u);
  EXPECT_EQ(roundUpPow2(1), 1u);
  EXPECT_EQ(roundUpPow2(10), 16u);
  EXPECT_EQ(roundUpPow2(7), 8u);
  EXPECT_EQ(roundUpPow2(129), 256u);
}

TEST(Numeric, log2Int) {
  EXPECT_EQ(log2Int(5), 2);
  EXPECT_EQ(log2Int(1), 0);
  EXPECT_EQ(log2Int(0.1), -4);
  EXPECT_EQ(log2Int(8.512223), 3);
  EXPECT_EQ(log2Int(1025), 10);
}

TEST(Numeric, SolveQuadratic) {
  double x1, x2;
  EXPECT_EQ(solve_quadratic(1., 2., 1., x1, x2), true);
  ASSERT_NEAR(x1, -1, 1e-8);
  ASSERT_NEAR(x2, -1, 1e-8);
  EXPECT_EQ(solve_quadratic(1., -2., 1., x1, x2), true);
  ASSERT_NEAR(x1, 1, 1e-8);
  ASSERT_NEAR(x2, 1, 1e-8);
  EXPECT_EQ(solve_quadratic(1., 0., -1., x1, x2), true);
  ASSERT_NEAR(x1, -1, 1e-8);
  ASSERT_NEAR(x2, 1, 1e-8);
  EXPECT_EQ(solve_quadratic(1., 1., -1., x1, x2), true);
  ASSERT_NEAR(x1, -1.618033988749895, 1e-8);
  ASSERT_NEAR(x2, 0.618033988749895, 1e-8);
  EXPECT_EQ(solve_quadratic(11., 3.447, -41., x1, x2), true);
  ASSERT_NEAR(x1, -2.09364385173371, 1e-8);
  ASSERT_NEAR(x2, 1.78028021537008, 1e-8);
  EXPECT_EQ(solve_quadratic(1, 0, 1., x1, x2), false);
}

TEST(Numeric, CatmullRomSpline) {
  ASSERT_NEAR(catmullRomSpline(-1.5, 0.5, 1.5, 2.5, 0.), 0.5, 1e-2);
  ASSERT_NEAR(catmullRomSpline(-1.5, 0.5, 1.5, 2.5, 1.), 1.5, 1e-2);
  ASSERT_NEAR(catmullRomSpline(-1.5, 0.5, 1.0, 2.5, 0.), 0.5, 1e-2);
  ASSERT_NEAR(catmullRomSpline(-1.5, 0.5, 1.0, 2.5, 1.), 1.0, 1e-2);
  ASSERT_NEAR(catmullRomSpline(std::sin(0.1), std::sin(0.2), std::sin(0.3),
                               std::sin(0.4), 0.),
              std::sin(0.2), 1e-3);
  ASSERT_NEAR(catmullRomSpline(std::sin(0.1), std::sin(0.2), std::sin(0.3),
                               std::sin(0.4), 1.),
              std::sin(0.3), 1e-3);
  ASSERT_NEAR(catmullRomSpline(std::sin(0.1), std::sin(0.2), std::sin(0.3),
                               std::sin(0.4), .5),
              std::sin(0.25), 1e-3);
  ASSERT_NEAR(catmullRomSpline(std::sin(0.1), std::sin(0.2), std::sin(0.3),
                               std::sin(0.4), .75),
              std::sin(0.275), 1e-3);
  ASSERT_NEAR(catmullRomSpline(-100., 0., 100., 200, 0.5), 50., 1e-3);
  ASSERT_NEAR(catmullRomSpline(1., 4., 9., 16., 0.), 4., 1e-3);
  ASSERT_NEAR(catmullRomSpline(1., 4., 9., 16., 1.), 9., 1e-3);
  ASSERT_NEAR(catmullRomSpline(1., 4., 9., 16., .5), SQR(2.5), 1e-4);
  HaltonSequence rng;
  for (size_t i = 0; i < 1000; i++) {
    double s = rng.randomDouble();
    ASSERT_NEAR(catmullRomSpline(1., 4., 9., 16., s), SQR(2. + s), 1e-5);
  }
  for (size_t i = 0; i < 1000; i++) {
    double s = rng.randomDouble();
    ASSERT_NEAR(catmullRomSpline(std::sin(0.1), std::sin(0.2), std::sin(0.3),
                                 std::sin(0.4), s),
                std::sin(0.2 + 0.1 * s), 1e-4);
  }
}
