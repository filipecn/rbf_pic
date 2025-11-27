#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(TrilinearHatKernel, Operator) {
  TrilinearHatKernel<Point2d> k;
  ASSERT_NEAR(k.phi(Point2d(), Point2d(0.5)), .25, 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(1.5)), 0., 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, 1.5)), 0., 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, 0.3)), 0.35, 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, -0.3)), 0.35, 1e-8);
}

TEST(QuadraticBSplineKernel, Operator) {
  QuadraticBSplineKernel<Point2d> k;
  ASSERT_NEAR(k.phi(Point2d(), Point2d(0.5)), .25, 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(1.5)), 0., 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, 1.5)), 0., 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, 0.3)), 0.33, 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(.5, -0.3)), 0.33, 1e-8);
  ASSERT_NEAR(k.phi(Point2d(), Point2d(-1.3, -0.3)), 0.0132, 1e-8);
}
