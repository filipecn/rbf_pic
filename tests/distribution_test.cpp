#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Distribution, linear) {
  {
    BBox2d box = BBox2d::squareBox(1.);
    auto points = Distribution<2>::linear(box, Definitions::Side::LEFT, 0);
    EXPECT_EQ(points.size(), 2u);
    EXPECT_EQ(points[0], Point2d(0., 0.));
    EXPECT_EQ(points[1], Point2d(0., 1.));
    points = Distribution<2>::linear(box, Definitions::Side::RIGHT, 0);
    EXPECT_EQ(points.size(), 2u);
    EXPECT_EQ(points[0], Point2d(1., 0.));
    EXPECT_EQ(points[1], Point2d(1., 1.));
    points = Distribution<2>::linear(box, Definitions::Side::BOTTOM, 0);
    EXPECT_EQ(points.size(), 2u);
    EXPECT_EQ(points[0], Point2d(0., 0.));
    EXPECT_EQ(points[1], Point2d(1., 0.));
    points = Distribution<2>::linear(box, Definitions::Side::TOP, 0);
    EXPECT_EQ(points.size(), 2u);
    EXPECT_EQ(points[0], Point2d(0., 1.));
    EXPECT_EQ(points[1], Point2d(1., 1.));
    points = Distribution<2>::linear(box, Definitions::Side::LEFT, 4);
    EXPECT_EQ(points.size(), 6u);
    for(size_t i = 0; i < 6; i++)
      EXPECT_EQ(points[i], Point2d(0., i * 0.2));
    points = Distribution<2>::linear(box, Definitions::Side::RIGHT, 4);
    EXPECT_EQ(points.size(), 6u);
    for(size_t i = 0; i < 6; i++)
      EXPECT_EQ(points[i], Point2d(1., i * 0.2));
    points = Distribution<2>::linear(box, Definitions::Side::BOTTOM, 4);
    EXPECT_EQ(points.size(), 6u);
    for(size_t i = 0; i < 6; i++)
      EXPECT_EQ(points[i], Point2d(i * 0.2, 0.));
    points = Distribution<2>::linear(box, Definitions::Side::TOP, 4);
    EXPECT_EQ(points.size(), 6u);
    for(size_t i = 0; i < 6; i++)
      EXPECT_EQ(points[i], Point2d(i * 0.2, 1.));
  }
}
