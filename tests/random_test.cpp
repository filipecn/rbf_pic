#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;
TEST(HaltonSequence, randomDouble) {
  {
    HaltonSequence rng;
    std::vector<double> s;
    for (size_t i = 0; i < 1000; i++)
      s.emplace_back(rng.randomDouble());
    for (size_t i = 0; i < 1000; i++)
      for (size_t j = i + 1; j < 1000; j++)
        EXPECT_EQ(s[i] != s[j], true);
  }
  {
    HaltonSequence rng;
    std::vector<double> s;
    for (size_t i = 0; i < 1000; i++)
      s.emplace_back(rng.randomDouble(-10., 20.));
    for (size_t i = 0; i < 1000; i++) {
      EXPECT_EQ(-10. <= s[i] && s[i] <= 20., true);
      for (size_t j = i + 1; j < 1000; j++)
        EXPECT_EQ(s[i] != s[j], true);
    }
  }
}

TEST(BoxSampler, sample) {
  {
    BoxSampler sampler;
    std::vector<Point2d> s;
    BBox2D region(Point2d(-1., 3.), Point2d(2., 4.));
    for (size_t i = 0; i < 1000; i++)
      s.emplace_back(sampler.sample(region));
    for (size_t i = 0; i < 1000; i++) {
      EXPECT_EQ(region.inside(s[i]), true);
      for (size_t j = i + 1; j < 1000; j++)
        EXPECT_EQ(s[i] != s[j], true);
    }
  }
  {
    BoxSampler sampler;
    std::vector<Point2d> s;
    BBox2D region = BBox2D::make_unit_bbox();
    for (size_t i = 0; i < 1000; i++)
      s.emplace_back(sampler.sample());
    for (size_t i = 0; i < 1000; i++) {
      EXPECT_EQ(region.inside(s[i]), true);
      for (size_t j = i + 1; j < 1000; j++)
        EXPECT_EQ(s[i] != s[j], true);
    }
  }
}
