#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(LinearVector, Operators) {
  LinearVector lv = std::vector<double>(10, 2.0);
  EXPECT_EQ(lv.size() == 10u, true);
  for (size_t i = 0; i < lv.size(); i++) {
    ASSERT_NEAR(lv[i], 2.0, 1e-8);
    lv[i] = 1. * i;
    ASSERT_NEAR(lv[i], 1. * i, 1e-8);
  }
  LinearVector lv2(std::vector<double>(100, 1.0));
  EXPECT_EQ(lv2.size() == 100u, true);
  lv2.resize(1000, 2.0);
  EXPECT_EQ(lv2.size() == 1000u, true);
  for (size_t i = 0; i < lv2.size(); i++) {
    if (i < 100u)
      ASSERT_NEAR(lv2[i], 1., 1e-8);
    else
      ASSERT_NEAR(lv2[i], 2., 1e-8);
  }
}
