#include <fstream>
#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  int ret = RUN_ALL_TESTS();

  return ret;
}
