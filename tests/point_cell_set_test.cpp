#include <gtest/gtest.h>
#include <furoo.h>

using namespace furoo;

TEST(PointCellSet, Access) {
  PointCellSet<double,2> s;
  s.add(Point2d(0.5, 0.5));
  s.add(Point2d(1.5, 1.5));
  s.add(Point2d(2.5, 2.5));
  s.add(Point2d(3.5, 3.5));
  s.add(Point2d(4.5, 4.5));
  s.setCell(0, 0);
  s.setCell(1, 0);
  s.setCell(2, 0);
  s.setCell(3, 0);
  s.setCell(4, 0);
}