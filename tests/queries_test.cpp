#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Queries, BBoxBBox) {
  EXPECT_EQ(bbox_bbox_intersection(BBox2D::make_unit_bbox(),
                                   BBox2D(Point2d(.5, .0), Point2d(.75, 2.))),
            true);
}
