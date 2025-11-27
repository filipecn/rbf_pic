#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Octree, Constructors) {
  {
    Octree qt;
    EXPECT_EQ(qt.nodeCount(), 0u);
    EXPECT_EQ(qt.height(), 0u);
  }
  {
    Octree qt(BBox3D::make_unit_bbox(), [](Octree::Node &node) -> bool {
      if (node.level() < 3)
        return true;
      return false;
    });
    EXPECT_EQ(qt.nodeCount(), 1u + 8u + 8u * 8u + 8u * 8u * 8u);
    EXPECT_EQ(qt.height(), 3u);
  }
  {
    Octree qt(BBox3D::make_unit_bbox(), 3);
    EXPECT_EQ(qt.nodeCount(), 1u + 8u + 8u * 8u + 8u * 8u * 8u);
    EXPECT_EQ(qt.height(), 3u);
  }
}

TEST(Octree, Traversal) {
  Octree qt(BBox3D::make_unit_bbox(), 4);
  size_t k = 0;
  qt.traverse([&k](Octree::Node &node) -> bool {
    k++;
    double s = node.region().size(0);
    double es = 1. / pow(2., 1. * node.level());
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(1);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(2);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    return true;
  });
  EXPECT_EQ(k, qt.nodeCount());
  k = 0;
  qt.traverse([&k](const Octree::Node &node) -> bool {
    k++;
    double s = node.region().size(0);
    double es = 1. / pow(2., 1. * node.level());
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(1);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(2);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    return true;
  });
  EXPECT_EQ(k, qt.nodeCount());
}

TEST(Octree, Iterate) {
  Octree qt(BBox3D::make_unit_bbox(), 4);
  size_t k = 0;
  qt.iterateLeafs([&k](Octree::Leaf &l) {
    k++;
    EXPECT_EQ(l.edges.size(), 6u);
  });
  EXPECT_EQ(k, 8u * 8u * 8u * 8u);
}

TEST(Octree, Neighborhood) {
  {
    Octree qt(BBox3D::make_unit_bbox(), 4);
    auto leafs = qt.leafs();
    auto edges = qt.leafEdges();
    for (size_t i = 0; i < leafs.size(); i++) {
      EXPECT_EQ(leafs[i].edges.size(), 6u);
      for (size_t e = 0; e < leafs[i].edges.size(); e++) {
        auto &edge = edges[leafs[i].edges[e]];
        EXPECT_EQ(edge.neighborIndex(i), qt.neighborIndex(i, e));
      }
    }
  }
  {
    Octree qt(BBox3D::make_unit_bbox(), 1);
    auto leafs = qt.leafs();
    EXPECT_EQ(leafs.size(), 8u);
  }
}

TEST(Octree, faceCenterPosition) {
  Point3d expectedPosition[36] = {
      Point3d(0.5, 0.25, 0.25), // 0 - 1
      Point3d(0.5, 0.75, 0.25), // 2 - 3
      Point3d(0.5, 0.25, 0.75), // 4 - 5
      Point3d(0.5, 0.75, 0.75), // 6 - 7
      Point3d(0.25, 0.5, 0.25), // 0 - 2
      Point3d(0.75, 0.5, 0.25), // 1 - 3
      Point3d(0.25, 0.5, 0.75), // 4 - 6
      Point3d(0.75, 0.5, 0.75), // 5 - 7
      Point3d(0.25, 0.25, 0.5), // 0 - 4
      Point3d(0.75, 0.25, 0.5), // 1 - 5
      Point3d(0.25, 0.75, 0.5), // 2 - 6
      Point3d(0.75, 0.75, 0.5), // 3 - 7
      Point3d(1, 0.25, 0.25),   // 1 - x
      Point3d(1, 0.75, 0.25),   // 3 - x
      Point3d(1, 0.25, 0.75),   // 5 - x
      Point3d(1, 0.75, 0.75),   // 7 - x
      Point3d(0, 0.25, 0.25),   // x - 0
      Point3d(0, 0.75, 0.25),   // x - 2
      Point3d(0, 0.25, 0.75),   // x - 4
      Point3d(0, 0.75, 0.75),   // x - 6
      Point3d(0.25, 0, 0.25),   // 0 - x
      Point3d(0.75, 0, 0.25),   // 1 - x
      Point3d(0.25, 0, 0.75),   // 4 - x
      Point3d(0.75, 0, 0.75),   // 5 - x
      Point3d(0.25, 1, 0.25),   // x - 2
      Point3d(0.75, 1, 0.25),   // x - 3
      Point3d(0.25, 1, 0.75),   // x - 6
      Point3d(0.75, 1, 0.75),   // x - 7
      Point3d(0.25, 0.25, 0),   // x - 0
      Point3d(0.75, 0.25, 0),   // x - 2
      Point3d(0.25, 0.75, 0),   // x - 1
      Point3d(0.75, 0.75, 0),   // x - 3
      Point3d(0.25, 0.25, 1),   // 4 - x
      Point3d(0.75, 0.25, 1),   // 5 - x
      Point3d(0.25, 0.75, 1),   // 6 - x
      Point3d(0.75, 0.75, 1)    // 7 - x
  };
  Octree qt(BBox3D::make_unit_bbox(), 1);
  auto &edges = qt.leafEdges();
  for (size_t i = 0; i < edges.size(); i++)
    EXPECT_EQ(edges[i].faceCenterPosition == expectedPosition[i], true);
}

TEST(Octree, intersect) {
  Octree qt(BBox3D::make_unit_bbox(), 1);
  EXPECT_EQ(qt.intersect(Point3d(-1., 1., 0.)), -1);
  EXPECT_EQ(qt.intersect(Point3d(0.25, 0.25, 0.25)), 0);
  EXPECT_EQ(qt.intersect(Point3d(0.75, 0.25, 0.25)), 1);
  EXPECT_EQ(qt.intersect(Point3d(0.25, 0.75, 0.25)), 2);
  EXPECT_EQ(qt.intersect(Point3d(0.75, 0.75, 0.25)), 3);
  EXPECT_EQ(qt.intersect(Point3d(0.25, 0.25, 0.75)), 4);
  EXPECT_EQ(qt.intersect(Point3d(0.75, 0.25, 0.75)), 5);
  EXPECT_EQ(qt.intersect(Point3d(0.25, 0.75, 0.75)), 6);
  EXPECT_EQ(qt.intersect(Point3d(0.75, 0.75, 0.75)), 7);
}
