#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(QuadTree, Constructors) {
  {
    QuadTree qt;
    EXPECT_EQ(qt.nodeCount(), 0u);
    EXPECT_EQ(qt.height(), 0u);
  }
  {
    QuadTree qt(BBox2D::make_unit_bbox(), [](QuadTree::Node &node) -> bool {
      if (node.level() < 3)
        return true;
      return false;
    });
    EXPECT_EQ(qt.nodeCount(), 1u + 4u + 4u * 4u + 4u * 4u * 4u);
    EXPECT_EQ(qt.height(), 3u);
  }
  {
    QuadTree qt(BBox2D::make_unit_bbox(), 3);
    EXPECT_EQ(qt.nodeCount(), 1u + 4u + 4u * 4u + 4u * 4u * 4u);
    EXPECT_EQ(qt.height(), 3u);
  }
}

TEST(QuadTree, Traversal) {
  QuadTree qt(BBox2D::make_unit_bbox(), 4);
  size_t k = 0;
  qt.traverse([&k](QuadTree::Node &node) -> bool {
    // EXPECT_EQ(k, node.id);
    k++;
    double s = node.region().size(0);
    double es = 1. / pow(2., 1. * node.level());
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(1);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    return true;
  });
  EXPECT_EQ(k, qt.nodeCount());
  k = 0;
  qt.traverse([&k](const QuadTree::Node &node) -> bool {
    // EXPECT_EQ(k, node.id);
    k++;
    double s = node.region().size(0);
    double es = 1. / pow(2., 1. * node.level());
    EXPECT_EQ(IS_EQUAL(s, es), true);
    s = node.region().size(1);
    EXPECT_EQ(IS_EQUAL(s, es), true);
    return true;
  });
  EXPECT_EQ(k, qt.nodeCount());
}

TEST(QuadTree, Iterate) {
  QuadTree qt(BBox2D::make_unit_bbox(), 4);
  size_t k = 0;
  qt.iterateLeafs([&k](QuadTree::Leaf &l) {
    k++;
    EXPECT_EQ(l.edges.size(), 4u);
  });
  EXPECT_EQ(k, 4u * 4u * 4u * 4u);
}

TEST(QuadTree, Neighborhood) {
  {
    QuadTree qt(BBox2D::make_unit_bbox(), 4);
    auto leafs = qt.leafs();
    auto edges = qt.leafEdges();
    for (size_t i = 0; i < leafs.size(); i++) {
      EXPECT_EQ(leafs[i].edges.size(), 4u);
      for (size_t e = 0; e < leafs[i].edges.size(); e++) {
        auto &edge = edges[leafs[i].edges[e]];
        EXPECT_EQ(edge.neighborIndex(i), qt.neighborIndex(i, e));
        // ASSERT_NEAR(edge.neighborWeight(i), 0., 1e-8);
        // edge.neighborWeight(i) = 1. * e;
        // ASSERT_NEAR(edge.neighborWeight(i), 1. * e, 1e-8);
      }
    }
  }
  {
    QuadTree qt(BBox2D::make_unit_bbox(), 1);
    auto leafs = qt.leafs();
    EXPECT_EQ(leafs.size(), 4u);
  }
}

TEST(QuadTree, faceCenterPosition) {
  Point2d expectedPosition[12] = {
      Point2d(0.5, 0.75), Point2d(0.5, 0.25), Point2d(0.25, 0.5),
      Point2d(0.75, 0.5), Point2d(0.25, 0.),  Point2d(0.75, 0.),
      Point2d(0.25, 1.),  Point2d(0.75, 1.),  Point2d(1., 0.75),
      Point2d(1., 0.25),  Point2d(0., 0.75),  Point2d(0., 0.25)};
  QuadTree qt(BBox2D::make_unit_bbox(), 1);
  auto &edges = qt.leafEdges();
  for (size_t i = 0; i < edges.size(); i++)
    EXPECT_EQ(edges[i].faceCenterPosition == expectedPosition[i], true);
}

TEST(QuadTree, intersect) {
  QuadTree qt(BBox2D::make_unit_bbox(), 1);
  EXPECT_EQ(qt.intersect(Point2d(-1., 1.)), -1);
  EXPECT_EQ(qt.intersect(Point2d(0.75, 0.75)), 1);
  EXPECT_EQ(qt.intersect(Point2d(0.25, 0.75)), 0);
  EXPECT_EQ(qt.intersect(Point2d(0.75, 0.25)), 3);
  EXPECT_EQ(qt.intersect(Point2d(0.25, 0.25)), 2);
}
