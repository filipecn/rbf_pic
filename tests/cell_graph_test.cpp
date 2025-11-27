#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(CellGraph2, RefineOneCell) {
  { // refine 1
    CellGraph2 graph(BBox2d::squareBox());
    auto children = graph.refine(1);
    EXPECT_EQ(graph.nodeCount(), 4u);
    EXPECT_EQ(graph.cellCount(), 4u);
    for (size_t i = 0; i < 4; i++)
      EXPECT_EQ(children[i], i + 1);
    {
      auto neighbors = graph.neighbors(1);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[0].id, 2);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[1].id, 3);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(2);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[0].id, 1);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 4);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(3);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[1].id, 1);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 4);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(4);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[1].id, 2);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 3);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
    }
  }
  { // REFINE 1 AND 1
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    auto children = graph.refine(1);
    EXPECT_EQ(graph.nodeCount(), 7u);
    EXPECT_EQ(children[0], 1u);
    EXPECT_EQ(children[1], 5u);
    EXPECT_EQ(children[2], 6u);
    EXPECT_EQ(children[3], 7u);
    {
      auto neighbors = graph.neighbors(5);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[0].id, 1);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 7);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[2].id, 2);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(2);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[3].id, 5);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[4].id, 7);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[0].id, 4);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(3);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[3].id, 6);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[4].id, 7);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 4);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(6);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[1].id, 1);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[2].id, 3);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 7);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(7);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[1].id, 5);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 6);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 2);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[3].id, 3);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
    }
  }
  { // REFINE 1 AND 2
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    auto children = graph.refine(2);
    EXPECT_EQ(graph.nodeCount(), 7u);
    EXPECT_EQ(children[0], 2u);
    EXPECT_EQ(children[1], 5u);
    EXPECT_EQ(children[2], 6u);
    EXPECT_EQ(children[3], 7u);
    {
      auto neighbors = graph.neighbors(1);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[3].id, 2);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[4].id, 6);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[0].id, 3);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(2);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[2].id, 1);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 6);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 5);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(4);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 3);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[3].id, 6);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[4].id, 7);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::TOP);
    }
    {
      auto neighbors = graph.neighbors(6);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[2].id, 1);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 2);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[3].id, 4);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 7);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(7);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[1].id, 5);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 6);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 4);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
    }
  }
  { // REFINE 1 AND 3
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    auto children = graph.refine(3);
    EXPECT_EQ(graph.nodeCount(), 7u);
    EXPECT_EQ(children[0], 3u);
    EXPECT_EQ(children[1], 5u);
    EXPECT_EQ(children[2], 6u);
    EXPECT_EQ(children[3], 7u);
    {
      auto neighbors = graph.neighbors(1);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 2);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[4].id, 5);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[3].id, 3);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(3);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 1);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[1].id, 6);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 5);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(4);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[0].id, 2);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[3].id, 5);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[4].id, 7);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::LEFT);
    }
    {
      auto neighbors = graph.neighbors(5);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 1);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 3);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 4);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[1].id, 7);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(7);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[1].id, 5);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 6);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 4);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::RIGHT);
    }
  }
  { // REFINE 1 AND 4
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    auto children = graph.refine(4);
    EXPECT_EQ(graph.nodeCount(), 7u);
    EXPECT_EQ(children[0], 4u);
    EXPECT_EQ(children[1], 5u);
    EXPECT_EQ(children[2], 6u);
    EXPECT_EQ(children[3], 7u);
    {
      auto neighbors = graph.neighbors(2);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 1);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[4].id, 5);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[3].id, 4);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(3);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[0].id, 1);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[4].id, 6);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[3].id, 4);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::RIGHT);
    }
    {
      auto neighbors = graph.neighbors(4);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 2);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 5);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 3);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 6);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(5);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[2].id, 2);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[0].id, 4);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 7);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::BOTTOM);
    }
    {
      auto neighbors = graph.neighbors(6);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[3].id, 0);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
      EXPECT_EQ(neighbors[1].id, 4);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[2].id, 3);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[0].id, 7);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
    }
  }
  {
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    graph.refine(1);
    graph.refine(7);
    EXPECT_EQ(graph.nodeCount(), 10u);
    {
      auto neighbors = graph.neighbors(2);
      EXPECT_EQ(neighbors.size(), 6u);
      EXPECT_EQ(neighbors[2].id, 0);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[1].id, 0);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[3].id, 5);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[4].id, 8);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[5].id, 10);
      EXPECT_EQ(neighbors[5].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[0].id, 4);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::BOTTOM);
    }
    graph.refine(5);
    {
      auto neighbors = graph.neighbors(13);
      EXPECT_EQ(neighbors.size(), 4u);
      EXPECT_EQ(neighbors[0].id, 12);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[1].id, 11);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[2].id, 2);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[3].id, 8);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::BOTTOM);
    }
    graph.refine(2);
    {
      auto neighbors = graph.neighbors(15);
      EXPECT_EQ(neighbors.size(), 5u);
      EXPECT_EQ(neighbors[0].id, 16);
      EXPECT_EQ(neighbors[0].side, Definitions::Side::RIGHT);
      EXPECT_EQ(neighbors[1].id, 2);
      EXPECT_EQ(neighbors[1].side, Definitions::Side::TOP);
      EXPECT_EQ(neighbors[3].id, 8);
      EXPECT_EQ(neighbors[3].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[4].id, 10);
      EXPECT_EQ(neighbors[4].side, Definitions::Side::LEFT);
      EXPECT_EQ(neighbors[2].id, 4);
      EXPECT_EQ(neighbors[2].side, Definitions::Side::BOTTOM);
    }
  }
}

TEST(CellGraph2, Siblings) {
  CellGraph2 graph(BBox2d::squareBox());
  graph.refine(1);
  {
    auto siblings = graph.nodeSiblings(1);
    EXPECT_EQ(siblings.size(), 4u);
    EXPECT_EQ(siblings[0], 1u);
    EXPECT_EQ(siblings[1], 2u);
    EXPECT_EQ(siblings[2], 3u);
    EXPECT_EQ(siblings[3], 4u);
    EXPECT_TRUE(graph.canCoarse({4, 3, 1, 2}));
  }
  graph.refine(2);
  {
    auto siblings = graph.nodeSiblings(1);
    EXPECT_EQ(siblings.size(), 3u);
    EXPECT_EQ(siblings[0], 1u);
    EXPECT_EQ(siblings[1], 3u);
    EXPECT_EQ(siblings[2], 4u);
  }
  {
    auto siblings = graph.nodeSiblings(2);
    EXPECT_EQ(siblings.size(), 4u);
    EXPECT_EQ(siblings[0], 2u);
    EXPECT_EQ(siblings[1], 5u);
    EXPECT_EQ(siblings[2], 6u);
    EXPECT_EQ(siblings[3], 7u);
    EXPECT_TRUE(graph.canCoarse({2, 5, 6, 7}));
  }
  graph.refine(4);
  {
    auto siblings = graph.nodeSiblings(4);
    EXPECT_EQ(siblings.size(), 4u);
    EXPECT_EQ(siblings[0], 4u);
    EXPECT_EQ(siblings[1], 8u);
    EXPECT_EQ(siblings[2], 9u);
    EXPECT_EQ(siblings[3], 10u);
    EXPECT_TRUE(graph.canCoarse({4, 8, 9, 10}));
    EXPECT_FALSE(graph.canCoarse({4, 8, 9, 7}));
    EXPECT_FALSE(graph.canCoarse({4, 8, 9, 6}));
    EXPECT_FALSE(graph.canCoarse({4, 8, 9, 3}));
  }
}

TEST(CellGraph2, Coarse) {
  {
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    graph.coarse(1);
    EXPECT_EQ(graph.nodeCount(), 1u);
    {
      auto n = graph.neighbors(1);
      EXPECT_EQ(n.size(), 4u);
      for (size_t i = 0; i < 4u; i++)
        EXPECT_EQ(n[i].id, 0);
    }
  }
  {
    CellGraph2 graph(BBox2d::squareBox());
    graph.refine(1);
    graph.refine(2);
    graph.refine(4);
    graph.refine(6);
    graph.coarse(12);
    EXPECT_EQ(graph.nodeCount(), 10u);
    {
      auto n = graph.neighbors(12);
      EXPECT_EQ(n.size(), 4u);
      EXPECT_EQ(n[0].id, 1);
      EXPECT_EQ(n[1].id, 2);
      EXPECT_EQ(n[2].id, 4);
      EXPECT_EQ(n[3].id, 7);
    }
  }
}

TEST(CellGraph2, Refine) {
  CellGraph2 graph(BBox2d::squareBox());
  graph.refine(1, [](const CellGraph2::Node &node) -> bool {
    return node.level() < 3;
  });
  EXPECT_EQ(graph.nodeCount(), 64u);
  EXPECT_EQ(graph.edgeCount(), 144u);
}

TEST(CellGraph2, CellId) {
  {
    CellGraph2 graph(BBox2d::squareBox(10));
    graph.refine(1);
    graph.refine(2);
    graph.refine(4);
    graph.refine(6);
    graph.refine(8);
    EXPECT_EQ(graph.cellId(Point2d(2, 2)), 3);
    EXPECT_EQ(graph.cellId(Point2d(6, 6)), 12);
    EXPECT_EQ(graph.cellId(Point2d(6, 4)), 4);
    EXPECT_EQ(graph.cellId(Point2d(8, 6)), 7);
    EXPECT_EQ(graph.cellId(Point2d(9, 7)), 7);
    EXPECT_EQ(graph.cellId(Point2d(-1)), -1);
  }
}

TEST(cellGraph2, RefineCoarse) {
  CellGraph2 cellGraph2(BBox2d::squareBox());
  cellGraph2.refine(1);
  cellGraph2.coarse(3);
  cellGraph2.refine(3);
  cellGraph2.cellId(Point2d(0.25, 0.25));
}

TEST(cellGraph2, IterateStar) {
  CellGraph2 graph(BBox2d::squareBox(1));
  graph.refine(1);
  graph.refine(2);
  graph.refine(7);
  graph.refine(9);
  size_t count = 0;
  graph.iterateNeighborFaces(Point2d(0.75, 0.25),
                             0, Definitions::Orientation::HORIZONTAL,
                             [&](size_t id) {
                               UNUSED_VARIABLE(id);
                               count++;
                             });
  EXPECT_EQ(count, 5u);
  count = 0;
  graph.iterateNeighborFaces(Point2d(0.75, 0.25),
                             0, Definitions::Orientation::VERTICAL,
                             [&](size_t id) {
                               UNUSED_VARIABLE(id);
                               count++;
                             });
  EXPECT_EQ(count, 2u);
  count = 0;
  graph.iterateNeighborFaces(Point2d(0.75, 0.25),
                             1, Definitions::Orientation::HORIZONTAL,
                             [&](size_t id) {
                               UNUSED_VARIABLE(id);
                               count++;
                             });
  EXPECT_EQ(count, 11u);
  count = 0;
  graph.iterateNeighborFaces(Point2d(0.75, 0.25),
                             1, Definitions::Orientation::VERTICAL,
                             [&](size_t id) {
                               UNUSED_VARIABLE(id);
                               count++;
                             });
  EXPECT_EQ(count, 11u);
}

TEST(cellGraph2, IterateCellNeighborhood) {
  CellGraph2 graph(BBox2d::squareBox());
  graph.refine(1);
  graph.refine(1);
  graph.refine(4);
  graph.iterateNeighborCells(7,
                             2,
                             [](size_t i) { std::cerr << i << std::endl; });

}