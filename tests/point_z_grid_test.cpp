#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(PointZGrid, SearchZ) {
  {
    PointZGrid2d zgrid(16u);
    BBox2d region(Point2d(0., 0.), Point2d(1, 1));
    zgrid.setDomainRegion(region);
    std::vector<Point2d> points;
    for (size_t i = 0; i < 16; i++)
      for (size_t j = 0; j < 16; j++) {
        Point2d p(i / 16. + 1. / 32., j / 16. + 1. / 32.);
        zgrid.add(p);
        points.emplace_back(p);
      }
    EXPECT_EQ(zgrid.size(), 16u * 16u);
    for (size_t i = 0; i < 16u * 16u; i++)
      EXPECT_EQ(points[i], zgrid[i]);
    std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
      if (a.x() < b.x())
        return true;
      if (a.x() > b.x())
        return false;
      return a.y() < b.y();
    };
    for (size_t level = 0; level < 3; level++)
      for (size_t i = 0; i < 16; i += 16 / (1 << level))
        for (size_t j = 0; j < 16; j += 16 / (1 << level)) {
          BBox2d b(Point2d(i, j) / 16.,
                   Point2d(i + (16 / (1 << level)), j + (16 / (1 << level))) /
                       16.);
          std::cerr << "REGION " << b << " level " << level << std::endl;
          std::vector<Point2d> searchResult;
          zgrid.searchZ(b, level, [&](size_t id) {
            searchResult.emplace_back(zgrid[id]);
          });
          // brute force check
          for (auto r : searchResult)
            EXPECT_TRUE(b.contains(r));
          std::vector<Point2d> trueSearch;
          for (auto r : points)
            if (b.contains(r))
              trueSearch.emplace_back(r);
          std::sort(trueSearch.begin(), trueSearch.end(), comp);
          std::sort(searchResult.begin(), searchResult.end(), comp);
          // for (size_t r = 0; r < searchResult.size(); r++)
          //   std::cerr << searchResult[r] << " ";
          // std::cerr << std::endl;
          // for (size_t r = 0; r < trueSearch.size(); r++)
          //   std::cerr << trueSearch[r] << " ";
          // std::cerr << std::endl;
          // std::cerr << std::endl;
          // std::cerr << std::endl;
          EXPECT_EQ(trueSearch.size(), searchResult.size());
          for (size_t r = 0; r < searchResult.size(); r++)
            EXPECT_EQ(trueSearch[r], searchResult[r]);
        }
  }
  return;
  {
    BoxSampler rng;
    PointZGrid2d zgrid(4u);
    BBox2d region(Point2d(0., 0.), Point2d(1, 1));
    zgrid.setDomainRegion(region);
    std::vector<Point2d> points;
    for (size_t j = 0; j < 400; j++) {
      auto p = rng.sample(region);
      zgrid.add(p);
      points.emplace_back(p);
    }
    EXPECT_EQ(zgrid.size(), 400u);
    for (size_t i = 0; i < 400; i++)
      EXPECT_EQ(points[i], zgrid[i]);
    std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
      if (a.x() < b.x())
        return true;
      if (a.x() > b.x())
        return false;
      return a.y() < b.y();
    };
    for (size_t level = 0; level < 3; level++)
      for (size_t i = 0; i < 16; i += 16 / (1 << level))
        for (size_t j = 0; j < 16; j += 16 / (1 << level)) {
          BBox2d b(Point2d(i, j) / 16.,
                   Point2d(i + (16 / (1 << level)), j + (16 / (1 << level))) /
                       16.);
          std::cerr << "REGION " << b << std::endl;
          std::vector<Point2d> searchResult;
          zgrid.searchZ(b, level + 1, [&](size_t id) {
            searchResult.emplace_back(zgrid[id]);
          });
          // brute force check
          for (auto r : searchResult)
            EXPECT_TRUE(b.contains(r));
          std::vector<Point2d> trueSearch;
          for (auto r : points)
            if (b.contains(r))
              trueSearch.emplace_back(r);
          std::sort(trueSearch.begin(), trueSearch.end(), comp);
          std::sort(searchResult.begin(), searchResult.end(), comp);
          for (size_t r = 0; r < searchResult.size(); r++)
            std::cerr << searchResult[r] << " ";
          std::cerr << std::endl;
          for (size_t r = 0; r < trueSearch.size(); r++)
            std::cerr << trueSearch[r] << " ";
          std::cerr << std::endl;
          std::cerr << std::endl;
          std::cerr << std::endl;
          EXPECT_EQ(trueSearch.size(), searchResult.size());
          for (size_t r = 0; r < searchResult.size(); r++)
            EXPECT_EQ(trueSearch[r], searchResult[r]);
        }
  }
}

TEST(PointZGrid, Constructors) {
  {
    PointZGrid2d grid(8);
    grid.add(Point2d(1, 1));
    // grid.search(BBox2d(Point2d(0), Point2d(20)), [](size_t id) {
    //   UNUSED_VARIABLE(id);
    // });
  }
}

TEST(PointZGrid, Add) {
  {
    PointZGrid2d grid(16);
    grid.setDomainRegion(BBox2d::squareBox());
    for (size_t i = 0; i < 40; i++)
      for (size_t j = 0; j < 40; j++)
        grid.add(Point2d(0.007 * i, 0.007 * j));
    EXPECT_EQ(grid.size(), 40u * 40u);
    std::set<size_t> ids;
    grid.iteratePoints([&](size_t i, Point2d p) {
      UNUSED_VARIABLE(p);
      ids.insert(i);
    });
    EXPECT_EQ(ids.size(), grid.size());
  }
}

TEST(PointZGrid, Remove) {
  {
    PointZGrid2d grid(16);
    grid.setDomainRegion(BBox2d::squareBox());
    for (size_t i = 0; i < 5; i++)
      for (size_t j = 0; j < 5; j++)
        grid.add(Point2d(0.2 * i, 0.2 * j));
    grid.size();
    grid.iteratePoints([&](size_t id, Point2d p) {
      UNUSED_VARIABLE(p);
      std::cerr << id << " ";
    });
    std::cerr << std::endl;
    for (size_t i = 0; i < 3; i++)
      grid.remove(i);
    std::set<size_t> ids;
    grid.iteratePoints([&](size_t id, Point2d p) {
      ids.insert(id);
      std::cerr << id << " ";
      UNUSED_VARIABLE(p);
    });
    std::cerr << std::endl;
    EXPECT_EQ(ids.size(), 5u * 5u - 3u);
    EXPECT_EQ(grid.size(), 5u * 5u - 3u);
    ids.clear();
    grid.search(BBox2d::squareBox(), [&](size_t id) {
      std::cerr << id << " ";
      ids.insert(id);
    });
    std::cerr << std::endl;
    EXPECT_EQ(ids.size(), 5u * 5u - 3u);
  }
}

TEST(PointZGrid, knn) {
  {
    BoxSampler rng;
    PointZGrid2d zgrid(4u);
    BBox2d region(Point2d(1., 1.), Point2d(3, 3));
    std::vector<Point2d> points;
    for (size_t j = 0; j < 400; j++) {
      auto p = rng.sample(region);
      zgrid.add(p);
      points.emplace_back(p);
    }
    EXPECT_EQ(zgrid.size(), 400u);
    for (size_t i = 0; i < 400; i++)
      EXPECT_EQ(points[i], zgrid[i]);
    for (size_t i = 0; i < 10000; i++) {
      Point2d center = rng.sample(BBox2d(Point2d(0, 0), Point2d(4, 1)));
      std::function<bool(Point2d, Point2d)> comp2 = [center](Point2d a,
                                                             Point2d b) {
        return distance2(a, center) < distance2(b, center);
      };
      std::vector<Point2d> searchResult;
      zgrid.iterateClosestPoints(center, 13,
                                 [&](size_t id, Point2d p) {
                                   UNUSED_VARIABLE(p);
                                   searchResult.emplace_back(zgrid[id]);
                                 },
                                 nullptr);
      EXPECT_EQ(searchResult.size(), 13u);
      // brute force check
      std::sort(points.begin(), points.end(), comp2);
      std::sort(searchResult.begin(), searchResult.end(), comp2);
      for (size_t r = 0; r < 13; r++)
        EXPECT_DOUBLE_EQ(distance2(center, points[r]),
                         distance2(center, searchResult[r]));
    }
  }
  {
    BoxSampler rng;
    PointZGrid2d zgrid(4u);
    BBox2d region(Point2d(0., 0.), Point2d(4, 4));
    std::vector<Point2d> points;
    for (size_t j = 0; j < 400; j++) {
      auto p = rng.sample(region);
      zgrid.add(p);
      points.emplace_back(p);
    }
    EXPECT_EQ(zgrid.size(), 400u);
    for (size_t i = 0; i < 400; i++)
      EXPECT_EQ(points[i], zgrid[i]);
    for (size_t i = 0; i < 10000; i++) {
      Point2d center = rng.sample(region);
      std::function<bool(Point2d, Point2d)> comp2 = [center](Point2d a,
                                                             Point2d b) {
        return distance2(a, center) < distance2(b, center);
      };
      std::vector<Point2d> searchResult;
      zgrid.iterateClosestPoints(center, 13,
                                 [&](size_t id, Point2d p) {
                                   UNUSED_VARIABLE(p);
                                   searchResult.emplace_back(zgrid[id]);
                                 },
                                 nullptr);
      EXPECT_EQ(searchResult.size(), 13u);
      // brute force check
      std::sort(points.begin(), points.end(), comp2);
      std::sort(searchResult.begin(), searchResult.end(), comp2);
      for (size_t r = 0; r < 13; r++)
        EXPECT_DOUBLE_EQ(distance2(center, points[r]),
                         distance2(center, searchResult[r]));
    }
  }
}

TEST(PointZGrid, Search_And_Operators) {
  {
    PointZGrid2d zgrid(64);
    zgrid.setDomainRegion(BBox2d::squareBox(4));
    zgrid.add(Point2d(0.5, 0.5));
    zgrid.add(Point2d(3.5, 3.5));
    zgrid.add(Point2d(1.5, 1.5));
    {
      std::vector<size_t> r;
      zgrid.search(BBox2d::squareBox(4),
                   [&](size_t id) { r.emplace_back(id); });
      EXPECT_EQ(r.size(), 3u);
    }
    {
      std::vector<size_t> r;
      zgrid.search(BBox2d(Point2d(), Point2d(1)),
                   [&](size_t id) { r.emplace_back(id); });
      EXPECT_EQ(r.size(), 1u);
    }
  }
  {
    BoxSampler rng;
    PointZGrid2d zgrid(4u);
    BBox2d region(Point2d(0., 0.), Point2d(4, 4));
    std::vector<Point2d> points;
    for (size_t j = 0; j < 400; j++) {
      auto p = rng.sample(region);
      zgrid.add(p);
      points.emplace_back(p);
    }
    EXPECT_EQ(zgrid.size(), 400u);
    for (size_t i = 0; i < 400; i++)
      EXPECT_EQ(points[i], zgrid[i]);
    std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
      if (a.x() < b.x())
        return true;
      if (a.x() > b.x())
        return false;
      return a.y() < b.y();
    };
    for (size_t i = 0; i < 10000; i++) {
      BBox2d b(rng.sample(region), rng.sample(region));
      std::vector<Point2d> searchResult;
      zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
      // brute force check
      for (auto r : searchResult)
        EXPECT_TRUE(b.contains(r));
      std::vector<Point2d> trueSearch;
      for (auto r : points)
        if (b.contains(r))
          trueSearch.emplace_back(r);
      EXPECT_EQ(trueSearch.size(), searchResult.size());
      std::sort(trueSearch.begin(), trueSearch.end(), comp);
      std::sort(searchResult.begin(), searchResult.end(), comp);
      for (size_t r = 0; r < searchResult.size(); r++)
        EXPECT_EQ(trueSearch[r], searchResult[r]);
    }
    for (size_t i = 0; i < 10000; i++) {
      BBox2d b(rng.sample(region), rng.sample(region));
      std::vector<Point2d> searchResult;
      zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
      // brute force check
      for (auto r : searchResult)
        EXPECT_TRUE(b.contains(r));
      std::vector<Point2d> trueSearch;
      for (auto r : points)
        if (b.contains(r))
          trueSearch.emplace_back(r);
      ASSERT_FATAL(trueSearch.size() == searchResult.size());
      std::sort(trueSearch.begin(), trueSearch.end(), comp);
      std::sort(searchResult.begin(), searchResult.end(), comp);
      for (size_t r = 0; r < searchResult.size(); r++)
        EXPECT_EQ(trueSearch[r], searchResult[r]);
    }
  }
}

TEST(PointZGrid2, Constructors) {
  {
    PointZGrid2 zg(16u, 1.);
    EXPECT_EQ(zg.nbits(), 4u);
    EXPECT_EQ(zg.maxDepth(), 4u);
  }
  {
    PointZGrid2 zg(64u, 1.);
    EXPECT_EQ(zg.nbits(), 6u);
    EXPECT_EQ(zg.maxDepth(), 6u);
  }
  {
    PointZGrid2 zg(2u, 8u, Transform2D());
    EXPECT_EQ(zg.nbits(), 3u);
    EXPECT_EQ(zg.maxDepth(), 3u);
  }
}

TEST(PointZGrid2, PointIterator) {
  {
    PointZGrid2 zgrid(16u, 32u, Transform2D());
    for (int i = 0; i < 16; i++)
      for (int j = 0; j < 16; j++)
        zgrid.add(Point2d({1. * i, 1. * j}));
    PointZGrid2::PointIterator it(zgrid);
    int i = 0, j = 0;
    size_t k = 0;
    while (it.next()) {
      EXPECT_EQ(Point2d({1. * i, 1. * j}), it.getWorldPosition());
      EXPECT_EQ(k++, it.getId());
      j++;
      if (j >= 16) {
        j = 0;
        i++;
      }
      EXPECT_EQ(
          it.pointElement()->zcode,
          mortonCode(it.getWorldPosition().x(), it.getWorldPosition().y()));
      ++it;
    }
  }
}

TEST(PointZGrid2, Update) {
  PointZGrid2 grid(8);
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
      grid.add(Point2d(i, j));
  grid.update();
  std::vector<Point2d> order;
  int ddx = 0, ddy = 0;
  int dx = 0 + ddx, dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 0 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);

  ddx = 4;
  ddy = 0;
  dx = 0 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 0 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);

  ddx = 0;
  ddy = 4;
  dx = 0 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 0 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);

  ddx = 4;
  ddy = 4;
  dx = 0 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 0 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 0 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);
  dx = 2 + ddx;
  dy = 2 + ddy;
  order.emplace_back(0 + dx, 0 + dy);
  order.emplace_back(1 + dx, 0 + dy);
  order.emplace_back(0 + dx, 1 + dy);
  order.emplace_back(1 + dx, 1 + dy);

  EXPECT_EQ(order.size(), 8u * 8u);
  int i = 0;
  for (PointZGrid2::PointIterator it(grid); it.next(); ++it) {
    auto position = it.getWorldPosition();
    EXPECT_EQ(order[i], position);
    i++;
  }
}

TEST(PointZGrid2, SearchTree) {
  {
    PointZGrid2 zgrid(4u);
    PointZGrid2::SearchTree tree(zgrid);
    tree.traverse([&](QuadTree::Node &node) -> bool {
      EXPECT_EQ(
          mortonCode(node.region().lower().x(), node.region().lower().y()),
          tree.ids[node.id]);
      return true;
    });
    std::vector<Point2d> points;
    for (size_t i = 1; i < 10; i++)
      for (size_t j = 1; j < 10; j++) {
        zgrid.add(Point2d(j * .4, i * .4));
        points.emplace_back(j * .4, i * .4);
      }
    EXPECT_EQ(zgrid.size(), 9u * 9u);
    zgrid.update();
    size_t k = 0;
    zgrid.searchTree().traverse([&](const QuadTree::Node &node) -> bool {
      EXPECT_EQ(
          mortonCode(node.region().lower().x(), node.region().lower().y()),
          tree.ids[node.id]);
      k++;
      return true;
    });
    EXPECT_EQ(k, 21u);
    for (size_t i = 0; i < 4; i++)
      for (size_t j = 0; j < 4; j++) {
        BBox2D b(Point2d(1. * j, 1. * i), Point2d(1. * (j + 1), 1. * (i + 1)));
        k = 0;
        std::vector<Point2d> searchResult;
        zgrid.searchTree().iteratePoints(b, [&](PointZGrid2::PointElement *p) {
          searchResult.emplace_back(zgrid[p->id]);
        });
        // brute force check
        for (auto r : searchResult)
          EXPECT_EQ(b.inside(r), true);
        std::vector<Point2d> trueSearch;
        for (auto r : points) {
          if (b.inside(r))
            trueSearch.emplace_back(r);
        }
        EXPECT_EQ(trueSearch.size(), searchResult.size());
        std::sort(trueSearch.begin(), trueSearch.end(),
                  [](Point2d a, Point2d b) {
                    if (a.x() < b.x())
                      return true;
                    if (a.x() > b.x())
                      return false;
                    return a.y() < b.y();
                  });
        std::sort(searchResult.begin(), searchResult.end(),
                  [](Point2d a, Point2d b) {
                    if (a.x() < b.x())
                      return true;
                    if (a.x() > b.x())
                      return false;
                    return a.y() < b.y();
                  });
        for (size_t r = 0; r < searchResult.size(); r++)
          EXPECT_EQ(trueSearch[r], searchResult[r]);
      }
    zgrid.add(Point2d(0., 0.));
    zgrid.add(Point2d(4., 4.));
    zgrid.add(Point2d(0., 3.));
    zgrid.add(Point2d(3., 0.));
    zgrid.update();
    k = 0;
    size_t expectedIds[6] = {19, 20, 21, 22, 23, 24};
    std::vector<size_t> res;
    zgrid.searchTree().iteratePoints(
        BBox2D(Point2d(0.6, 1.), Point2d(3., 1.4)),
        [&](PointZGrid2::PointElement *p) { res.emplace_back(p->id); });
    EXPECT_EQ(res.size(), 6u);
    std::sort(res.begin(), res.end());
    for (k = 0; k < res.size(); k++)
      EXPECT_EQ(res[k], expectedIds[k]);
  }
  return;
  {
    BoxSampler rng;
    Transform2D t = Transform2D::scale(0.1, 0.2);
    t.applyTranslate(Vector2d(1., 0.));
    PointZGrid2 zgrid(4u, 4u, t);
    BBox2D region(t(Point2d(0., 0.)), t(Point2d(2, 4)));
    std::vector<Point2d> points;
    for (size_t j = 1; j < 40; j++)
      for (size_t i = 1; i < 20; i++) {
        zgrid.add(Point2d(i * .01 + 0.1, j * .02));
        points.emplace_back(i * .01 + 0.1, j * .02);
      }
    zgrid.update();
    zgrid.searchTree().traverse([&](QuadTree::Node &n) -> bool {
      EXPECT_EQ(mortonCode(n.region().lower().x(), n.region().lower().y()),
                zgrid.searchTree().ids[n.id]);
      return true;
    });
    for (size_t i = 0; i < 10000; i++) {
      BBox2D b(rng.sample(region), rng.sample(region));
      std::vector<Point2d> searchResult;
      zgrid.searchTree().iteratePoints(b, [&](PointZGrid2::PointElement *p) {
        searchResult.emplace_back(zgrid[p->id]);
      });
      // brute force check
      FILE *fp = fopen("data", "w+");
      // brute force check
      for (auto r : searchResult) {
        EXPECT_EQ(b.inside(r), true);
        fprintf(fp, "%f %f\n", r.x(), r.y());
      }
      fclose(fp);
      fp = fopen("data2", "w+");
      std::vector<Point2d> trueSearch;
      for (auto r : points) {
        if (b.inside(r)) {
          trueSearch.emplace_back(r);
        }
        fprintf(fp, "%f %f\n", r.x(), r.y());
      }
      fclose(fp);

      EXPECT_EQ(trueSearch.size(), searchResult.size());
      ASSERT_FATAL(trueSearch.size() == searchResult.size());
      std::sort(trueSearch.begin(), trueSearch.end(), [](Point2d a, Point2d b) {
        if (a.x() < b.x())
          return true;
        if (a.x() > b.x())
          return false;
        return a.y() < b.y();
      });
      std::sort(searchResult.begin(), searchResult.end(),
                [](Point2d a, Point2d b) {
                  if (a.x() < b.x())
                    return true;
                  if (a.x() > b.x())
                    return false;
                  return a.y() < b.y();
                });
      for (size_t r = 0; r < searchResult.size(); r++)
        EXPECT_EQ(trueSearch[r], searchResult[r]);
    }
  }
  {
    Transform2D t = Transform2D::translate(Vector2d(0.2, 0.2));
    t.applyScale(0.1, 0.1);
    PointZGrid2 zgrid(256, 256, Transform2D());
    BBox2D region = Transform2D::scale(256, 256)(BBox2D::make_unit_bbox());
    BoxSampler rng;
    std::vector<Point2d> points;
    for (size_t i = 0; i < 100000; i++) {
      points.emplace_back(rng.sample(region));
      zgrid.add(points[i]);
    }
    EXPECT_EQ(zgrid.size(), 100000u);
    zgrid.update();
    for (size_t i = 0; i < 1000; i++) {
      BBox2D b(rng.sample(region), rng.sample(region));
      std::cout << "for bbox " << i << " " << b;
      std::vector<Point2d> searchResult;
      zgrid.searchTree().iteratePoints(b, [&](PointZGrid2::PointElement *p) {
        searchResult.emplace_back(zgrid[p->id]);
      });
      FILE *fp = fopen("data", "w+");
      // brute force check
      for (auto r : searchResult) {
        EXPECT_EQ(b.inside(r), true);
        fprintf(fp, "%f %f\n", r.x(), r.y());
      }
      fclose(fp);
      fp = fopen("data2", "w+");
      std::vector<Point2d> trueSearch;
      for (auto r : points) {
        if (b.inside(r)) {
          trueSearch.emplace_back(r);
        }
        fprintf(fp, "%f %f\n", r.x(), r.y());
      }
      fclose(fp);
      EXPECT_EQ(trueSearch.size(), searchResult.size());
      ASSERT_FATAL(trueSearch.size() == searchResult.size());
      std::sort(trueSearch.begin(), trueSearch.end(), [](Point2d a, Point2d b) {
        if (a.x() < b.x())
          return true;
        if (a.x() > b.x())
          return false;
        return a.y() < b.y();
      });
      std::sort(searchResult.begin(), searchResult.end(),
                [](Point2d a, Point2d b) {
                  if (a.x() < b.x())
                    return true;
                  if (a.x() > b.x())
                    return false;
                  return a.y() < b.y();
                });
      for (size_t r = 0; r < searchResult.size(); r++)
        EXPECT_EQ(trueSearch[r], searchResult[r]);
    }
  }
}

TEST(PointZGrid2, Search_And_Operators) {
  BoxSampler rng;
  Transform2D t = Transform2D::scale(0.1, 0.2);
  t.applyTranslate(Vector2d(1., 0.));
  PointZGrid2 zgrid(4u, 4u, t);
  BBox2D region(t(Point2d(0., 0.)), t(Point2d(4, 4)));
  std::vector<Point2d> points;
  for (size_t j = 0; j < 400; j++) {
    auto p = rng.sample(region);
    zgrid.add(p);
    points.emplace_back(p);
  }
  EXPECT_EQ(zgrid.size(), 400u);
  for (size_t i = 0; i < 400; i++)
    EXPECT_EQ(points[i], zgrid[i]);
  zgrid.update();
  std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
    if (a.x() < b.x())
      return true;
    if (a.x() > b.x())
      return false;
    return a.y() < b.y();
  };
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;
    zgrid.searchTree().iteratePoints(b, [&](PointZGrid2::PointElement *p) {
      searchResult.emplace_back(zgrid[p->id]);
    });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    EXPECT_EQ(trueSearch.size(), searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
  // for (size_t i = 0; i < 400; i++)
  //  zgrid[i] = rng.sample(region);
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;
    zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    ASSERT_FATAL(trueSearch.size() == searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
}

TEST(PointZGrid2, Search_And_Operators_No_Tree) {
  BoxSampler rng;
  Transform2D t = Transform2D::scale(0.1, 0.2);
  t.applyTranslate(Vector2d(1., 0.));
  PointZGrid2 zgrid(4u, 4u, t);
  BBox2D region(t(Point2d(0., 0.)), t(Point2d(4, 4)));
  std::vector<Point2d> points;
  for (size_t j = 0; j < 400; j++) {
    auto p = rng.sample(region);
    zgrid.add(p);
    points.emplace_back(p);
  }
  EXPECT_EQ(zgrid.size(), 400u);
  for (size_t i = 0; i < 400; i++)
    EXPECT_EQ(points[i], zgrid[i]);
  zgrid.update();
  std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
    if (a.x() < b.x())
      return true;
    if (a.x() > b.x())
      return false;
    return a.y() < b.y();
  };
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;

    zgrid.searchTree().iteratePoints(b, [&](PointZGrid2::PointElement *p) {
      searchResult.emplace_back(zgrid[p->id]);
    });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    EXPECT_EQ(trueSearch.size(), searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
  // for (size_t i = 0; i < 400; i++)
  //  zgrid[i] = rng.sample(region);
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;
    zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    ASSERT_FATAL(trueSearch.size() == searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
}

TEST(PointZGrid2, setActive) {
  PointZGrid2 zgrid(16);
  BoxSampler rng;
  std::vector<Point2d> points;
  for (size_t j = 0; j < 400; j++) {
    auto p = rng.sample();
    zgrid.add(p);
    points.emplace_back(p);
  }
  EXPECT_EQ(zgrid.size(), 400u);
  for (size_t i = 0; i < 400; i++)
    zgrid.remove(i);
  EXPECT_EQ(zgrid.size(), 0u);
  for (size_t i = 0; i < 200; i++)
    zgrid.add(rng.sample());
  EXPECT_EQ(zgrid.size(), 200u);
}

TEST(PointZGrid2, iteratePoints) {
  PointZGrid2 zgrid(16);
  BBox2D region = Transform2D::scale(16, 16)(BBox2D::make_unit_bbox());
  BoxSampler rng;
  for (size_t j = 0; j < 200; j++)
    zgrid.add(rng.sample(region));
  zgrid.update();
  std::vector<size_t> indices;
  zgrid.iteratePoints([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(p);
    indices.emplace_back(i);
  });
  EXPECT_EQ(indices.size(), 200u);
  for (size_t i = 0; i < 100; i++)
    zgrid.remove(i);
  zgrid.update();
  size_t j = 0;
  zgrid.iteratePoints([&](unsigned int i, Point2d p) {
    static unsigned int k = 0;
    UNUSED_VARIABLE(i);
    UNUSED_VARIABLE(p);
    while (indices[k] < 100)
      k++;
    EXPECT_EQ(indices[k], i);
    k++;
    j++;
  });
  EXPECT_EQ(j, 100u);
}

TEST(PointZGrid2, setPosition) {
  BoxSampler rng;
  Transform2D t = Transform2D::scale(0.1, 0.2);
  t.applyTranslate(Vector2d(1., 0.));
  PointZGrid2 zgrid(4u, 4u, t);
  BBox2D region(t(Point2d(0., 0.)), t(Point2d(4, 4)));
  std::vector<Point2d> points;
  for (size_t j = 0; j < 400; j++) {
    auto p = rng.sample(region);
    zgrid.add(p);
    points.emplace_back(p);
  }
  zgrid.update();
  std::function<bool(Point2d, Point2d)> comp = [](Point2d a, Point2d b) {
    if (a.x() < b.x())
      return true;
    if (a.x() > b.x())
      return false;
    return a.y() < b.y();
  };
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;
    zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    ASSERT_FATAL(trueSearch.size() == searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
  // now change positions
  for (size_t i = 0; i < 400; i++) {
    auto p = rng.sample(region);
    zgrid.setPosition(i, p);
    points[i] = p;
  }
  for (size_t i = 0; i < 10000; i++) {
    BBox2D b(rng.sample(region), rng.sample(region));
    std::vector<Point2d> searchResult;
    zgrid.search(b, [&](size_t id) { searchResult.emplace_back(zgrid[id]); });
    // brute force check
    for (auto r : searchResult)
      EXPECT_EQ(b.inside(r), true);
    std::vector<Point2d> trueSearch;
    for (auto r : points)
      if (b.inside(r))
        trueSearch.emplace_back(r);
    ASSERT_FATAL(trueSearch.size() == searchResult.size());
    std::sort(trueSearch.begin(), trueSearch.end(), comp);
    std::sort(searchResult.begin(), searchResult.end(), comp);
    for (size_t r = 0; r < searchResult.size(); r++)
      EXPECT_EQ(trueSearch[r], searchResult[r]);
  }
}

TEST(PointZGrid2, extremeCases) {
  PointZGrid2 zgrid(2u, 2u, BBox2D(Point2d(), Point2d(2, 2)));
  zgrid.add(Point2d(0., 0.));
  zgrid.add(Point2d(2., 2.));
  zgrid.add(Point2d(1.99999, 1.99999));
  zgrid.add(Point2d(1.5, 0.));
  zgrid.add(Point2d(1.5, 2.));
  zgrid.add(Point2d(0., 1.5));
  zgrid.add(Point2d(2., 1.5));
  zgrid.add(Point2d(1., 1.));
  PointZGrid2::PointIterator it(zgrid);
  while (it.next()) {
    std::cout << it.getWorldPosition();
    std::cout << it.getId() << std::endl;
    printBits(it.pointElement()->zcode);
    std::cout << std::endl;
    ++it;
  }
  zgrid.update();
  auto tree = zgrid.searchTree();
  std::cout << "tree height " << tree.height() << std::endl;
  std::cout << tree.leafCount() << std::endl;
  std::cout << "traverse\n";
  tree.traverse([&](const QuadTree::Node &n) -> bool {
    std::cout << n.id << " ---> ";
    printBits(tree.ids[n.id]);
    std::cout << std::endl;
    return true;
  });

  zgrid.search(BBox2D(Point2d(1., 1.), Point2d(2., 2.)),
               [](unsigned id) { std::cout << id << std::endl; });
  std::cout << "cel\n";
  zgrid.search(BBox2D(Point2d(0., 0.), Point2d(1., 1.)),
               [](unsigned id) { std::cout << id << std::endl; });
  std::cout << "cel\n";
  zgrid.search(BBox2D(Point2d(0., 1.), Point2d(1., 2.)),
               [](unsigned id) { std::cout << id << std::endl; });
  std::cout << "cel\n";
  zgrid.search(BBox2D(Point2d(1., 1.), Point2d(2., 2.)),
               [](unsigned id) { std::cout << id << std::endl; });
  std::cout << "cel\n";
  zgrid.search(BBox2D(Point2d(1., 0.), Point2d(2., 1.)),
               [](unsigned id) { std::cout << id << std::endl; });
}

TEST(PointZGrid2, knn) {
  {
    PointZGrid2 zgrid(4u, 4u, BBox2D(Point2d(), Point2d(4, 4)));
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++)
        zgrid.add(Point2d(i, j));
    zgrid.update();
    for (PointZGrid2::PointIterator it(zgrid); it.next(); ++it)
      std::cerr << it.pointElement()->zcode << " ";
    std::cerr << std::endl;
    zgrid.iterateClosestPoints(Point2d(0, 0), 2u,
                               [](size_t i, Point2d p) {
                                 std::cerr << i << " ------------------ " << p
                                           << std::endl;
                               },
                               nullptr);
  }
  {
    PointZGrid2 zgrid(4u, 4u, BBox2D(Point2d(), Point2d(4, 4)));
    for (int j = 2; j < 4; j++)
      for (int i = 2; i < 4; i++)
        zgrid.add(Point2d(i, j));
    zgrid.update();
    for (PointZGrid2::PointIterator it(zgrid); it.next(); ++it)
      std::cerr << it.pointElement()->zcode << " ";
    std::cerr << std::endl;
    zgrid.iterateClosestPoints(Point2d(0, 0), 2,
                               [](size_t i, Point2d p) {
                                 std::cerr << i << " ------------------ " << p
                                           << std::endl;
                               },
                               nullptr);
  }
  {
    PointZGrid2 zgrid(4u, 4u, BBox2D(Point2d(), Point2d(4, 4)));
    for (int i = 0; i < 4; i++) {
      zgrid.add(Point2d(i, 0));
      zgrid.add(Point2d(i, 3));
    }
    zgrid.update();
    for (PointZGrid2::PointIterator it(zgrid); it.next(); ++it)
      std::cerr << it.pointElement()->zcode << " ";
    std::cerr << std::endl;
    zgrid.iterateClosestPoints(Point2d(2, 1.5), 2,
                               [](size_t i, Point2d p) {
                                 std::cerr << i << " ------------------ " << p
                                           << std::endl;
                               },
                               nullptr);
  }
  {
    PointZGrid2 zgrid(8u, 8u, BBox2D(Point2d(), Point2d(8, 8)));
    for (int i = 4; i < 7; i++)
      zgrid.add(Point2d(i, 3));
    for (int i = 0; i < 3; i++)
      zgrid.add(Point2d(4, i));
    zgrid.update();
    for (PointZGrid2::PointIterator it(zgrid); it.next(); ++it)
      std::cerr << it.pointElement()->zcode << " ";
    std::cerr << std::endl;
    zgrid.iterateClosestPoints(Point2d(3, 3), 5,
                               [](size_t i, Point2d p) {
                                 std::cerr << i << " ------------------ " << p
                                           << std::endl;
                               },
                               nullptr);
  }
}
