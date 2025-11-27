#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(BBox, Constructors) {
  {
    BBox2d b(Point2d(-1.0, -1.0), Point2d(1.0, 1.0));
    EXPECT_EQ(b.lower(), Point2d(-1.0, -1.0));
    EXPECT_EQ(b.upper(), Point2d(1.0, 1.0));
    BBox2d c(Point2d(1.0, 1.0), Point2d(-1.0, -1.0));
    EXPECT_EQ(c.lower(), Point2d(-1.0, -1.0));
    EXPECT_EQ(c.upper(), Point2d(1.0, 1.0));
  }
  {
    BBox3d b(Point3d(-1.0), Point3d(1.0));
    EXPECT_EQ(b.lower(), Point3d(-1.0));
    EXPECT_EQ(b.upper(), Point3d(1.0));
    BBox3D c(Point3d(1.0), Point3d(-1.0));
    EXPECT_EQ(c.lower(), Point3d(-1.0));
    EXPECT_EQ(c.upper(), Point3d(1.0));
  }
}

TEST(BBox, Statics) {
  {
    auto b = BBox2d::squareBox();
    EXPECT_EQ(b.lower(), Point2d());
    EXPECT_EQ(b.upper(), Point2d(1));
  }
  {
    auto b = BBox2d::squareBox(10);
    EXPECT_EQ(b.lower(), Point2d(0, 0));
    EXPECT_EQ(b.upper(), Point2d(10, 10));
    EXPECT_EQ(b.center(), Point2d(5, 5));
  }
  {
    auto children = BBox2d::halfSplit(BBox2d::squareBox(4));
    EXPECT_EQ(children.size(), 4u);
    EXPECT_EQ(children[0].lower(), Point2d(0));
    EXPECT_EQ(children[0].upper(), Point2d(2));
    EXPECT_EQ(children[1].lower(), Point2d(2, 0));
    EXPECT_EQ(children[1].upper(), Point2d(4, 2));
    EXPECT_EQ(children[2].lower(), Point2d(0, 2));
    EXPECT_EQ(children[2].upper(), Point2d(2, 4));
    EXPECT_EQ(children[3].lower(), Point2d(2));
    EXPECT_EQ(children[3].upper(), Point2d(4));
  }
  {
    auto children = BBox3d::halfSplit(BBox3d::squareBox(4));
    EXPECT_EQ(children.size(), 8u);
    EXPECT_EQ(children[0].lower(), Point3d(0, 0, 0));
    EXPECT_EQ(children[0].upper(), Point3d(2, 2, 2));
    EXPECT_EQ(children[1].lower(), Point3d(2, 0, 0));
    EXPECT_EQ(children[1].upper(), Point3d(4, 2, 2));
    EXPECT_EQ(children[2].lower(), Point3d(0, 2, 0));
    EXPECT_EQ(children[2].upper(), Point3d(2, 4, 2));
    EXPECT_EQ(children[3].lower(), Point3d(2, 2, 0));
    EXPECT_EQ(children[3].upper(), Point3d(4, 4, 2));
    EXPECT_EQ(children[4].lower(), Point3d(0, 0, 2));
    EXPECT_EQ(children[4].upper(), Point3d(2, 2, 4));
    EXPECT_EQ(children[5].lower(), Point3d(2, 0, 2));
    EXPECT_EQ(children[5].upper(), Point3d(4, 2, 4));
    EXPECT_EQ(children[6].lower(), Point3d(0, 2, 2));
    EXPECT_EQ(children[6].upper(), Point3d(2, 4, 4));
    EXPECT_EQ(children[7].lower(), Point3d(2, 2, 2));
    EXPECT_EQ(children[7].upper(), Point3d(4, 4, 4));
  }
}

TEST(BBox2D, Constructors) {
  BBox2D b(Point2d(-1.0, -1.0), Point2d(1.0, 1.0));
  EXPECT_EQ(b.lower() == Point2d(-1.0, -1.0), true);
  EXPECT_EQ(b.upper() == Point2d(1.0, 1.0), true);
  BBox2D c(Point2d(1.0, 1.0), Point2d(-1.0, -1.0));
  EXPECT_EQ(c.lower() == Point2d(-1.0, -1.0), true);
  EXPECT_EQ(c.upper() == Point2d(1.0, 1.0), true);
  EXPECT_EQ(c.lower() == c[0], true);
  EXPECT_EQ(c.upper() == c[1], true);
  BBox2D r(Point2d(10., 5.), 5.);
  EXPECT_EQ(r.lower(), Point2d(5., 0.));
  EXPECT_EQ(r.upper(), Point2d(15., 10.));
}

TEST(BBox, contains) {
  {
    auto b = BBox2d::squareBox(1);
    EXPECT_TRUE(b.contains(Point2d(0, 0)));
    EXPECT_TRUE(b.contains(b));
  }

}
TEST(BBox, Intersect) {
  {
    BBox2d b = BBox2d::squareBox();
    Vector2d n;
    Point2d p;
    {
      auto i = b.intersect(Point2d(0.1, 0.1), Vector2d(-1., 0.), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0., 0.1));
      EXPECT_EQ(n, Vector2d(1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.1, 0.2),
                      Vector2d(-1., -1.).normalized(),
                      p,
                      &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0., 0.1));
      EXPECT_EQ(n, Vector2d(1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.1, 0.2), Vector2d(1., -1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0.3, 0.));
      EXPECT_EQ(n, Vector2d(0., 1.));
    }
    {
      auto i =
          b.intersect(Point2d(0.1, 0.2), Vector2d(1., 1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0.9, 1.));
      EXPECT_EQ(n, Vector2d(0., -1.));
    }
    {
      auto i =
          b.intersect(Point2d(0.2, 0.1), Vector2d(1., 1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(1., 0.9));
      EXPECT_EQ(n, Vector2d(-1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.2, 0.1), Vector2d(1., 0.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(1., 0.1));
      EXPECT_EQ(n, Vector2d(-1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.2, 0.1), Vector2d(-1., 0.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0., 0.1));
      EXPECT_EQ(n, Vector2d(1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.2, 0.1), Vector2d(1., 2.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0.65, 1));
      EXPECT_EQ(n, Vector2d(0., -1.));
    }
    {
      auto i =
          b.intersect(Point2d(0.5, 0.5), Vector2d(-2., 1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(0., 0.75));
      EXPECT_EQ(n, Vector2d(1., 0.));
    }
    {
      auto i =
          b.intersect(Point2d(0.5, 0.5),
                      Vector2d(-1., -2.).normalized(),
                      p,
                      &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(.25, 0.));
      EXPECT_EQ(n, Vector2d(0., 1.));
    }
    {
      auto i =
          b.intersect(Point2d(0.5, 0.5),
                      Vector2d(-1., -1.).normalized(),
                      p,
                      &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(.0, 0.));
      EXPECT_EQ(n, Vector2d(1., 1.).normalized());
    }
    BBox2d bb(Point2d(-3.723, .2434235), Point2d(45.3425, -224.1399314));
    {
      auto i =
          bb.intersect(Point2d(10., -20.),
                       Vector2d(-1., 0.).normalized(),
                       p,
                       &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(-3.723, -20.));
      EXPECT_EQ(n, Vector2d(1., 0.).normalized());
    }
    {
      auto i =
          bb.intersect(Point2d(10., -20.),
                       Vector2d(0., 1.).normalized(),
                       p,
                       &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(10., .2434235));
      EXPECT_EQ(n, Vector2d(0., -1.).normalized());
    }
    {
      auto i = bb.intersect(Point2d(-3.723, -224.1399314) + Vector2d(1., 1.),
                            Vector2d(-1., -1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(-3.723, -224.1399314));
      EXPECT_EQ(n, Vector2d(1., 1.).normalized());
    }
    {
      auto i =
          b.intersect(Point2d(0.9, 0.5), Vector2d(1., 1.).normalized(), p, &n);
      EXPECT_EQ(i, true);
      EXPECT_EQ(p, Point2d(1., 0.6));
      EXPECT_EQ(n, Vector2d(-1., 0.).normalized());
    }
  }
}

TEST(BBox2D, Inside) {
  BBox2D b(Point2d(-1.0, -1.0), Point2d(1.0, 1.0));
  EXPECT_EQ(b.inside(Point2d(0.0, 0.0)), true);
  EXPECT_EQ(b.inside(Point2d(-1.0, 0.0)), true);
  EXPECT_EQ(b.inside(Point2d(1.0, 1.0)), true);
  EXPECT_EQ(b.inside(Point2d(1.0, 0.0)), true);
  EXPECT_EQ(b.inside(Point2d(1.000001, 0.0)), false);
}

TEST(BBox2D, Dimensions) {
  BBox2D b(Point2d(-1.0, -10.0), Point2d(1.0, 10.0));
  ASSERT_NEAR(b.size(0), 2.0, 1e-8);
  ASSERT_NEAR(b.size(1), 20.0, 1e-8);
  EXPECT_EQ(b.extends() == Vector2d(2.0, 20.0), true);
  EXPECT_EQ(b.center() == Point2d(0.0, 0.0), true);
  EXPECT_EQ(b.maxExtent(), 1);
}

TEST(BBox2D, Intersect) {
  BBox2D b = BBox2D::make_unit_bbox();
  Vector2d n;
  Point2d p;
  {
    auto i = b.intersect(Point2d(0.1, 0.1), Vector2d(-1., 0.), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0., 0.1));
    EXPECT_EQ(n, Vector2d(1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.1, 0.2), Vector2d(-1., -1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0., 0.1));
    EXPECT_EQ(n, Vector2d(1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.1, 0.2), Vector2d(1., -1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0.3, 0.));
    EXPECT_EQ(n, Vector2d(0., 1.));
  }
  {
    auto i =
        b.intersect(Point2d(0.1, 0.2), Vector2d(1., 1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0.9, 1.));
    EXPECT_EQ(n, Vector2d(0., -1.));
  }
  {
    auto i =
        b.intersect(Point2d(0.2, 0.1), Vector2d(1., 1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(1., 0.9));
    EXPECT_EQ(n, Vector2d(-1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.2, 0.1), Vector2d(1., 0.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(1., 0.1));
    EXPECT_EQ(n, Vector2d(-1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.2, 0.1), Vector2d(-1., 0.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0., 0.1));
    EXPECT_EQ(n, Vector2d(1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.2, 0.1), Vector2d(1., 2.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0.65, 1));
    EXPECT_EQ(n, Vector2d(0., -1.));
  }
  {
    auto i =
        b.intersect(Point2d(0.5, 0.5), Vector2d(-2., 1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(0., 0.75));
    EXPECT_EQ(n, Vector2d(1., 0.));
  }
  {
    auto i =
        b.intersect(Point2d(0.5, 0.5), Vector2d(-1., -2.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(.25, 0.));
    EXPECT_EQ(n, Vector2d(0., 1.));
  }
  {
    auto i =
        b.intersect(Point2d(0.5, 0.5), Vector2d(-1., -1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(.0, 0.));
    EXPECT_EQ(n, Vector2d(1., 1.).normalized());
  }
  BBox2D bb(Point2d(-3.723, .2434235), Point2d(45.3425, -224.1399314));
  {
    auto i =
        bb.intersect(Point2d(10., -20.), Vector2d(-1., 0.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(-3.723, -20.));
    EXPECT_EQ(n, Vector2d(1., 0.).normalized());
  }
  {
    auto i =
        bb.intersect(Point2d(10., -20.), Vector2d(0., 1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(10., .2434235));
    EXPECT_EQ(n, Vector2d(0., -1.).normalized());
  }
  {
    auto i = bb.intersect(Point2d(-3.723, -224.1399314) + Vector2d(1., 1.),
                          Vector2d(-1., -1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(-3.723, -224.1399314));
    EXPECT_EQ(n, Vector2d(1., 1.).normalized());
  }
  {
    auto i =
        b.intersect(Point2d(0.9, 0.5), Vector2d(1., 1.).normalized(), p, n);
    EXPECT_EQ(i, true);
    EXPECT_EQ(p, Point2d(1., 0.6));
    EXPECT_EQ(n, Vector2d(-1., 0.).normalized());
  }
}

TEST(BBox2D, PointCoordinates) {
  BBox2D bbox(Point2d(-5, -5), Point2d(5, 5));
  EXPECT_EQ(bbox.pointCoordinates(Point2d(0, 0)), Point2d(0.5, 0.5));
  EXPECT_EQ(bbox.pointCoordinates(Point2d(5, 5)), Point2d(1.0, 1.0));
  EXPECT_EQ(bbox.pointCoordinates(Point2d(-5, -5)), Point2d(0.0, 0.0));
}

TEST(BBox3D, Constructors) {
  BBox3D b(Point3d(-1.0), Point3d(1.0));
  EXPECT_EQ(b.lower() == Point3d(-1.0), true);
  EXPECT_EQ(b.upper() == Point3d(1.0), true);
  BBox3D c(Point3d(1.0), Point3d(-1.0));
  EXPECT_EQ(c.lower() == Point3d(-1.0), true);
  EXPECT_EQ(c.upper() == Point3d(1.0), true);
  EXPECT_EQ(c.lower() == c[0], true);
  EXPECT_EQ(c.upper() == c[1], true);
  BBox3D r(Point3d(10., 5., 7.), 5.);
  EXPECT_EQ(r.lower(), Point3d(5., 0., 2.));
  EXPECT_EQ(r.upper(), Point3d(15., 10., 12.));
}

TEST(BBox3D, Inside) {
  BBox3D b(Point3d(-1.0), Point3d(1.0));
  EXPECT_EQ(b.inside(Point3d(0.0)), true);
  EXPECT_EQ(b.inside(Point3d(-1., 0., 0.)), true);
  EXPECT_EQ(b.inside(Point3d(1., 1., 1.)), true);
  EXPECT_EQ(b.inside(Point3d(1.0, 0.0, 0.)), true);
  EXPECT_EQ(b.inside(Point3d(1.000001, 0.0, 0.)), false);
}

TEST(BBox3D, Dimensions) {
  BBox3D b(Point3d(-1.0, -10.0, 10.), Point3d(1.0, 10.0, 100.));
  ASSERT_NEAR(b.size(0), 2.0, 1e-8);
  ASSERT_NEAR(b.size(1), 20.0, 1e-8);
  ASSERT_NEAR(b.size(2), 90.0, 1e-8);
  EXPECT_EQ(b.extends() == Vector3d(2.0, 20.0, 90.0), true);
  EXPECT_EQ(b.center() == Point3d(0.0, 0.0, 55.0), true);
  EXPECT_EQ(b.maxExtent(), 2);
}

TEST(BBox3D, Intersect) {NOT_IMPLEMENTED(); }
