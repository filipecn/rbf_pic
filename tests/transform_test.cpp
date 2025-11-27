#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Transform, toUnitBox) {
  {
    auto t = Transform2::toUnitBox(
        BBox2d(Point2d(0.375, 0.25), Point2d(0.625, 0.375)));
    EXPECT_EQ(Point2d(0.,0.), t(Point2d(0.375, 0.25)));
    EXPECT_EQ(Point2d(1.,1.), t(Point2d(0.625, 0.375)));
  }
}

TEST(Transform, Constructors) {
  {
    Transform2 t;
    EXPECT_TRUE(t.matrix().isIdentity());
    EXPECT_TRUE(t.inverseMatrix().isIdentity());
    Matrix3d m;
    m.setIdentity();
    m(0, 0) = 2.0;
    m(1, 1) = 4.0;
    Transform2 t2(m, inverse(m));
    Point2d p = t2(Point2d(1.0, 1.0));
    EXPECT_EQ(p, Point2d(2.0, 4.0));
    EXPECT_EQ(Point2d(1.0, 1.0), t2.inverse()(p));
  }
  {
    Transform3 t;

    EXPECT_EQ(t.matrix().isIdentity(), true);
    EXPECT_EQ(t.inverseMatrix().isIdentity(), true);
    Matrix4d m;
    m.setIdentity();
    m(0, 0) = 3.0;
    m(1, 1) = 4.0;
    m(2, 2) = 5.0;
    Transform3D t2(m, inverse(m));
    Point3d p = t2(Point3d(1.0));
    EXPECT_EQ(p == Point3d(3.0, 4.0, 5.), true);
    EXPECT_EQ(Point3d(1.0) == t2.inverse()(p), true);
  }
}

TEST(Transform, Translate) {
  {
    Transform2 t = Transform2::translate(Vector2d(10., -20.));
    EXPECT_EQ(t(Point2d(-2.0, 1.0)), Point2d(8.0, -19.0));
    EXPECT_EQ(t.inverse()(Point2d(-20., -20.)), Point2d(-30., 0.));
    t = Transform2::translate(Vector2d(2.));
    EXPECT_EQ(t(Point2d(1.0, 1.0)), Point2d(3.0, 3.0));
    ASSERT_NEAR(t(Vector2d(1.0, 0.0)).length2(), 1.0, 1e-8);
  }
  {
    Transform3 t = Transform3::translate(Vector3d(10.0, -30.0, 5.));
    EXPECT_EQ(t(Point3d(-3.0, 1.0, 10.)), Point3d(7.0, -29.0, 15.));
    EXPECT_EQ(t.inverse()(Point3d(-30., -30., -50.)), Point3d(-40., 0., -55.));
    t = Transform3::translate(Vector3d(3.0, 10., 5.));
    EXPECT_EQ(t(Vector3d(1.0)), Vector3d(1.0, 1., 1.0));
    ASSERT_NEAR(t(Vector3d(1.0, 0.0, 0.)).length2(), 1., 1e-8);
    ASSERT_NEAR(t(Vector3d(0.0, 1.0, 0.)).length2(), 1., 1e-8);
    ASSERT_NEAR(t(Vector3d(0.0, 0.0, 1.)).length2(), 1., 1e-8);
  }
}

TEST(Transform, Scale) {
  {
    Transform2 t = Transform2::scale(Vector2d(10., -20.));
    EXPECT_EQ(t(Point2d(-2.0, 1.0)), Point2d(-20.0, -20.0));
    EXPECT_EQ(t.inverse()(Point2d(-20., -20.)), Point2d(-2., 1.));
    t = Transform2::scale(Vector2d(2.));
    EXPECT_EQ(t(Vector2d(1.0, 1.0)), Vector2d(2.0, 2.0));
    ASSERT_NEAR(t(Vector2d(1.0, 0.0)).length2(), 4.0, 1e-8);
  }
  {
    Transform3 t = Transform3::scale(Vector3d(10.0, -30.0, 5.));
    EXPECT_EQ(t(Point3d(-3.0, 1.0, 10.)), Point3d(-30.0, -30.0, 50.));
    EXPECT_EQ(t.inverse()(Point3d(-30., -30., -50.)), Point3d(-3., 1., -10.));
    t = Transform3::scale(Vector3d(3.0, 10., 5.));
    EXPECT_EQ(t(Vector3d(1.0)), Vector3d(3.0, 10., 5.0));
    ASSERT_NEAR(t(Vector3d(1.0, 0.0, 0.)).length2(), 9., 1e-8);
    ASSERT_NEAR(t(Vector3d(0.0, 1.0, 0.)).length2(), 100., 1e-8);
    ASSERT_NEAR(t(Vector3d(0.0, 0.0, 1.)).length2(), 25., 1e-8);
  }
}

TEST(Transform2D, Constructors) {
  Transform2D t;

  EXPECT_EQ(t.getMatrix().isIdentity(), true);
  EXPECT_EQ(t.inverse().getMatrix().isIdentity(), true);
  Matrix3d m;
  m.setIdentity();
  m(0, 0) = 2.0;
  m(1, 1) = 4.0;
  Transform2D t2(m, inverse(m));
  Point2d p = t2(Point2d(1.0, 1.0));
  EXPECT_EQ(p == Point2d(2.0, 4.0), true);
  EXPECT_EQ(Point2d(1.0, 1.0) == t2.inverse()(p), true);
  Transform2D b(BBox2D(Point2d(-1.0, -1.0), Point2d(2.0, 2.0)));
  EXPECT_EQ(Point2d(-1.0, -1.0), b.inverse()(Point2d(.0, .0)));
  EXPECT_EQ(Point2d(.0, .0), b(Point2d(-1.0, -1.0)));
  EXPECT_EQ(Point2d(2.0, 2.0), b.inverse()(Point2d(1.0, 1.0)));
  EXPECT_EQ(Point2d(1.0, 1.0), b(Point2d(2.0, 2.0)));
}

TEST(Transform2D, Translate) {
  Transform2D t;

  t.applyTranslate(Vector2d(10., -2.));
  EXPECT_EQ(t(Vector2d(0.0, 0.0)) == Vector2d(0.0, 0.0), true);
  EXPECT_EQ(t(Point2d(0.0, 0.0)) == Point2d(10.0, -2.0), true);
}

TEST(Transform2D, Scale) {
  Transform2D t;
  t.applyScale(10.0, -20.0);
  EXPECT_EQ(t(Point2d(-2.0, 1.0)) == Point2d(-20.0, -20.0), true);
  EXPECT_EQ(t.inverse()(Point2d(-20., -20.)) == Point2d(-2., 1.), true);
  t = Transform2D::scale(2.0, 2.0);
  EXPECT_EQ(t(Vector2d(1.0, 1.0)) == Vector2d(2.0, 2.0), true);
  ASSERT_NEAR(t(Vector2d(1.0, 0.0)).length2(), 4.0, 1e-8);
}

TEST(Transform2D, Rotate) {
  Transform2D t90 = Transform2D::rotate(90.f);
  Transform2D tt = t90; // test copy contructor

  EXPECT_EQ(t90(Vector2d(1.0, 0.0)) == Vector2d(0.0, 1.0), true);
  EXPECT_EQ(tt(Vector2d(1.0, 0.0)) == Vector2d(0.0, 1.0), true);
  Transform2D r90;
  r90.applyRotate(90.0);
  EXPECT_EQ(r90(Vector2d(1.0, 0.0)) == Vector2d(0.0, 1.0), true);
  r90.applyRotate(-90.0);
  EXPECT_EQ(r90(Vector2d(1.0, 0.0)) == Vector2d(1.0, 0.0), true);
}

TEST(Transform2D, Operators) {
  Transform2D t;

  t *= Transform2D::rotate(90.0);
  EXPECT_EQ(t(Vector2d(1.0, 0.0)) == Vector2d(0.0, 1.0), true);
  Transform2D r = t * t;
  EXPECT_EQ(r(Vector2d(1.0, 0.0)) == Vector2d(-1.0, 0.0), true);
  Transform2D s = Transform2D::scale(10.0, -20.0);
  BBox2D b = s(BBox2D(Point2d(-1.0, -1.0), Point2d(2.0, 2.0)));
  EXPECT_EQ(b.lower() == Point2d(-10.0, -40.0), true);
  EXPECT_EQ(b.upper() == Point2d(20.0, 20.0), true);
}

TEST(Transform2D, Inverse) {
  Transform2D t = Transform2D::scale(1.5, 1.5) * Transform2D::rotate(45.0) *
                  Transform2D::translate(Vector2d(-1.0, 2.0));
  Transform2D it = t.inverse();
  Transform2D r = Transform2D::translate(Vector2d(1.0, -2.0)) *
                  Transform2D::rotate(-45.0) *
                  Transform2D::scale(1.0 / 1.5, 1.0 / 1.5);
  Point2d p = t(Point2d(10.0, -10.0));

  EXPECT_EQ(it(p) == r(p), true);
}

TEST(Transform3D, Constructors) {
  Transform3D t;

  EXPECT_EQ(t.getMatrix().isIdentity(), true);
  EXPECT_EQ(t.inverse().getMatrix().isIdentity(), true);
  Matrix4d m;
  m.setIdentity();
  m(0, 0) = 3.0;
  m(1, 1) = 4.0;
  m(2, 2) = 5.0;
  Transform3D t2(m, inverse(m));
  Point3d p = t2(Point3d(1.0));
  EXPECT_EQ(p == Point3d(3.0, 4.0, 5.), true);
  EXPECT_EQ(Point3d(1.0) == t2.inverse()(p), true);
  Transform3D b(BBox3D(Point3d(-1.0), Point3d(3.0)));
  EXPECT_EQ(Point3d(-1.0), b.inverse()(Point3d(.0)));
  EXPECT_EQ(Point3d(.0), b(Point3d(-1.0)));
  EXPECT_EQ(Point3d(3.0), b.inverse()(Point3d(1.0)));
  EXPECT_EQ(Point3d(1.0), b(Point3d(3.0)));
}

TEST(Transform3D, Translate) {
  Transform3D t;

  t.applyTranslate(Vector3d(10., -3., 5.));
  EXPECT_EQ(t(Vector3d(0.0)) == Vector3d(0.0), true);
  EXPECT_EQ(t(Point3d(0.0)) == Point3d(10.0, -3.0, 5.), true);
}

TEST(Transform3D, Scale) {
  Transform3D t;
  t.applyScale(10.0, -30.0, 5.);
  EXPECT_EQ(t(Point3d(-3.0, 1.0, 10.)), Point3d(-30.0, -30.0, 50.));
  EXPECT_EQ(t.inverse()(Point3d(-30., -30., -50.)), Point3d(-3., 1., -10.));
  t = Transform3D::scale(3.0, 10., 5.);
  EXPECT_EQ(t(Vector3d(1.0)), Vector3d(3.0, 10., 5.0));
  ASSERT_NEAR(t(Vector3d(1.0, 0.0, 0.)).length2(), 9., 1e-8);
  ASSERT_NEAR(t(Vector3d(0.0, 1.0, 0.)).length2(), 100., 1e-8);
  ASSERT_NEAR(t(Vector3d(0.0, 0.0, 1.)).length2(), 25., 1e-8);
}

TEST(Transform3D, Rotate) { NOT_IMPLEMENTED(); }

TEST(Transform3D, Operators) {
  Transform3D s = Transform3D::scale(10.0, -30.0, 5.);
  BBox3D b = s(BBox3D(Point3d(-1.0), Point3d(3.0, 3.0, 3.)));
  EXPECT_EQ(b.lower(), Point3d(-10.0, -90.0, -5.));
  EXPECT_EQ(b.upper(), Point3d(30.0, 30.0, 15.));
}

TEST(Transform3D, Inverse) {
  Transform3D t = Transform3D::scale(1.5, 1.5, 1.5) *
                  Transform3D::rotate(45.0, Vector3d(0., 1., 0.)) *
                  Transform3D::translate(Vector3d(-1.0, 3.0, 2.0));
  Transform3D it = t.inverse();
  Transform3D r = Transform3D::translate(Vector3d(1.0, -3.0, -2.)) *
                  Transform3D::rotate(-45.0, Vector3d(0., 1., 0.)) *
                  Transform3D::scale(1.0 / 1.5, 1.0 / 1.5, 1.0 / 1.5);
  Point3d p = t(Point3d(10.0, -10.0, 10.));

  EXPECT_EQ(it(p), r(p));
}
