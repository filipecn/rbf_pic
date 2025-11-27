#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Point, Constructors) {
  {
    Point3d v;
    ASSERT_NEAR(0.0, v[0], 1e-8);
    ASSERT_NEAR(0.0, v[1], 1e-8);
    ASSERT_NEAR(0.0, v[2], 1e-8);
  }
  {
    Point2d v;
    ASSERT_NEAR(0.0, v[0], 1e-8);
    ASSERT_NEAR(0.0, v[1], 1e-8);
  }
  {
    Point3d v({0.1, 0.2, 0.3});
    ASSERT_NEAR(0.1, v[0], 1e-8);
    ASSERT_NEAR(0.2, v[1], 1e-8);
    ASSERT_NEAR(0.3, v[2], 1e-8);
  }
  {
    Point2d v({1.1, 1.2});
    ASSERT_NEAR(1.1, v[0], 1e-8);
    ASSERT_NEAR(1.2, v[1], 1e-8);
  }
  {
    Point3d v(0.1);
    ASSERT_NEAR(0.1, v[0], 1e-8);
    ASSERT_NEAR(0.1, v[1], 1e-8);
    ASSERT_NEAR(0.1, v[2], 1e-8);
  }
  {
    Point2d v(1.1);
    ASSERT_NEAR(1.1, v[0], 1e-8);
    ASSERT_NEAR(1.1, v[1], 1e-8);
  }
  {
    Point3i v;
    EXPECT_EQ(0, v[0]);
    EXPECT_EQ(0, v[1]);
    EXPECT_EQ(0, v[2]);
  }
  {
    Point2i v;
    EXPECT_EQ(0, v[0]);
    EXPECT_EQ(0, v[1]);
  }
  {
    Point3i v({1, 2, 3});
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
    EXPECT_EQ(3, v[2]);
  }
  {
    Point2i v({1, 2});
    EXPECT_EQ(1, v[0]);
    EXPECT_EQ(2, v[1]);
  }
  {
    Point3i v(8);
    EXPECT_EQ(8, v[0]);
    EXPECT_EQ(8, v[1]);
    EXPECT_EQ(8, v[2]);
  }
  {
    Point2i v(5);
    EXPECT_EQ(5, v[0]);
    EXPECT_EQ(5, v[1]);
  }

  Point<int, 10> a({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  Point<int, 10> b = a;
  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(a[i], b[i]);
  }
}

TEST(Point, Boolean) {
  Point3i a({1, 2, 3});
  Point3i b({1, 2, 3});
  Point3i c({1, 0, 0});
  Point3i d({10, 10, 10});

  EXPECT_EQ(a == b, true);
  EXPECT_EQ(a != b, false);
  EXPECT_EQ(a == c, false);
  EXPECT_EQ(a != c, true);
  EXPECT_EQ(c <= a, true);
  EXPECT_EQ(a >= c, true);
  EXPECT_EQ(d > a, true);
  EXPECT_EQ(b > a, false);
  EXPECT_EQ(a < d, true);
  EXPECT_EQ(c < a, false);
}

TEST(Point, Arithmetic) {
  EXPECT_EQ(Point2i({1, 1}) - Point2i({2, -3}), Vector2i({-1, 4}));
  EXPECT_EQ(Point2i({0, 0}) + Point2i({2, -3}), Point2i({2, -3}));
  // EXPECT_EQ( Point2i( { 1, 2 } ) * Point2i( { 2, -3 } ), Point2i( { 2, -6 } )
  // );
  EXPECT_EQ(Point2i({1, 2}) * 2, Point2i({2, 4}));
  // EXPECT_EQ( Point2i( { 10, 20 } ) / Point2i( { 2, -2 } ), Point2i( { 5, -10
  // } ) );
  Point2d a({5.0, 10.0});
  EXPECT_EQ(a / 2.0, Point2d({2.5, 5.0}));
  EXPECT_EQ(a * 3.0, Point2d({15.0, 30.0}));
  a /= 2.0;
  EXPECT_EQ(a, Point2d({2.5, 5.0}));
  a -= 2.0;
  EXPECT_EQ(a, Point2d({0.5, 3.0}));
  EXPECT_EQ(Point2i(10, 5) / 2, Point2i(5, 2));
}
//
// TEST( Point, Gets ) {
//  EXPECT_EQ( Point3i( { 1, 2, 3 } ).xy(), Point2i( { 1, 2 } ) );
//  EXPECT_EQ( Point3i( { 1, 2, 3 } ).xy( 1, 0 ), Point2i( { 2, 1 } ) );
//  Point< double, 10 > a( { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 }
//  );
//  EXPECT_EQ( a.doubleXY( 5, 9 ), Point2d( 6.0, 10.0 ) );
//  a[5] = 60.0;
//  EXPECT_EQ( a.doubleXY( 0, 5 ), Point2d( 1.0, 60.0 ) );
//  EXPECT_EQ( Point3i( 1, 2, 3 ).doubleXY(), Point2d( 1.0, 2.0 ) );
//  EXPECT_EQ( Point3i( 1, 2, 3 ).doubleXY( 1, 0 ), Point2d( 2.0, 1.0 ) );
//  Point< char, 4 > b( { 'a', 'b', 'c', 'd' } );
//  EXPECT_EQ( b.doubleXYZ(), Point3d( { 97.0, 98.0, 99.0 } ) );
//  EXPECT_EQ( b.doubleXYZ( 3, 0, 2 ), Point3d( { 100.0, 97.0, 99.0 } ) );
// }
//
TEST(Point, Analytic) {
  //  Point< double, 10 > a(
  //      { 1.0, 254.0, 3.0, 4.0, 50.0, 6.0, 7.0, -823.0, 9.0, 10.0 } );
  //  ASSERT_NEAR( a.max(), 254.0, 1e-8 );
  ASSERT_NEAR(Point2i({3, 4}).distance2(Point2i({0, 0})), 25.0, 1e-8);
  ASSERT_NEAR(Point2i({3, 4}).distance(Point2i({0, 0})), 5.0, 1e-8);
  ASSERT_NEAR(distance2(Point2i({3, 4}), Point2i({0, 0})), 25.0, 1e-8);
  ASSERT_NEAR(distance(Point2i({3, 4}), Point2i({0, 0})), 5.0, 1e-8);
  //  EXPECT_EQ( Point2i( { 3, 4 } ).normalized(), Point2d( { 3.0, 4.0 } ) / 5.0
  //  );
  //  EXPECT_EQ( Point2i( { 3, 4 } ).left(), Point2i( { -4, 3 } ) );
  //  EXPECT_EQ( Point2i( { 3, 4 } ).right(), Point2i( { 4, -3 } ) );
}

// TEST( Point, Functions ) {
// Point2d p;
//  EXPECT_EQ( 2 * Point2i( { 1, 5 } ), Point2i( { 2, 10 } ) );
//  EXPECT_EQ( round( Point2d( { -10.9, 2.01 } ) ), Point2i( { -11, 2 } ) );
//  EXPECT_EQ( round( Point3d( { -10.499, 2.51, 1.0 } ) ), Point3i( { -10, 3, 1
//  } ) );
//  EXPECT_EQ( ceil( Point2d( { -10.9, 2.01 } ) ), Point2i( { -9, 3 } ) );
//  EXPECT_EQ( ceil( Point3d( { -10.499, 2.51, 1.0 } ) ), Point3i( { -9, 3, 2 }
//  ) );
//  EXPECT_EQ( floor( Point2d( { -10.9, 2.01 } ) ), Point2i( { -10, 2 } ) );
//  EXPECT_EQ( floor( Point3d( { -10.499, 2.51, 1.0 } ) ), Point3i( { -10, 2, 1
//  } ) );
//  EXPECT_EQ( min( Point3i( { 0, 4, -1 } ), Point3i( { -1, 8, 0 } ) ), Point3i(
//  { -1, 4, -1 } ) );
//  EXPECT_EQ( max( Point3i( { 0, 4, -1 } ), Point3i( { -1, 8, 0 } ) ), Point3i(
//  { 0, 8, 0 } ) );
//  EXPECT_EQ( min( Point2i( { 0, 4 } ), Point2i( { -1, 8 } ) ), Point2i( { -1,
//  4 } ) );
//  EXPECT_EQ( max( Point2i( { 0, 4 } ), Point2i( { -1, 8 } ) ), Point2i( { 0, 8
//  } ) );
//  Point2d a( { 3.0, 4.0 } );
//  normalize( a );
//  EXPECT_EQ( a, Point2d( { 3.0, 4.0 } ) / 5.0 );
// }
