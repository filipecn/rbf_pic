#include <cmath>
#include <furoo.h>
#include <gtest/gtest.h>
#include <vector>

using namespace furoo;

TEST(RegularGrid, Constructors) {
  {
    RegularGrid grid;
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 0u);
    EXPECT_EQ(size.y(), 0u);
    EXPECT_EQ(grid.cellCount(), 0u);
  }
  {
    RegularGrid grid(10);
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++)
        EXPECT_EQ(grid.fieldPosition(i, j) == Point2d(i, j), true);
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    EXPECT_EQ(grid.cellCount(), 100u);
  }
  {
    RegularGrid grid(10, 25.);
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    double h = 25. / 10.;
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++) {
        EXPECT_EQ(grid.fieldPosition(i, j) ==
                      Point2d(i * h + .5 * h, j * h + .5 * h),
                  true);
        EXPECT_EQ(grid.fieldPosition(i, j) ==
                      grid.transformToWorld(Point2d(i, j)),
                  true);
        EXPECT_EQ(grid.transformToGrid(
                      Point2d(i * h + .5 * h, j * h + .5 * h)) == Point2d(i, j),
                  true);
      }
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    EXPECT_EQ(grid.cellCount(), 100u);
  }
  {
    RegularGrid grid(20, Transform2D::scale(1.5, 2.));
    auto size = grid.gridSize();
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 20; j++) {
        EXPECT_EQ(grid.fieldPosition(i, j) == Point2d(i * 1.5, j * 2.), true);
      }
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 20u);
    EXPECT_EQ(grid.cellCount(), 400u);
  }
  {
    RegularGrid grid(20, 40, Transform2D::translate(Vector2d(1.5, 2.)));
    auto size = grid.gridSize();
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 40; j++) {
        EXPECT_EQ(grid.fieldPosition(i, j) == Point2d(i + 1.5, j + 2.), true);
      }
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 40u);
    EXPECT_EQ(grid.cellCount(), 800u);
  }
  {
    RegularGrid grid(20, 40, 0.1, 0.2);
    auto size = grid.gridSize();
    double hx = 0.1 / 20;
    double hy = 0.2 / 40;
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 40; j++) {
        EXPECT_EQ(grid.fieldPosition(i, j) ==
                      Point2d(i * hx + hx * .5, j * hy + hy * .5),
                  true);
      }
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 40u);
    EXPECT_EQ(grid.cellCount(), 800u);
  }
  {
    RegularGrid grid(10, 10, BBox2D(Point2d(-2., 4.), Point2d(2., 6.)));
    EXPECT_EQ(grid.transformToWorld(Point2d(1, 0)), Point2d(-1.6, 4.));
    EXPECT_EQ(grid.transformToWorld(Point2d(0, 0)), Point2d(-2., 4.));
    EXPECT_EQ(grid.transformToWorld(Point2d(0, 1)), Point2d(-2., 4.2));
  }
}

TEST(RegularGrid3, Constructors) {
  {
    RegularGrid3 grid;
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 0u);
    EXPECT_EQ(size.y(), 0u);
    EXPECT_EQ(size.z(), 0u);
    EXPECT_EQ(grid.cellCount(), 0u);
  }
  {
    RegularGrid3 grid(10);
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    EXPECT_EQ(size.z(), 10u);
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++)
        for (size_t k = 0; k < 10; k++)
          EXPECT_EQ(grid.fieldPosition(i, j, k), Point3d(i, j, k));
    EXPECT_EQ(grid.cellCount(), 1000u);
  }
  {
    RegularGrid3 grid(10, 25.);
    auto size = grid.gridSize();
    EXPECT_EQ(size.x(), 10u);
    EXPECT_EQ(size.y(), 10u);
    EXPECT_EQ(size.z(), 10u);
    double h = 25. / 10.;
    for (size_t i = 0; i < 10; i++)
      for (size_t j = 0; j < 10; j++)
        for (size_t k = 0; k < 10; k++) {
          EXPECT_EQ(grid.fieldPosition(i, j, k),
                    Point3d(i * h + .5 * h, j * h + .5 * h, k * h + .5 * h));
          EXPECT_EQ(grid.fieldPosition(i, j, k),
                    grid.transformToWorld(Point3d(i, j, k)));
          EXPECT_EQ(grid.transformToGrid(Point3d(i * h + .5 * h, j * h + .5 * h,
                                                 k * h + .5 * h)),
                    Point3d(i, j, k));
        }
    EXPECT_EQ(grid.cellCount(), 1000u);
  }
  {
    RegularGrid3 grid(20, Transform3D::scale(1.5, 2., 10.));
    auto size = grid.gridSize();
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 20; j++)
        for (size_t k = 0; k < 20; k++)
          EXPECT_EQ(grid.fieldPosition(i, j, k),
                    Point3d(i * 1.5, j * 2., k * 10.));
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 20u);
    EXPECT_EQ(size.z(), 20u);
    EXPECT_EQ(grid.cellCount(), 8000u);
  }
  {
    RegularGrid3 grid(20, 40, 100,
                      Transform3D::translate(Vector3d(1.5, 2., 5.)));
    auto size = grid.gridSize();
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 40; j++)
        for (size_t k = 0; k < 100; k++)
          EXPECT_EQ(grid.fieldPosition(i, j, k),
                    Point3d(i + 1.5, j + 2., k + 5.));
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 40u);
    EXPECT_EQ(size.z(), 100u);
    EXPECT_EQ(grid.cellCount(), 80000u);
  }
  {
    RegularGrid3 grid(20, 40, 100, 0.1, 0.2, 0.3);
    auto size = grid.gridSize();
    double hx = 0.1 / 20;
    double hy = 0.2 / 40;
    double hz = 0.3 / 100;
    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 40; j++)
        for (size_t k = 0; k < 100; k++)
          EXPECT_EQ(
              grid.fieldPosition(i, j, k),
              Point3d(i * hx + hx * .5, j * hy + hy * .5, k * hz + hz * .5));
    EXPECT_EQ(size.x(), 20u);
    EXPECT_EQ(size.y(), 40u);
    EXPECT_EQ(size.z(), 100u);
    EXPECT_EQ(grid.cellCount(), 80000u);
  }
  {
    RegularGrid3 grid(10, 10, 10,
                      BBox3D(Point3d(-2., 4., -10), Point3d(2., 6., 10.)));
    EXPECT_EQ(grid.transformToWorld(Point3d(1, 0, 0)), Point3d(-1.6, 4., -10.));
    EXPECT_EQ(grid.transformToWorld(Point3d(0, 0, 0)), Point3d(-2., 4., -10.));
    EXPECT_EQ(grid.transformToWorld(Point3d(0, 1, 0)), Point3d(-2., 4.2, -10.));
    EXPECT_EQ(grid.transformToWorld(Point3d(0, 0, 1)), Point3d(-2., 4., -8.));
  }
}
