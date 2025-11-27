#include <cmath>
#include <functional>
#include <furoo.h>
#include <gtest/gtest.h>
#include <vector>

using namespace furoo;

TEST(MLS, Quadratic) {
  Point2d target(0.5234375, 0.609375);
  std::vector<Point2d> points = {
      Point2d(0.499, 0.595275),
      Point2d(0.506, 0.595275),
      Point2d(0.513, 0.595275),
      Point2d(0.520, 0.595275),
      Point2d(0.527, 0.595275),
      Point2d(0.534, 0.595275),
      Point2d(0.541, 0.595275),
      Point2d(0.506, 0.588275),
      Point2d(0.513, 0.588275),
      Point2d(0.520, 0.588275),
      Point2d(0.527, 0.588275),
      Point2d(0.534, 0.588275),
      Point2d(0.541, 0.588275),
  };
  std::vector<double> f(13, -0.1635);
  MLS<2> mls(Definitions::PolynomialType::QUADRATIC);
  EXPECT_NEAR(-0.1635, mls.interpolateAt(target, points, f), 1e-5);

}
TEST(MLS, ConstantBase) {
  MLS<2> mls(Definitions::PolynomialType::LINEAR);
  ShepardInterpolant2 interp;
  std::vector<Point2d> points;
  std::vector<double> f;
  BoxSampler sampler;
  Point2d target;

  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  points.emplace_back(Point2d({0.0, 0.0}));
  points.emplace_back(Point2d({0.1, 0.0}));
  points.emplace_back(Point2d({0.0, 0.1}));
  points.emplace_back(Point2d({0.1, 0.1}));
  for (auto p : points)
    f.emplace_back(_f(p));
  for (auto p : points) {
    double sol = _f(p);
    double approx = interp.interpolateAt(p, points, f);
    double mlsApprox = mls.interpolateAt(p, points, f);
    ASSERT_NEAR(approx, sol, 1e-8);
    EXPECT_DOUBLE_EQ(approx, mlsApprox);
  }
  for (size_t i = 0; i < 100; i++) {
    auto p = sampler.sample() * 0.2 - 0.1;
    double sol = _f(p);
    double approx = interp.interpolateAt(p, points, f);
    double mlsApprox = mls.interpolateAt(p, points, f);
    ASSERT_NEAR(approx, sol, 2e-2);
    EXPECT_NEAR(approx, mlsApprox, 2e-2);
  }
}

TEST(RBFInterpolant, GaussianInterpolation) {
  RBFInterpolant<Point2d> rbf;
  std::vector<Point2d> points;
  std::vector<double> f;
  Point2d target;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::exp(p[0]) * std::cos(p[1]);
  };

  points.emplace_back(Point2d({0.0, 0.0}));
  points.emplace_back(Point2d({0.25, 0.0}));
  points.emplace_back(Point2d({0.25, 0.25}));
  points.emplace_back(Point2d({0.0, 0.25}));
  target = Point2d({0.1, 0.1});

  for (auto p : points) {
    f.emplace_back(_f(p));
    std::cout << _f(p) << '\n';
  }

  // std::cout << "RBFError ";
  // std::cout << rbf.interpolateAt(target, points, f) - _f(target) << '\n';
  ASSERT_NEAR(rbf.interpolateAt(target, points, f), _f(target), 1e-2);

  BoxSampler sampler;
  BilinearInterpolant bl;
  double maxerr = 0;
  for (size_t i = 0; i < 10; i++) {
    Point2d target =
        sampler.sample(BBox2D(Point2d(0., 0.), Point2d(0.25, 0.25)));
    double solution = _f(target);
    double rbfE, blE, blS;
    double rbfS = rbf.interpolateAt(target, points, f);
    blS = bl.interpolateAt(target, points, f);

    maxerr = max(maxerr, std::abs(rbfS - _f(target)));

    rbfE = rbfS - solution;
    blE = blS - solution;

    ASSERT_NEAR(rbfS, solution, 1e-2);
    EXPECT_LT(std::fabs(rbfE), std::fabs(blE));
  }
  std::cout << "err: " << maxerr << '\n';
}

TEST(BilinearInterpolant, interpolateAt) {
  BilinearInterpolant interp;
  std::vector<Point2d> points;
  std::vector<double> f;
  Point2d target;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };

  points.emplace_back(Point2d({0.0, 0.0}));
  points.emplace_back(Point2d({0.1, 0.0}));
  points.emplace_back(Point2d({0.0, 0.1}));
  points.emplace_back(Point2d({0.1, 0.1}));
  target = Point2d({0.05, 0.05});

  for (auto p : points) {
    f.emplace_back(_f(p));
  }
  std::cout << "Bilinear Error ";
  std::cout << interp.interpolateAt(target, points, f) - _f(target) << '\n';
  ASSERT_NEAR(interp.interpolateAt(target, points, f), _f(target), 5e-3);
}

TEST(BilinearInterpolant, weights) {
  BoxSampler sampler;
  BilinearInterpolant bil;
  std::vector<double> W(4, 0.);
  std::function<double(Point2d)> f = [](Point2d p) -> double {
    return std::sin(0.1 * p.x()) * std::cos(0.1 * p.y());
  };
  for (int i = 0; i < 1000; i++) {
    auto target = sampler.sample();
    bil.weights(target, std::vector<Point2d>(), W);
    double w[4];
    w[0] = (1 - target.x()) * (1 - target.y());
    w[1] = target.x() * (1 - target.y());
    w[2] = target.x() * target.y();
    w[3] = (1 - target.x()) * target.y();
    double wsum = W[0] + W[1] + W[2] + W[3];
    ASSERT_NEAR(wsum, 1., 1e-8);
    for (int j = 0; j < 4; j++)
      ASSERT_NEAR(W[j], w[j], 1e-8);
    ASSERT_NEAR((W[0] * f(Point2d(0, 0)) + W[1] * f(Point2d(1, 0)) +
        W[2] * f(Point2d(1, 1)) + W[3] * f(Point2d(0, 1))) /
        (wsum),
                f(target), 1e-2);
  }
}

TEST(BicubicInterpolant, interpolateAt) {
  BilinearInterpolant bil;
  BicubicInterpolant bic;
  std::vector<double> f;
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> fxy;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fx = [](Point2d p) {
    return -std::sin(p[0]) * std::cos(p[1]);
  };
  std::function<double(Point2d)> _fy = [](Point2d p) {
    return -std::cos(p[0]) * std::sin(p[1]);
  };
  std::function<double(Point2d)> _fxy = [](Point2d p) {
    return std::sin(p[0]) * std::sin(p[1]);
  };
  std::vector<Point2d> points = {Point2d(0., 0.), Point2d(0.1, 0.),
                                 Point2d(0.1, 0.1), Point2d(0., 0.1)};
  for (auto p : points) {
    f.emplace_back(_f(p));
    fx.emplace_back(_fx(p));
    fy.emplace_back(_fy(p));
    fxy.emplace_back(_fxy(p));
  }
  BoxSampler sampler;
  for (size_t i = 0; i < 1000; i++) {
    Point2d target = sampler.sample(BBox2D::make_unit_bbox());
    ASSERT_NEAR(bic.interpolateAt(target, f, fx, fy, fxy, Vector2d(0.1, 0.1)),
                _f(target * 0.1), 1e-6);
    EXPECT_GT(
        fabs(bil.interpolateAt(target, points, f) - _f(target * 0.1)),
        fabs(bic.interpolateAt(target, f, fx, fy, fxy, Vector2d(0.1, 0.1)) -
            _f(target * 0.1)));
  }
}

typedef Point<double, 1> Point1d;
TEST(TrilinearHatKernel, Linear1D) {
  KernelInterpolant<Point1d> k(new TrilinearHatKernel<Point1d>(0.5));
  std::vector<Point1d> points;
  std::vector<double> f;
  auto region = BBox2D::make_unit_bbox();
  BoxSampler sampler;

  for (size_t i = 0; i < 100; i++) {
    points.emplace_back(Point1d(sampler.sample(region)[0]));
  }
  for (auto p : points) {
    f.emplace_back(p[0]);
  }

  for (size_t i = 0; i < 100; i++) {
    EXPECT_EQ(points[i], f[i]);
  }

  std::cout << k.interpolateAt(Point1d(0.0), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.1), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.2), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.3), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.4), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.5), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.6), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.7), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.8), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.9), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(1.0), points, f) << ' ' << '\n';
}

TEST(QuadraticBSplineKernel, Linear1D) {
  KernelInterpolant<Point1d> k(new QuadraticBSplineKernel<Point1d>(1));
  std::vector<Point1d> points;
  std::vector<double> f;

  points.emplace_back(Point1d(0.1));
  // points.emplace_back(Point1d(.5));
  points.emplace_back(Point1d(1.1));

  f.emplace_back(0.1);
  // f.emplace_back(0.5);
  f.emplace_back(1.1);

  std::cout << k.interpolateAt(Point1d(0.0), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.1), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.2), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.3), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.4), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.5), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.6), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.7), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.8), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(0.9), points, f) << ' ' << '\n';
  std::cout << k.interpolateAt(Point1d(1.0), points, f) << ' ' << '\n';
}

TEST(ShepardInterpolant2, interpolateAt) {
  ShepardInterpolant2 interp;
  std::vector<Point2d> points;
  std::vector<double> f;
  BoxSampler sampler;
  Point2d target;

  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0]) * std::cos(p[1]);
  };
  points.emplace_back(Point2d({0.0, 0.0}));
  points.emplace_back(Point2d({0.1, 0.0}));
  points.emplace_back(Point2d({0.0, 0.1}));
  points.emplace_back(Point2d({0.1, 0.1}));
  for (auto p : points) {
    f.emplace_back(_f(p));
    std::cout << _f(p) << '\n';
  }
  std::cout << '\n';
  for (auto p : points) {
    double sol = _f(p);
    double approx = interp.interpolateAt(p, points, f);
    ASSERT_NEAR(approx, sol, 1e-8);
  }
  for (size_t i = 0; i < 1e2; i++) {
    auto p = sampler.sample() * 0.2 - 0.1;
    double sol = _f(p);
    double approx = interp.interpolateAt(p, points, f);
    std::cout << approx << '\n';
    // std::cout << approx << " " << sol << '\n';
    // std::cout << "Error: " << abs(approx - sol) << '\n';
    ASSERT_NEAR(approx, sol, 2e-2);
  }
}
TEST(MonotonicCatmullRom, interpolate) {
  double a = 0.0;
  double b = 0.0;
  double c = 1.0;
  double d = 1.0;

  for (int i = 0; i <= 10; ++i) {
    double result = MonotonicCatmullRom::interpolate(a, b, c, d, i * 0.1);
    EXPECT_TRUE(result >= b && result <= c);

    if (i == 0) {
      EXPECT_DOUBLE_EQ(b, result);
    } else if (i == 10) {
      EXPECT_DOUBLE_EQ(c, result);
    }
  }

  a = 0.0;
  b = 1.0;
  c = 2.0;
  d = 3.0;

  for (int i = 0; i <= 10; ++i) {
    double result = MonotonicCatmullRom::interpolate(a, b, c, d, i * 0.1);
    EXPECT_TRUE(result >= b && result <= c);

    if (i == 0) {
      EXPECT_DOUBLE_EQ(b, result);
    } else if (i == 10) {
      EXPECT_DOUBLE_EQ(c, result);
    }
  }

  a = 0.0;
  b = 1.0;
  c = 2.0;
  d = 0.0;

  for (int i = 0; i <= 10; ++i) {
    double result = MonotonicCatmullRom::interpolate(a, b, c, d, i * 0.1);
    EXPECT_TRUE(result >= b && result <= c);

    if (i == 0) {
      EXPECT_DOUBLE_EQ(b, result);
    } else if (i == 10) {
      EXPECT_DOUBLE_EQ(c, result);
    }
  }

  a = 0.0;
  b = 2.0;
  c = 1.0;
  d = 3.0;

  for (int i = 0; i <= 10; ++i) {
    double result = MonotonicCatmullRom::interpolate(a, b, c, d, i * 0.1);
    EXPECT_TRUE(result >= c && result <= b);

    if (i == 0) {
      EXPECT_DOUBLE_EQ(b, result);
    } else if (i == 10) {
      EXPECT_DOUBLE_EQ(c, result);
    }
  }
}
