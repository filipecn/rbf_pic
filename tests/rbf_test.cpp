#include <functional>
#include <furoo.h>
#include <gtest/gtest.h>
#include <vector>

using namespace furoo;

TEST(GaussianKernel, phi) {
  GaussianKernel<Point2d> gk;
  double exp1 = std::exp(-1.0);
  ASSERT_NEAR(gk.getH(), 1.0, 1e-8);
  ASSERT_NEAR(gk.phi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), 1.0, 1e-8);
  ASSERT_NEAR(gk.dphi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), -2.0, 1e-8);
  ASSERT_NEAR(gk.d2phi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), -2.0, 1e-8);

  ASSERT_NEAR(gk.phi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), exp1, 1e-8);
  ASSERT_NEAR(gk.dphi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), -2.0 * exp1, 1e-8);
  ASSERT_NEAR(gk.d2phi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), 2.0 * exp1, 1e-8);
}

TEST(IMQKernel, phis) {
  IMQKernel<Point2d> imq;
  double sq = std::sqrt(2.0);
  ASSERT_NEAR(imq.getH(), 1.0, 1e-8);
  ASSERT_NEAR(imq.phi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), 1.0, 1e-8);
  ASSERT_NEAR(imq.dphi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), -1.0, 1e-8);
  ASSERT_NEAR(imq.d2phi(Point2d(0.0, 0.0), Point2d(0.0, 0.0)), -1.0, 1e-8);

  ASSERT_NEAR(imq.phi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), 1 / sq, 1e-8);
  ASSERT_NEAR(imq.dphi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), -0.5 / sq, 1e-8);
  ASSERT_NEAR(imq.d2phi(Point2d(0.0, 0.0), Point2d(1.0, 0.0)), 0.25 / sq, 1e-8);
}

TEST(RBF, GaussianKernel) {
  RBF<Point2d> rbf(new GaussianKernel<Point2d>());
  std::vector<Point2d> points;

  points.emplace_back(Point2d(0.5, 0.0));
  points.emplace_back(Point2d(1.5, 0.0));
  points.emplace_back(Point2d(-0.5, 0.0));
  points.emplace_back(Point2d(0.5, 1.0));
  points.emplace_back(Point2d(0.5, -1.0));

  LinearVector w = rbf.laplacian(points);

  ASSERT_NEAR(w[0], -6.820262746251077, 1e-8);
  ASSERT_NEAR(w[1], 1.705065686562769, 1e-8);
  ASSERT_NEAR(w[2], 1.705065686562769, 1e-8);
  ASSERT_NEAR(w[3], 1.705065686562769, 1e-8);
  ASSERT_NEAR(w[4], 1.705065686562769, 1e-8);

  w = rbf.gradientX(points);
  ASSERT_NEAR(w[0], 0.0, 1e-8);
  ASSERT_NEAR(w[1], 0.5, 1e-8);
  ASSERT_NEAR(w[2], -0.5, 1e-8);
  ASSERT_NEAR(w[3], 0.0, 1e-8);
  ASSERT_NEAR(w[4], 0.0, 1e-8);

  w = rbf.gradientY(points);
  ASSERT_NEAR(w[0], 0.0, 1e-8);
  ASSERT_NEAR(w[1], 0.0, 1e-8);
  ASSERT_NEAR(w[2], 0.0, 1e-8);
  ASSERT_NEAR(w[3], 0.5, 1e-8);
  ASSERT_NEAR(w[4], -0.5, 1e-8);
}

TEST(RBF, GradientLinearFunc) {
  RBF<Point2d> rbf(new GaussianKernel<Point2d>());
  std::vector<Point2d> points;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return p[0] + p[1]; // Linear function
    return std::cos(p[0] * 3) * std::sin(p[1] * 3);
  };
  std::function<double(Point2d)> _gx = [](Point2d p) {
    UNUSED_VARIABLE(p);
    return 1.0; // Linear function
  };
  std::function<double(Point2d)> _gy = [](Point2d p) {
    UNUSED_VARIABLE(p);
    return 1.0; // Linear function
  };

  points.emplace_back(Point2d(0.0, 0.0));
  points.emplace_back(Point2d(0.5, 0.0));
  points.emplace_back(Point2d(1.0, 0.0));
  points.emplace_back(Point2d(1.0, 0.5));
  points.emplace_back(Point2d(1.0, 1.0));
  points.emplace_back(Point2d(0.5, 1.0));
  points.emplace_back(Point2d(0.0, 1.0));
  points.emplace_back(Point2d(0.0, 0.5));

  std::vector<Point2d> testPoints;
  testPoints.emplace_back(Point2d({0.5, 0.5}));
  testPoints.emplace_back(Point2d({0.4, 0.0}));
  testPoints.emplace_back(Point2d({0.0, 0.1}));
  testPoints.emplace_back(Point2d({0.45, 0.12}));
  testPoints.emplace_back(Point2d({0.35, 0.97}));
  testPoints.emplace_back(Point2d({0.78, 0.81}));
  testPoints.emplace_back(Point2d({0.78, 0.14}));

  for (auto p : testPoints) {
    LinearVector wx = rbf.gradientXAt(p, points);
    LinearVector wy = rbf.gradientYAt(p, points);
    // std::cout << "-------" << '\n';
    // std::cout << wx;
    // std::cout << wy;
    // std::cout << p;
    double gx = 0;
    double gy = 0;
    for (size_t i = 0; i < points.size(); i++) {
      //   std::cout << wx[i] << " * " << _f(points[i]) << '\n';
      gx += wx[i] * _f(points[i]);
      gy += wy[i] * _f(points[i]);
    }
    // std::cout << '\n';
    // std::cout << "Gx: " << gx << "\tGy:" << gy << '\n';
    // std::cout << "Gx: " << _gx(p) << "\tGy:" << _gy(p) << '\n';
    ASSERT_NEAR(gx, _gx(p), 1e-2);
    ASSERT_NEAR(gy, _gy(p), 1e-2);
    // std::cout << "OK" << '\n';
  }
}

TEST(RBF, GradientSphereFunc) {
  RBF<Point2d> rbf(new GaussianKernel<Point2d>(), 0.125);
  std::vector<Point2d> points;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return p[0] * p[0] + p[1] * p[1]; // x-only
  };
  std::function<double(Point2d)> _gx = [](Point2d p) {
    return 2 * p[0]; // x-only
  };
  std::function<double(Point2d)> _gy = [](Point2d p) {
    return 2 * p[1]; // x-only
  };

  points.emplace_back(Point2d(0.0, 0.0));
  points.emplace_back(Point2d(0.5, 0.0));
  points.emplace_back(Point2d(1.0, 0.0));
  points.emplace_back(Point2d(1.0, 0.5));
  points.emplace_back(Point2d(1.0, 1.0));
  points.emplace_back(Point2d(0.5, 1.0));
  points.emplace_back(Point2d(0.0, 1.0));
  points.emplace_back(Point2d(0.0, 0.5));

  std::vector<Point2d> testPoints;
  testPoints.emplace_back(Point2d({0.5, 0.5}));
  testPoints.emplace_back(Point2d({0.4, 0.0}));
  testPoints.emplace_back(Point2d({0.0, 0.1}));
  testPoints.emplace_back(Point2d({0.45, 0.12}));
  testPoints.emplace_back(Point2d({0.35, 0.97}));
  testPoints.emplace_back(Point2d({0.78, 0.81}));
  testPoints.emplace_back(Point2d({0.78, 0.14}));

  for (auto p : testPoints) {
    LinearVector wx = rbf.gradientXAt(p, points);
    LinearVector wy = rbf.gradientYAt(p, points);
    // std::cout << "-------" << '\n';
    // std::cout << wx;
    // std::cout << wy;
    // std::cout << p;
    double gx = 0;
    double gy = 0;
    for (size_t i = 0; i < points.size(); i++) {
      //   std::cout << wx[i] << " * " << _f(points[i]) << '\n';
      gx += wx[i] * _f(points[i]);
      gy += wy[i] * _f(points[i]);
    }
    // std::cout << '\n';
    // std::cout << "Gx: " << gx << "\tGy:" << gy << '\n';
    // std::cout << "Gx: " << _gx(p) << "\tGy:" << _gy(p) << '\n';
    ASSERT_NEAR(gx, _gx(p), 1e-2);
    ASSERT_NEAR(gy, _gy(p), 1e-2);
    // std::cout << "OK" << '\n';
  }
}

TEST(RBF, GradientTrigonometricFunc) {
  RBF<Point2d> rbf(new GaussianKernel<Point2d>());
  std::vector<Point2d> points;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::cos(p[0] * 3) * std::sin(p[1] * 3);
  };
  std::function<double(Point2d)> _gx = [](Point2d p) {
    return -3 * std::sin(3 * p[0]) * std::sin(3 * p[1]);
  };
  std::function<double(Point2d)> _gy = [](Point2d p) {
    return 3 * std::cos(3 * p[0]) * std::cos(3 * p[1]);
  };

  points.emplace_back(Point2d(0.0, 0.0));
  points.emplace_back(Point2d(0.25, 0.0));
  points.emplace_back(Point2d(0.5, 0.0));
  points.emplace_back(Point2d(0.5, 0.25));
  points.emplace_back(Point2d(0.5, 0.5));
  points.emplace_back(Point2d(0.25, 0.5));
  points.emplace_back(Point2d(0.0, 0.5));
  points.emplace_back(Point2d(0.0, 0.25));
  for (size_t i = 0; i < points.size(); i++) {
    points[i] /= 4;
  }

  std::vector<Point2d> testPoints;
  testPoints.emplace_back(Point2d({0.25, 0.25}));
  testPoints.emplace_back(Point2d({0.2, 0.0}));
  testPoints.emplace_back(Point2d({0.0, 0.05}));
  testPoints.emplace_back(Point2d({0.27, 0.06}));
  testPoints.emplace_back(Point2d({0.17, 0.48}));
  testPoints.emplace_back(Point2d({0.34, 0.40}));
  testPoints.emplace_back(Point2d({0.34, 0.07}));
  for (size_t i = 0; i < testPoints.size(); i++) {
    testPoints[i] /= 4;
  }

  for (auto p : testPoints) {
    LinearVector wx = rbf.gradientXAt(p, points);
    LinearVector wy = rbf.gradientYAt(p, points);
    // std::cout << "-------" << '\n';
    // std::cout << wx;
    // std::cout << wy;
    // std::cout << p;
    double gx = 0;
    double gy = 0;
    for (size_t i = 0; i < points.size(); i++) {
      //   std::cout << wx[i] << " * " << _f(points[i]) << '\n';
      gx += wx[i] * _f(points[i]);
      gy += wy[i] * _f(points[i]);
    }
    // std::cout << '\n';
    // std::cout << "Gx: " << gx << "\tGy:" << gy << '\n';
    // std::cout << "Gx: " << _gx(p) << "\tGy:" << _gy(p) << '\n';
    ASSERT_NEAR(gx, _gx(p), 1e-2);
    ASSERT_NEAR(gy, _gy(p), 1e-2);
    // std::cout << "OK" << '\n';
  }
}

TEST(RBF, GradientExponentialFunc) {
  RBF<Point2d> rbf(new GaussianKernel<Point2d>());
  std::vector<Point2d> points;
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::exp(p[0]) + std::exp(p[1]); // x-only
  };
  std::function<double(Point2d)> _gx = [](Point2d p) {
    return std::exp(p[0]); // x-only
  };
  std::function<double(Point2d)> _gy = [](Point2d p) {
    return std::exp(p[1]); // x-only
  };

  points.emplace_back(Point2d(0.0, 0.0));
  points.emplace_back(Point2d(0.5, 0.0));
  points.emplace_back(Point2d(1.0, 0.0));
  points.emplace_back(Point2d(1.0, 0.5));
  points.emplace_back(Point2d(1.0, 1.0));
  points.emplace_back(Point2d(0.5, 1.0));
  points.emplace_back(Point2d(0.0, 1.0));
  points.emplace_back(Point2d(0.0, 0.5));
  for (size_t i = 0; i < points.size(); i++) {
    points[i] /= 4;
  }

  std::vector<Point2d> testPoints;
  testPoints.emplace_back(Point2d({0.5, 0.5}));
  testPoints.emplace_back(Point2d({0.4, 0.0}));
  testPoints.emplace_back(Point2d({0.0, 0.1}));
  testPoints.emplace_back(Point2d({0.45, 0.12}));
  testPoints.emplace_back(Point2d({0.35, 0.97}));
  testPoints.emplace_back(Point2d({0.78, 0.81}));
  testPoints.emplace_back(Point2d({0.78, 0.14}));
  for (size_t i = 0; i < testPoints.size(); i++) {
    testPoints[i] /= 4;
  }

  for (auto p : testPoints) {
    LinearVector wx = rbf.gradientXAt(p, points);
    LinearVector wy = rbf.gradientYAt(p, points);
    // std::cout << "-------" << '\n';
    // std::cout << wx;
    // std::cout << wy;
    // std::cout << p;
    double gx = 0;
    double gy = 0;
    for (size_t i = 0; i < points.size(); i++) {
      //   std::cout << wx[i] << " * " << _f(points[i]) << '\n';
      gx += wx[i] * _f(points[i]);
      gy += wy[i] * _f(points[i]);
    }
    // std::cout << '\n';
    // std::cout << "Gx: " << gx << "\tGy:" << gy << '\n';
    // std::cout << "Gx: " << _gx(p) << "\tGy:" << _gy(p) << '\n';
    ASSERT_NEAR(gx, _gx(p), 1e-2);
    ASSERT_NEAR(gy, _gy(p), 1e-2);
    // std::cout << "OK" << '\n';
  }
}

TEST(DifferentialRBF, gradientAt) {
  std::function<double(Point2d)> _f = [](Point2d p) {
    return std::exp(p[0]) + std::exp(p[1]); // x-only
  };
  std::function<double(Point2d)> _gx = [](Point2d p) {
    return std::exp(p[0]); // x-only
  };
  std::function<double(Point2d)> _gy = [](Point2d p) {
    return std::exp(p[1]); // x-only
  };
  DifferentialRBF<2> rbf;
  //  rbf.setBasePolynomial(Definitions::PolynomialType::CONSTANT);
  //  rbf.setKernel(new GaussianKernel<Point2d>(1.));
  double hx = 0.005;
  Point2d center(0.5, 0.5);
  double gradX = 0., gradY = 0.;
  { // Y
    std::vector<Point2d> points;
    points.push_back({center.x(), center.y() - hx});
    points.push_back({center.x(), center.y() + hx});
    std::vector<double> f;
    for (auto p : points)
      f.emplace_back(_f(p));
    gradY = rbf.gradientAt(center, 1, points, f);
    EXPECT_NEAR(gradY, _gy(center), hx * hx);
  }
  { // X
    std::vector<Point2d> points;
    points.push_back({center.x() - hx, center.y()});
    points.push_back({center.x() + hx, center.y()});
    std::vector<double> f;
    for (auto p : points)
      f.push_back(_f(p));
    gradX = rbf.gradientAt(center, 0, points, f);
    EXPECT_NEAR(gradX, _gx(center), hx * hx);
  }
  { // x-----x----o--------x---x
    std::vector<Point2d> points;
    points.push_back({0.5 - 0.5 * hx, 0.5});
    points.push_back({0.5 - 2.1 * hx, 0.5});
    points.push_back({0.5 - 4 * hx, 0.5});
    points.push_back({0.5 - 5.3 * hx, 0.5});
    points.push_back({0.5, 0.5});
    points.push_back({0.5 + 0.8 * hx, 0.5});
    points.push_back({0.5 + 4 * hx, 0.5});
    points.push_back({0.5 + 4.2 * hx, 0.5});
    points.push_back({0.5 + 4.9 * hx, 0.5});
    std::vector<double> f;
    for (auto p : points)
      f.push_back(_f(p));
    gradX = rbf.gradientAt(Point2d(0.5 + 2 * hx, 0.5), 0, points, f);
    EXPECT_NEAR(gradX, _gx(Point2d(0.5 + 2 * hx, 0.5)), hx * hx);
  }
  {
    std::vector<Point2d> points;
    points.push_back({0.5, 0.5});
    points.push_back({0.5 + 3 * hx / 2., 0.5 + hx / 2});
    std::vector<double> f;
    for (auto p : points)
      f.push_back(_f(p));
    gradX = rbf.gradientAt(Point2d(0.5 + hx, 0.5 + hx / 2), 0, points, f);
    EXPECT_NEAR(gradX, _gx(Point2d(0.5 + hx, 0.5 + hx / 2)), hx * hx);
  }
  {
    std::vector<Point2d> points;
    points.push_back({0.5, 0.5});
    points.push_back({0.5 + 3 * hx / 2., 0.5});
    std::vector<double> f;
    for (auto p : points)
      f.push_back(_f(p));
    gradX = rbf.gradientAt(Point2d(0.5 + hx, 0.5), 0, points, f);
    EXPECT_NEAR(gradX, _gx(Point2d(0.5 + hx, 0.5)), hx * hx);
  }
}

TEST(RBF, DivergenceFromPoints) {
  {
    std::vector<Point2d> points;
    points.push_back({0.2, 0.7});
    points.push_back({0.8, 0.8});
    points.push_back({0.9, 0.1});
    points.push_back({0.2, 0.2});
    std::vector<double> values(4, 1.);
    DifferentialRBF<2> rbf;
    rbf.setKernel(new CubicKernel<Point2d>());
//    rbf.gradientAt(Point2d(0.5, 0.5), 0, points, values);
  }
  {
    std::vector<Point2d> points;
    points.push_back({0.5, 0.5});
    points.push_back({0.2, 0.7});
    points.push_back({0.8, 0.8});
    points.push_back({0.9, 0.1});
    points.push_back({0.2, 0.2});
    std::vector<double> values(4, 1.);
    DifferentialRBF<2> rbf;
    rbf.setKernel(new CubicKernel<Point2d>());
    rbf.laplacianWeights(points);
  }
  return;
  {
    double h = 0.5;
    BBox2d region(Point2d(), Point2d(h, h));
    BoxSampler sampler;
    auto falloff = [](double x) { return -x * x + 2 * x; };
    //    auto falloff = [](double x) { return x; };
    auto fx = [](Point2d p) { return 0.; };
    auto fy = [](Point2d p) { return -0.981; };
    auto center = region.center();
    size_t n = 10;
    std::vector<Point2d> points;
    std::vector<double> vx, vy;
    for (size_t i = 0; i < n; i++) {
      auto point = sampler.sample(region);
      points.emplace_back(point);
      vx.emplace_back(fx(point));
      vy.emplace_back(falloff(point.y() / h) * fy(point));
    }
    DifferentialRBF<2> rbf;
    rbf.setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    std::cerr << rbf.gradientAt(center, 0, points, vx) +
                     rbf.gradientAt(center, 1, points, vy)
              << std::endl;
    RBFInterpolant<Point2d> interpolatorX, interpolatorY;
    interpolatorX.setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    interpolatorY.setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    interpolatorX.computeWeights(center, points, vx);
    interpolatorY.computeWeights(center, points, vy);
    auto kernelX = interpolatorX.kernel();
    auto kernelY = interpolatorX.kernel();
    std::vector<double> valuesX, valuesY;
    for (auto point : points) {
      valuesX.emplace_back(kernelX->dphi(point, center) * (point - center).x());
      valuesY.emplace_back(kernelY->dphi(point, center) * (point - center).y());
    }
    double valueX = interpolatorX.evaluate(center, valuesX);
    double valueY = interpolatorY.evaluate(center, valuesY);
    std::cerr << valueX + valueY << std::endl;
  }
  {
    std::vector<Point2d> points;
    std::vector<double> vx, vy;
    points.emplace_back(0.0781345, 0.0220575);
    points.emplace_back(0.0921603, 0.0290434);
    points.emplace_back(0.0851462, 0.0290452);
    points.emplace_back(0.0781332, 0.0290471);
    points.emplace_back(0.0921622, 0.0220543);
    points.emplace_back(0.0851481, 0.0220562);
    vx.emplace_back(0.0268937);
    vx.emplace_back(0.0320649);
    vx.emplace_back(0.0292424);
    vx.emplace_back(0.0266301);
    vx.emplace_back(0.0324446);
    vx.emplace_back(0.0296221);
    vy.emplace_back(-0.0375422);
    vy.emplace_back(-0.0403764);
    vy.emplace_back(-0.0400018);
    vy.emplace_back(-0.0396287);
    vy.emplace_back(-0.038189);
    vy.emplace_back(-0.0378144);
    auto center = Point<double, 2>(0.0859375, 0.0234375);
    DifferentialRBF<2> rbf;
    rbf.setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    std::cerr << rbf.gradientAt(center, 0, points, vx) +
                     rbf.gradientAt(center, 1, points, vy)
              << std::endl;
  }
}