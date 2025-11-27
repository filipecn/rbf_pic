#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(TemporalIntegrator, EulerMethod) {
  EulerTemporalIntegrator<2> eti;
  Point2d p(0.5, 0.0);
  double dt = 0.2;
  std::function<double(double)> f = [](double t) {
    return SQR(t + 1.) - .5 * std::exp(t);
  };

  for (size_t t = 0; t < 10; t++) {
    Vector2d v(p.x() - SQR(t * dt) + 1, 0.);
    auto newP = eti.integrate(p, v, dt);
    ASSERT_NEAR(newP.x(), f(t * dt), 5e-1);
    p = newP;
  }
}

TEST(TemporalIntegrator, RungeKutta2) {
  return;
  RungeKutta2TemporalIntegrator<2> eti;
  Point2d p(0.5, 0.0);
  double dt = 0.2;
  std::function<double(double)> f = [](double t) {
    return SQR(t + 1.) - .5 * std::exp(t);
  };

  for (size_t t = 0; t < 10; t++) {
    std::cout << p << '\n';
    Vector2d v(p.x() - SQR(t * dt) + 1, 0.);
    auto newP = eti.integrate(p, v, dt);
    ASSERT_NEAR(newP.x(), f(t * dt), 5e-1);
    p = newP;
  }
}
