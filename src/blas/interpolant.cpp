#include <blas/interpolation.h>

namespace furoo {

double BilinearInterpolant::interpolateAt(Point2d target,
                                          const std::vector<Point2d> &points,
                                          const std::vector<double> &f) {
  // ASSERT_FATAL(points.size() == 4);
  LinearMatrix A(4);
  LinearVector b(4);
  GaussJordanLinearSolver solver;
  for (size_t i = 0; i < 4; i++) {
    Point2d p = points.at(i);
    A(i, 0) = 1;
    A(i, 1) = p[0];
    A(i, 2) = p[1];
    A(i, 3) = p[0] * p[1];
    b[i] = f[i];
  }
  solver.solve(&A, &b, &b);

  return b[0] + b[1] * target[0] + b[2] * target[1] +
         b[3] * target[0] * target[1];
}

void BilinearInterpolant::weights(const std::vector<Point2d> &points,
                                  const std::vector<double> &f,
                                  std::vector<double> &w) {
  UNUSED_VARIABLE(points);
  NOT_IMPLEMENTED();
  // w[0] = (1 - target.x()) * (1 - target.y());
  // w[1] = target.x() * (1 - target.y());
  // w[2] = target.x() * target.y();
  // w[3] = (1 - target.x()) * target.y();
}

BicubicInterpolant::BicubicInterpolant() {
  static std::vector<double> coef = {
      1,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  1,  0, 0,  0,  0,  0,  0,  0, -3, 0,  0,  3,  0,  0,
      0,  0,  -2, 0,  0,  -1, 0, 0,  0,  0,  2,  0,  0, -2, 0,  0,  0,  0,  1,
      0,  0,  1,  0,  0,  0,  0, 0,  0,  0,  0,  1,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  1,  0,  0,
      0,  0,  0,  0,  0,  -3, 0, 0,  3,  0,  0,  0,  0, -2, 0,  0,  -1, 0,  0,
      0,  0,  2,  0,  0,  -2, 0, 0,  0,  0,  1,  0,  0, 1,  -3, 3,  0,  0,  -2,
      -1, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      -3, 3,  0,  0,  -2, -1, 0, 0,  9,  -9, 9,  -9, 6, 3,  -3, -6, 6,  -6, -3,
      3,  4,  2,  1,  2,  -6, 6, -6, 6,  -4, -2, 2,  4, -3, 3,  3,  -3, -2, -1,
      -1, -2, 2,  -2, 0,  0,  1, 1,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0, 2,  -2, 0,  0,  1,  1, 0,  0,  -6, 6,  -6, 6,
      -3, -3, 3,  3,  -4, 4,  2, -2, -2, -2, -1, -1, 4, -4, 4,  -4, 2,  2,  -2,
      -2, 2,  -2, -2, 2,  1,  1, 1,  1};
  m = LinearMatrix(16, 16, coef);
}

double BicubicInterpolant::interpolateAt(Point2d target,
                                         const std::vector<Point2d> &points,
                                         const std::vector<double> &f) {
  UNUSED_VARIABLE(target);
  UNUSED_VARIABLE(points);
  UNUSED_VARIABLE(f);
  return 0.;
}

double BicubicInterpolant::interpolateAt(
    Point2d target, const std::vector<double> &f, const std::vector<double> &fx,
    const std::vector<double> &fy, const std::vector<double> &fxy, Vector2d h,
    double *gx, double *gy) const {
  LinearVector x(16);
  double hh = h[0] * h[1];
  for (size_t i = 0; i < 4; i++) {
    x[i] = f[i];
    x[i + 4] = fx[i] * h[0];
    x[i + 8] = fy[i] * h[1];
    x[i + 12] = fxy[i] * hh;
  }
  LinearVector cl = m * x;
  LinearMatrix c(16);
  size_t l = 0;
  for (size_t i = 0; i < 4; i++)
    for (size_t j = 0; j < 4; j++)
      c(i, j) = cl[l++];
  double t = target.x();
  double u = target.y();
  double r = 0.;
  for (int i = 3; i >= 0; i--)
    r = t * r + ((c(i, 3) * u + c(i, 2)) * u + c(i, 1)) * u + c(i, 0);
  if (gx && gy) {
    for (int i = 3; i >= 0; i--) {
      *gy = t * (*gy) + (3.0 * c(i, 3) * u + 2.0 * c(i, 2)) * u + c(i, 1);
      *gx = u * (*gx) + (3.0 * c(3, i) * t + 2.0 * c(2, i)) * t + c(1, i);
    }
    *gx /= h[0];
    *gy /= h[1];
  }
  return r;
}

void BicubicInterpolant::weights(const std::vector<Point2d> &points,
                                 const std::vector<double> &f,
                                 std::vector<double> &w) {
  UNUSED_VARIABLE(points);
  UNUSED_VARIABLE(w);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Monotonic (clamped) Catmull-Rom Interpolation

void MonotonicCatmullRom::getFraction(double x, int low, int high, int *index,
                                      double *f) {
  double s = std::floor(x);
  *index = static_cast<int>(s);

  int offset = -low;
  low += offset;
  high += offset;

  if (low == high) {
    *index = low;
    *f = 0.;
  } else if (*index < low) {
    *index = low;
    *f = 0.;
  } else if (*index > high - 1) {
    *index = high - 1;
    *f = 1.;
  } else {
    *f = x - s;
  }

  *index -= offset;
}

double MonotonicCatmullRom::interpolate(double v0, double v1, double v2,
                                        double v3, double f) {

  double d1 = (v2 - v0) / 2.;
  double d2 = (v3 - v1) / 2.;
  double D1 = v2 - v1;

  if (std::abs(D1) < 0.0001)
    d1 = d2 = 0;

  if (sign(D1) != sign(d1))
    d1 = 0;

  if (sign(D1) != sign(d2))
    d2 = 0;

  double a3 = d1 + d2 - 2 * D1;
  double a2 = 3 * D1 - 2 * d1 - d2;
  double a1 = d1;
  double a0 = v1;

  return a3 * cubic(f) + a2 * square(f) + a1 * f + a0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace furoo
