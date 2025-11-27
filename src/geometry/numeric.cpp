#include <geometry/numeric.h>
#include <geometry/vector.h>

namespace furoo {

std::vector<double> linspace(double a, double b, int n) {
  double step = (b - a) / static_cast<double>(n - 1);
  std::vector<double> r;
  for (int i = 0; i < n; i++)
    r.emplace_back(a + i * step);
  return r;
}

double normalize(double n) {
  UNUSED_VARIABLE(n);
  return 1.0;
}

double max(double a, double b) { return std::max(a, b); }

double min(double a, double b) { return std::min(a, b); }

int ceil2Int(double f) { return static_cast<int>(f + 1.0); }

int floor2Int(double f) {
  return (f >= 0.) ? static_cast<int>(f) : static_cast<int>(f - 1.);
}

int round2Int(double f) { return static_cast<int>(f + .5); }

double mod(int a, int b) {
  int n = a / b;
  a -= n * b;
  if (a < 0)
    a += b;
  return a;
}

double linearStep(double v, double a, double b) {
  return clamp((v - a) / (b - a), 0., 1.);
}

double smoothStep(double v, double a, double b) {
  double t = clamp((v - a) / (b - a), 0., 1.);
  return t * t * (3. - 2. * t);
}

int log2Int(double x) { return floor2Int(std::log2(x)); }

bool isPowerOf2(int v) { return (v & (v - 1)) == 0; }

unsigned int roundUpPow2(unsigned int v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return v + 1;
}

bool solve_quadratic(double A, double B, double C, double &t0, double &t1) {
  double delta = B * B - 4.0 * A * C;
  if (IS_ZERO(A) || delta < 0.)
    return false;
  double sDelta = sqrt(delta);
  double q;
  if (B < 0.)
    q = -0.5 * (B - sDelta);
  else
    q = -0.5f * (B + sDelta);
  t0 = q / A;
  t1 = C / q;
  if (t0 > t1)
    std::swap(t0, t1);
  return true;
}

double smooth(double a, double b) { return std::max(0.0, 1.0 - a / SQR(b)); }

double sharpen(const double &r2, const double &h) {
  return std::max(h * h / std::max(r2, static_cast<double>(1.0e-5)) - 1.0, 0.0);
}

unsigned int separateBy1(unsigned int n) {
  n = (n ^ (n << 8)) & 0x00ff00ff;
  n = (n ^ (n << 4)) & 0x0f0f0f0f;
  n = (n ^ (n << 2)) & 0x33333333;
  n = (n ^ (n << 1)) & 0x55555555;
  return n;
}

unsigned int mortonCode(unsigned int x, unsigned int y) {
  return (separateBy1(y) << 1) + separateBy1(x);
}

unsigned int separateBy2(unsigned int n) {
  n = (n ^ (n << 16)) & 0xff0000ff;
  n = (n ^ (n << 8)) & 0x0300f00f;
  n = (n ^ (n << 4)) & 0x030c30c3;
  n = (n ^ (n << 2)) & 0x09249249;
  return n;
}

unsigned int mortonCode(unsigned int x, unsigned int y, unsigned int z) {
  return (separateBy2(z) << 2) + (separateBy2(y) << 1) + separateBy2(x);
}

double catmullRomSpline(double f0, double f1, double f2, double f3, double t) {
  double d1 = (f2 - f0) / 2;
  double d2 = (f3 - f1) / 2;
  double D1 = f2 - f1;

  if (std::fabs(D1) < 1e-8) {
    d1 = d2 = 0;
  }

  if (sign(D1) != sign(d1)) {
    d1 = 0;
  }

  if (sign(D1) != sign(d2)) {
    d2 = 0;
  }

  double a3 = d1 + d2 - 2 * D1;
  double a2 = 3 * D1 - 2 * d1 - d2;
  double a1 = d1;
  double a0 = f1;

  return a3 * CUBE(t) + a2 * SQR(t) + a1 * t + a0;
}
} // namespace furoo
