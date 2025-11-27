template <typename T> T atanPI_2(T y, T x) {
  T angle = std::atan2(y, x);
  if (angle < 0.0)
    return PI_2 + angle;
  return angle;
}

template <typename T>
inline T square(T x) {
  return x * x;
}

template <typename T>
inline T cubic(T x) {
  return x * x * x;
}

template <typename T, typename S>
T bisect(T a, T b, S fa, S fb, S error, const std::function<S(T p)> &f) {
  if (fb < fa)
    return bisect(b, a, fb, fa, error, f);
  S zero = 0;
  if (fa > zero || IS_EQUAL_ERROR(fa, zero, error))
    return a;
  if (fb < zero || IS_EQUAL_ERROR(fb, zero, error))
    return b;
  T c = (a + b) / static_cast<T>(2);
  if (IS_EQUAL(a, c) || IS_EQUAL(b, c))
    return c;
  S fc = f(c);
  if (fc > zero)
    return bisect(a, c, fa, fc, error, f);
  return bisect(c, b, fc, fb, error, f);
}

template <typename T, typename S> S lerp(T t, const S &a, const S &b) {
  return (1.0 - t) * a + t * b;
}

template <typename T, typename S>
S bilerp(T x, T y, const S &f00, const S &f10, const S &f11, const S &f01) {
  return lerp(y, lerp(x, f00, f10), lerp(x, f01, f11));
}

template <typename T, typename S>
S trilerp(T tx, T ty, T tz, const S &f000, const S &f100, const S &f010,
          const S &f110, const S &f001, const S &f101, const S &f011,
          const S &f111) {
  return lerp(bilerp(f000, f100, f010, f110, tx, ty),
              bilerp(f001, f101, f011, f111, tx, ty), tz);
}

template <typename S, typename T>
S catmullRomSpline(const S &f0, const S &f1, const S &f2, const S &f3, T f) {
  S d1 = (f2 - f0) / 2;
  S d2 = (f3 - f1) / 2;
  S D1 = f2 - f1;

  S a3 = d1 + d2 - 2 * D1;
  S a2 = 3 * D1 - 2 * d1 - d2;
  S a1 = d1;
  S a0 = f1;

  return a3 * CUBE(f) + a2 * SQR(f) + a1 * f + a0;
}

template <typename T, typename S> S nearest(T t, const S &a, const S &b) {
  return (t < static_cast<T>(0.5)) ? a : b;
}

template <typename T> T clamp(const T &n, const T &l, const T &u) {
  return std::max(l, std::min(n, u));
}

template <typename T>
T bilinearInterpolation(T f00, T f10, T f11, T f01, T x, T y) {
  return f00 * (1.0 - x) * (1.0 - y) + f10 * x * (1.0 - y) +
         f01 * (1.0 - x) * y + f11 * x * y;
}
template <typename T> T cubicInterpolate(T p[4], T x) {
  return p[1] +
         0.5 * x * (p[2] - p[0] +
                    x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                         x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}
template <typename T> T bicubicInterpolate(T p[4][4], T x, T y) {
  T arr[4];
  arr[0] = cubicInterpolate(p[0], y);
  arr[1] = cubicInterpolate(p[1], y);
  arr[2] = cubicInterpolate(p[2], y);
  arr[3] = cubicInterpolate(p[3], y);
  return cubicInterpolate(arr, x);
}
template <typename T>
T trilinearInterpolate(double *p, T ***data, T b, const int dimensions[3]) {
  int i0 = p[0], j0 = p[1], k0 = p[2];
  int i1 = p[0] + 1, j1 = p[1] + 1, k1 = p[2] + 1;
  double x = p[0] - i0;
  double y = p[1] - j0;
  double z = p[2] - k0;
  T v000 = (i0 < 0 || j0 < 0 || k0 < 0 || i0 >= dimensions[0] ||
            j0 >= dimensions[1] || k0 >= dimensions[2])
               ? b
               : data[i0][j0][k0];
  T v001 = (i0 < 0 || j0 < 0 || k1 < 0 || i0 >= dimensions[0] ||
            j0 >= dimensions[1] || k1 >= dimensions[2])
               ? b
               : data[i0][j0][k1];
  T v010 = (i0 < 0 || j1 < 0 || k0 < 0 || i0 >= dimensions[0] ||
            j1 >= dimensions[1] || k0 >= dimensions[2])
               ? b
               : data[i0][j1][k0];
  T v011 = (i0 < 0 || j1 < 0 || k1 < 0 || i0 >= dimensions[0] ||
            j1 >= dimensions[1] || k1 >= dimensions[2])
               ? b
               : data[i0][j1][k1];
  T v100 = (i1 < 0 || j0 < 0 || k0 < 0 || i1 >= dimensions[0] ||
            j0 >= dimensions[1] || k0 >= dimensions[2])
               ? b
               : data[i1][j0][k0];
  T v101 = (i1 < 0 || j0 < 0 || k1 < 0 || i1 >= dimensions[0] ||
            j0 >= dimensions[1] || k1 >= dimensions[2])
               ? b
               : data[i1][j0][k1];
  T v110 = (i1 < 0 || j1 < 0 || k0 < 0 || i1 >= dimensions[0] ||
            j1 >= dimensions[1] || k0 >= dimensions[2])
               ? b
               : data[i1][j1][k0];
  T v111 = (i1 < 0 || j1 < 0 || k1 < 0 || i1 >= dimensions[0] ||
            j1 >= dimensions[1] || k1 >= dimensions[2])
               ? b
               : data[i1][j1][k1];
  return v000 * (1.0 - x) * (1.0 - y) * (1.0 - z) +
         v100 * x * (1.0 - y) * (1.0 - z) + v010 * (1.0 - x) * y * (1.0 - z) +
         v110 * x * y * (1.0 - z) + v001 * (1.0 - x) * (1.0 - y) * z +
         v101 * x * (1.0 - y) * z + v011 * (1.0 - x) * y * z + v111 * x * y * z;
}
template <typename T> T tricubicInterpolate(double *p, T ***data) {
  int x, y, z;
  int i, j, k;
  double dx, dy, dz;
  double u[4], v[4], w[4];
  T r[4], q[4];
  T vox = T(0);

  x = (int)p[0], y = (int)p[1], z = (int)p[2];
  dx = p[0] - (double)x, dy = p[1] - (double)y, dz = p[2] - (double)z;

  u[0] = -0.5 * CUBE(dx) + SQR(dx) - 0.5 * dx;
  u[1] = 1.5 * CUBE(dx) - 2.5 * SQR(dx) + 1;
  u[2] = -1.5 * CUBE(dx) + 2 * SQR(dx) + 0.5 * dx;
  u[3] = 0.5 * CUBE(dx) - 0.5 * SQR(dx);

  v[0] = -0.5 * CUBE(dy) + SQR(dy) - 0.5 * dy;
  v[1] = 1.5 * CUBE(dy) - 2.5 * SQR(dy) + 1;
  v[2] = -1.5 * CUBE(dy) + 2 * SQR(dy) + 0.5 * dy;
  v[3] = 0.5 * CUBE(dy) - 0.5 * SQR(dy);

  w[0] = -0.5 * CUBE(dz) + SQR(dz) - 0.5 * dz;
  w[1] = 1.5 * CUBE(dz) - 2.5 * SQR(dz) + 1;
  w[2] = -1.5 * CUBE(dz) + 2 * SQR(dz) + 0.5 * dz;
  w[3] = 0.5 * CUBE(dz) - 0.5 * SQR(dz);

  int ijk[3] = {x - 1, y - 1, z - 1};
  for (k = 0; k < 4; k++) {
    q[k] = 0;
    for (j = 0; j < 4; j++) {
      r[j] = 0;
      for (i = 0; i < 4; i++) {
        r[j] += u[i] * data[ijk[0]][ijk[1]][ijk[2]];
        ijk[0]++;
      }
      q[k] += v[j] * r[j];
      ijk[0] = x - 1;
      ijk[1]++;
    }
    vox += w[k] * q[k];
    ijk[0] = x - 1;
    ijk[1] = y - 1;
    ijk[2]++;
  }
  return (vox < T(0) ? T(0.0) : vox);
}
template <typename T>
T tricubicInterpolate(double *p, T ***data, T b, const int dimensions[3]) {
  int x, y, z;
  int i, j, k;
  double dx, dy, dz;
  double u[4], v[4], w[4];
  T r[4], q[4];
  T vox = T(0);

  x = (int)p[0], y = (int)p[1], z = (int)p[2];
  dx = p[0] - (double)x, dy = p[1] - (double)y, dz = p[2] - (double)z;

  u[0] = -0.5 * CUBE(dx) + SQR(dx) - 0.5 * dx;
  u[1] = 1.5 * CUBE(dx) - 2.5 * SQR(dx) + 1;
  u[2] = -1.5 * CUBE(dx) + 2 * SQR(dx) + 0.5 * dx;
  u[3] = 0.5 * CUBE(dx) - 0.5 * SQR(dx);

  v[0] = -0.5 * CUBE(dy) + SQR(dy) - 0.5 * dy;
  v[1] = 1.5 * CUBE(dy) - 2.5 * SQR(dy) + 1;
  v[2] = -1.5 * CUBE(dy) + 2 * SQR(dy) + 0.5 * dy;
  v[3] = 0.5 * CUBE(dy) - 0.5 * SQR(dy);

  w[0] = -0.5 * CUBE(dz) + SQR(dz) - 0.5 * dz;
  w[1] = 1.5 * CUBE(dz) - 2.5 * SQR(dz) + 1;
  w[2] = -1.5 * CUBE(dz) + 2 * SQR(dz) + 0.5 * dz;
  w[3] = 0.5 * CUBE(dz) - 0.5 * SQR(dz);

  int ijk[3] = {x - 1, y - 1, z - 1};
  for (k = 0; k < 4; k++) {
    q[k] = 0;
    for (j = 0; j < 4; j++) {
      r[j] = 0;
      for (i = 0; i < 4; i++) {
        if (ijk[0] < 0 || ijk[0] >= dimensions[0] || ijk[1] < 0 ||
            ijk[1] >= dimensions[1] || ijk[2] < 0 || ijk[2] >= dimensions[2])
          r[j] += u[i] * b;
        else
          r[j] += u[i] * data[ijk[0]][ijk[1]][ijk[2]];
        ijk[0]++;
      }
      q[k] += v[j] * r[j];
      ijk[0] = x - 1;
      ijk[1]++;
    }
    vox += w[k] * q[k];
    ijk[0] = x - 1;
    ijk[1] = y - 1;
    ijk[2]++;
  }
  return (vox < T(0) ? T(0.0) : vox);
}
