#ifndef FUROO_GEOMETRY_NUMERIC_H
#define FUROO_GEOMETRY_NUMERIC_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace furoo {

#ifndef INFINITY
#define INFINITY std::numeric_limits<double>::max()
#endif

#ifndef INT_INFINITY
#define INT_INFINITY std::numeric_limits<int>::max()
#endif

#define PI 3.14159265358979323846
#define PI_2 6.28318530718
#define INV_PI 0.31830988618379067154
#define INV_TWOPI 0.15915494309189533577
#define INV_FOURPI 0.07957747154594766788

#define SQR(A) ((A) * (A))

#define CUBE(A) ((A) * (A) * (A))

#define TO_DEGREES(A) ((A)*180.0 / PI)

#define TO_RADIANS(A) ((A)*PI / 180.0)

#define IS_ZERO(A) (std::fabs(A) < 1e-8)

#define IS_NOT_ZERO(A) (std::fabs(A) > 1e-8)

#define IS_EQUAL(A, B) (std::fabs((A) - (B)) < 1e-8)

#define IS_EQUAL_ERROR(A, B, C) (std::fabs((A) - (B)) < C)

#define IS_BETWEEN(A, B, C) ((A) > (B) && (A) < (C))

#define IS_BETWEEN_CLOSE(A, B, C) ((A) >= (B) && (A) <= (C))

template<typename T>
T sign(T x) {
  return x >= static_cast<T>(0.) ? 1 : -1;
}

std::vector<double> linspace(double a, double b, int n);
/** \brief computes the arc-tangent of y/x
 * \param y
 * \param x
 * \returns angle (in radians) between **0** and **2PI**
 */
template<typename T>
T atanPI_2(T y, T x);

/** \brief normalization (dummy function)
 * \param n
 * \returns n
 */
double normalize(double n);

/** \brief max function
 * \param a
 * \param b
 * \returns a if a >= b, b otherwise
 */
double max(double a, double b);

/** \brief min function
 * \param a
 * \param b
 * \returns a if a <= b, b otherwise
 */
double min(double a, double b);

/** \brief round
 * \param f **[in]**
 * \return ceil of **f**
 */
int ceil2Int(double f);

/** \brief round
 * \param f **[in]**
 * \return floor of **f**
 */
int floor2Int(double f);

/** \brief round
 * \param f **[in]**
 * \return next integer greater or equal to **f**
 */
int round2Int(double f);

/** \brief modulus
 * \param a **[in]**
 * \param b **[in]**
 * \return the remainder of a / b
 */
double mod(int a, int b);

template<typename T, typename S>
T bisect(T a, T b, S fa, S fb, S error, const std::function<S(T p)> &f);

/** \brief interpolation
 * \param t **[in]** parametric coordinate
 * \param a **[in]** lower bound **0**
 * \param b **[in]** upper bound **1**
 * \return linear interpolation between **a** and **b** at **t**.
 */
template<typename T, typename S>
S lerp(T t, const S &a, const S &b);

/** \brief bilinear interpolate
 * \param x **[in]** parametric coordinate in x
 * \param y **[in]** parametric coordinate in y
 * \param f00 **[in]** function value at **(0, 0)**
 * \param f10 **[in]** function value at **(1, 0)**
 * \param f11 **[in]** function value at **(1, 1)**
 * \param f01 **[in]** function value at **(0, 1)**
 * \return interpolated value at **(x,y)**
 */
template<typename T = double, typename S = double>
S bilerp(T x, T y, const S &f00, const S &f10, const S &f11, const S &f01);

template<typename T = double, typename S = double>
S trilerp(T tx, T ty, T tz, const S &f000, const S &f100, const S &f010,
          const S &f110, const S &f001, const S &f101, const S &f011,
          const S &f111);

double catmullRomSpline(double f0, double f1, double f2, double f3, double t);

/** \brief interpolates to the nearest value
 * \param t **[in]** parametric coordinate
 * \param a **[in]** lower bound
 * \param b **[in]** upper bound
 * \return interpolation between **a** and **b** at **t**.
 */
template<typename T = double, typename S = double>
S nearest(T t, const S &a, const S &b);

/** \brief clamp
 * \param n **[in]** value
 * \param l **[in]** low
 * \param u **[in]** high
 * \return clame **b** to be in **[l, h]**
 */
template<typename T>
T clamp(const T &n, const T &l, const T &u);

/** \brief smooth Hermit interpolation when **a** < **v** < **b**
 * \param v **[in]** coordinate
 * \param a **[in]** lower bound
 * \param b **[in]** upper bound
 * \return Hermit value between **0** and **1**
 */
double smoothStep(double v, double a, double b);

/** \brief linear interpolation
 * \param v **[in]** coordinate
 * \param a **[in]** lower bound
 * \param b **[in]** upper bound
 * \return linear value between **0** and **1**
 */
double linearStep(double v, double a, double b);

/** \brief log
 * \param x **[in]** value
 * \return integer base-2 logarithm of **x**
 */
int log2Int(double x);

/** \brief query
 * \param v **[in]** value
 * \return **true** if **v** is power of 2
 */
bool isPowerOf2(int v);

/** \brief round
 * \param v **[in]** value
 * \return the next number in power of 2
 */
unsigned int roundUpPow2(unsigned int v);
bool solve_quadratic(double A, double B, double C, double &t0, double &t1);
double trilinear_hat_function(double r);
double quadraticBSpline(double r);

template<typename T>
T bilinearInterpolation(T f00, T f10, T f11, T f01, T x, T y);

template<typename T>
T cubicInterpolate(T p[4], T x);

template<typename T>
T bicubicInterpolate(T p[4][4], T x, T y);

template<typename T>
T trilinearInterpolate(double *p, T ***data, T b, const int dimensions[3]);

template<typename T>
T tricubicInterpolate(double *p, T ***data);

template<typename T>
T tricubicInterpolate(double *p, T ***data, T b, const int dimensions[3]);
double smooth(double a, double b);
double sharpen(const double &r2, const double &h);

/** \brief split
 * \param n **[in]** number
 * \return the number with the bits of **n** splitted and separated by 1 bit.
 * Ex: n = 1101, return 101001
 */
unsigned int separateBy1(unsigned int n);

/** \brief morton code
 * \param x **[in]** number
 * \param y **[in]** number
 * \return the morton code of **x** and **y**. The morton code is the
 * combination of the two binary numbers with the bits interleaved.
 */
unsigned int mortonCode(unsigned int x, unsigned int y);

unsigned int separateBy2(unsigned int n);

unsigned int mortonCode(unsigned int x, unsigned int y, unsigned int z);

//!
//! \brief      Returns the square of \p x.
//!
//! \param[in]  x     The input.
//!
//! \tparam     T     Value type.
//!
//! \return     The squared value.
//!
template<typename T>
inline T square(T x);

//!
//! \brief      Returns the cubic of \p x.
//!
//! \param[in]  x     The input.
//!
//! \tparam     T     Value type.
//!
//! \return     The cubic of \p x.
//!
template<typename T>
inline T cubic(T x);

#include "geometry/numeric.inl"
} // namespace furoo

#endif
