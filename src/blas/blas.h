#ifndef FUROO_BLAS_BLAS_H
#define FUROO_BLAS_BLAS_H
#include <blas/matrix_interface.h>

namespace furoo {

struct Blas {
  static void set(const VectorInterface *b, VectorInterface *a);
  /** Performs matrix-vector multiplication
   * \param m **[in]** matrix
   * \param v **[in]** vector
   * \param r **[out]** result of m * v
   */
  static void mvm(const MatrixConstAccessor &m, const VectorInterface *v,
                  VectorInterface *r);
  /** Calculates a*x + y.
   * \param a **[in]**
   * \param x **[in]**
   * \param y **[in]**
   * \param r **[out]** result of ax + y
   */
  static void axpy(double a, const VectorInterface *x, VectorInterface *y,
                   VectorInterface *r);
  /** Computes b - Ax
   * \param A **[in]**
   * \param x **[in]**
   * \param b **[in]**
   * \param r **[out]** result of b - Ax
   */
  static void residual(const MatrixConstAccessor &A, const VectorInterface *x,
                       const VectorInterface *b, VectorInterface *r);
  static double norm(const VectorInterface *v);
  static double infnorm(const VectorInterface *v);
  static void dSolve(const MatrixConstAccessor &A, const VectorInterface *b,
                     VectorInterface *x);
  static double dot(const VectorInterface *a, const VectorInterface *b);
};

} // furoo namespace
#endif
