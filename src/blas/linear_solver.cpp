#include <blas/blas.h>
#include <blas/linear_solver.h>
#include <blas/linear_vector.h>
#include <common/debug.h>
#include <geometry/numeric.h>

namespace furoo {

IterativeLinearSolver::IterativeLinearSolver() {
  tolerance = 1e-6;
  maxNumberOfIterations = 1000;
}

void IterativeLinearSolver::setTolerance(double tol) { tolerance = tol; }

void IterativeLinearSolver::setMaxNumberOfIterations(size_t iter) {
  maxNumberOfIterations = iter;
}

bool GaussJordanLinearSolver::solve(MatrixInterface *_A, VectorInterface *_b,
                                    VectorInterface *_x) {
  UNUSED_VARIABLE(_x);
  auto &A = *_A;
  auto &b = *_b;
  size_t n = A.rowCount(), irow = 0, icol = 0;

  std::vector<size_t> ipiv(n, 0), indxr(n), indxc(n);
  for (size_t i = 0; i < n; i++) {
    double big = 0;
    for (size_t j = 0; j < n; j++) {
      if (ipiv[j] != 1) {
        for (size_t k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (std::fabs(A(j, k)) >= big) {
              big = std::fabs(A(j, k));
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (size_t l = 0; l < n; l++) {
        std::swap(A(irow, l), A(icol, l));
      }
      std::swap(b[irow], b[icol]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (A(icol, icol) == 0.0) {
      std::cerr << "matrix size " << A.columnCount() << " " << A.rowCount() << std::endl;
      throw("gaussj: Singular Matrix");
    }
    double pivinv = 1.0 / A(icol, icol);
    A(icol, icol) = 1.0;
    for (size_t l = 0; l < n; l++) {
      A(icol, l) *= pivinv;
    }
    b[icol] *= pivinv;
    for (size_t ll = 0; ll < n; ll++) {
      if (ll != icol) {
        double dum = A(ll, icol);
        A(ll, icol) = 0.0;
        for (size_t l = 0; l < n; l++) {
          A(ll, l) -= A(icol, l) * dum;
        }
        b[ll] -= b[icol] * dum;
      }
    }
  }
  for (int l = static_cast<int>(n) - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l]) {
      for (size_t k = 0; k < n; k++) {
        std::swap(A(k, indxr[l]), A(k, indxc[l]));
      }
    }
  }
  return true;
}

bool BiconjugateGradientLinearSolver::solve(MatrixInterface *_A,
                                            VectorInterface *_b,
                                            VectorInterface *_x) {
  MatrixConstAccessor A(_A);
  MatrixConstAccessor At(_A, true);
  LinearVector r(std::vector<double>(_b->size(), 0.));
  LinearVector rt(std::vector<double>(_b->size(), 0.));
  LinearVector p(std::vector<double>(_b->size(), 0.));
  LinearVector pt(std::vector<double>(_b->size(), 0.));
  LinearVector z(std::vector<double>(_b->size(), 0.));
  LinearVector zt(std::vector<double>(_b->size(), 0.));

  return _bicg(A, At, _x, _b, &r, &rt, &p, &pt, &z, &zt);
}

bool BiconjugateGradientLinearSolver::_bicg(
    const MatrixConstAccessor &A, const MatrixConstAccessor &At,
    VectorInterface *x, const VectorInterface *b, VectorInterface *r,
    VectorInterface *rt, VectorInterface *p, VectorInterface *pt,
    VectorInterface *z, VectorInterface *zt) {
  size_t n = x->size();
  double bnrm, znrm, xnrm, zm1nrm, dxnrm, bkden = 1.0, bk, err;
  if (IS_ZERO(Blas::infnorm(b)))
    return true;
  Blas::residual(A, x, b, r);
  for (size_t j = 0; j < n; j++) {
    (*rt)[j] = (*r)[j];
  }
  int itol = 1;
  /* {
   * std::cout << "r0\n";
   * for (size_t i = 0; i < n; i++)
   *   std::cout << (*r)[i] << " ";
   * std::cout << std::endl;
   * }*/
  if (itol == 1) {
    bnrm = Blas::norm(b);
    Blas::dSolve(A, r, z);
  } else if (itol == 2) {
    Blas::dSolve(A, b, z);
    bnrm = Blas::norm(z);
    Blas::dSolve(A, r, z);
  } else {
    Blas::dSolve(A, b, z);
    bnrm = (itol == 3) ? Blas::norm(z) : Blas::infnorm(z);
    Blas::dSolve(A, r, z);
    znrm = (itol == 3) ? Blas::norm(z) : Blas::infnorm(z);
  }

  /* {
   * std::cout << "z0\n";
   * for (size_t i = 0; i < n; i++)
   *   std::cout << (*z)[i] << " ";
   * std::cout << std::endl;
   * }*/
  size_t iter = 0;
  for (; iter < maxNumberOfIterations; ++iter) {
    Blas::dSolve(At, rt, zt);
    double bknum = Blas::dot(rt, z);
    if (iter == 0) {
      for (size_t j = 0; j < n; j++) {
        (*p)[j] = (*z)[j];
        (*pt)[j] = (*zt)[j];
      }
    } else {
      bk = bknum / bkden;
      Blas::axpy(bk, p, z, p);
      Blas::axpy(bk, pt, zt, pt);
    }
    bkden = bknum;
    Blas::mvm(A, p, z);
    double akden = Blas::dot(z, pt);
    double ak = bknum / akden;
    Blas::mvm(At, pt, zt);
    Blas::axpy(ak, p, x, x);
    Blas::axpy(-ak, z, r, r);
    Blas::axpy(-ak, zt, rt, rt);
    Blas::dSolve(A, r, z);
    if (itol == 1) {
      err = Blas::norm(r) / bnrm;
    } else if (itol == 2) {
      err = Blas::norm(z) / bnrm;
    } else {
      zm1nrm = znrm;
      znrm = (itol == 3) ? Blas::norm(z) : Blas::infnorm(z);
      if (fabs(zm1nrm - znrm) > tolerance * znrm) {
        dxnrm = fabs(ak) * ((itol == 3) ? Blas::norm(p) : Blas::infnorm(p));
        err = znrm / fabs(zm1nrm - znrm) * dxnrm;
      } else {
        err = znrm / bnrm;
        continue;
      }
      xnrm = (itol == 3) ? Blas::norm(x) : Blas::infnorm(x);
      if (err <= 0.5 * xnrm) {
        err /= xnrm;
      } else {
        err = znrm / bnrm;
        continue;
      }
    }
    if (err <= tolerance)
      break;
  }
  return !(iter == maxNumberOfIterations);
}

bool ConjugateGradientLinearSolver::solve(MatrixInterface *_A,
                                          VectorInterface *_b,
                                          VectorInterface *_x) {
  MatrixConstAccessor A(_A);
  LinearVector r(std::vector<double>(_b->size(), 0.));
  LinearVector q(std::vector<double>(_b->size(), 0.));
  LinearVector d(std::vector<double>(_b->size(), 0.));

  return _cg(A, _x, _b, &r, &d, &q);
}

bool ConjugateGradientLinearSolver::_cg(const MatrixConstAccessor &A,
                                        VectorInterface *x,
                                        const VectorInterface *b,
                                        VectorInterface *r, VectorInterface *q,
                                        VectorInterface *d) {
  if (IS_ZERO(Blas::infnorm(b)))
    return true;
  // r = b - A * x
  Blas::residual(A, x, b, r);
  // d = r
  Blas::set(r, d);
  // sigma = r '* r
  double sigma = Blas::dot(r, r);
  unsigned it = 0;
  while (it < maxNumberOfIterations) {
    // q = Ad
    Blas::mvm(A, d, q);
    // alpha = sigma / (d '* Ad)
    double alpha = sigma / Blas::dot(d, q);
    // x = alpha * d + x
    Blas::axpy(alpha, d, x, x);
    // r = r - alpha * Ad
    Blas::axpy(-alpha, q, r, r);
    // sigmaNew = r '* r
    double sigmaNew = Blas::dot(r, r);
    if (sigmaNew < SQR(tolerance))
      break;
    // d = r + (sigmaNew / sigma) * d
    Blas::axpy(sigmaNew / sigma, d, r, d);
    sigma = sigmaNew;
    ++it;
  }

  unsigned lastNumberOfIterations = it;
  double lastResidual = sigma;
  return lastResidual <= tolerance ||
         lastNumberOfIterations < maxNumberOfIterations;
}

} // namespace furoo
