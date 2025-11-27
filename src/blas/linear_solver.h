#ifndef FUROO_BLAS_LINEAR_SOLVER_H
#define FUROO_BLAS_LINEAR_SOLVER_H

#include <blas/matrix_interface.h>
#include <cmath>
#include <vector>

namespace furoo {

class LinearSolverInterface {
public:
  LinearSolverInterface() {}
  virtual bool solve(MatrixInterface *, VectorInterface *,
                     VectorInterface *) = 0;
};

class IterativeLinearSolver : public LinearSolverInterface {
protected:
  double tolerance;
  size_t maxNumberOfIterations;

public:
  IterativeLinearSolver();
  virtual bool solve(MatrixInterface *, VectorInterface *,
                     VectorInterface *) = 0;
  void setTolerance(double);
  void setMaxNumberOfIterations(size_t);
};

class GaussJordanLinearSolver : public LinearSolverInterface {
public:
  bool solve(MatrixInterface *, VectorInterface *, VectorInterface *) override;
};

class BiconjugateGradientLinearSolver : public IterativeLinearSolver {
private:
  bool _bicg(const MatrixConstAccessor &A, const MatrixConstAccessor &At,
             VectorInterface *x, const VectorInterface *b, VectorInterface *r,
             VectorInterface *rt, VectorInterface *p, VectorInterface *pt,
             VectorInterface *z, VectorInterface *zt);

public:
  bool solve(MatrixInterface *, VectorInterface *, VectorInterface *) override;
};

class ConjugateGradientLinearSolver : public IterativeLinearSolver {
public:
  bool solve(MatrixInterface *, VectorInterface *, VectorInterface *) override;

private:
  bool _cg(const MatrixConstAccessor &A, VectorInterface *x,
           const VectorInterface *b, VectorInterface *r, VectorInterface *q,
           VectorInterface *d);
};

} // furoo namespace

#endif
