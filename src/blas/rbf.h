#ifndef FUROO_BLAS_RBF_H
#define FUROO_BLAS_RBF_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <blas/linear_matrix.h>
#include <blas/linear_solver.h>
#include <blas/linear_vector.h>
#include <blas/rbf_kernels.h>
#include <common/code_profiler.h>
#include <common/definitions.h>
#include <geometry/point.h>
#include <memory>

namespace furoo {

/// Performs RBF based differential operations
/// \tparam D number of dimensions
template <int D> class DifferentialRBF {
public:
  /// DifferentialRBF class constructor.
  /// If a kernel is not provided, the default kernel is assumed to be
  /// Polyharmonic quintic kernel. It can be changed later, using the
  /// setKernel(.) method provided in this class
  /// \param kernel kernel to be used when computing differential values.
  /// QuinticKernel is the default kernel
  explicit DifferentialRBF(RadialKernel<Point<double, D>> &kernel =
    *defaultKernel());

  /// Change the current radial kernel used to compute differential weights.
  /// \param kernel pointer to the new kernel
  void setKernel(RadialKernel<Point<double, D>> &kernel);

  RadialKernel<Point<double, D>> *kernel();

  /// Change the polynomial base augmentation used in rbf to compute
  /// differential weights.
  /// Currently supported polynomials: LINEAR, QUADRATIC
  /// \param type new polynomial type to be used
  void setBasePolynomial(Definitions::PolynomialType type);

  ///
  /// \param target query position
  /// \param points list of stencil points
  void setStencil();

  /// Changes current target point of the stencil. If changed, the stencil
  /// should be recomputed as well
  /// \param target new target point
  void setTarget(const Point<double, D> target);

  /// Adds a new point to the stencil. The target point, where the differential
  /// values will be computed has to be added to the stencil as well
  /// \param point new point to be added to the computation
  void addStencilPoint(const Point<double, D> point);

  /// Compute divergent value at target position, given the stencil coordinates
  /// and their current values. Return the differential value properly.
  /// \param target query position
  /// \param points list of positions of sample points
  /// \param f list of function values on sample points \return divergent at
  /// target
  double divergentAt(Point<double, D> target,
                     const std::vector<Point<double, D>> &points,
                     const std::vector<Vector<double, D>> &f);

  /// Compute gradient value at target position and dimension coordinate using
  /// a pre-computed stencil matrix (via setStencil).
  /// NOTE: the target is supposed to belong to the stencil
  /// \param dimension  0 = X, 1 = Y, 2 = Z
  /// \param f list of function values on sample points
  /// \return gradient component at target on direction <dimension>
  double gradientAt(int dimension, const std::vector<double> &f);

  /// Compute gradient value at target position and dimension coordinate. Since
  /// gradient is supposed to be a vector, the dimension parameter is necessary
  /// so the method can return the expected value. The stencil is supposed to be
  /// given at points, and their current values in f. If target is part of this
  /// stencil, the accuracy rises
  /// \param target query position
  /// \param dimension  0 = X, 1 = Y, 2 = Z
  /// \param points list of positions of sample points
  /// \param f list of function values on sample points
  /// \return gradient component at target on direction <dimension>
  double gradientAt(Point<double, D> target, int dimension,
                    const std::vector<Point<double, D>> &points,
                    const std::vector<double> &f);

  /// Compute laplacian weights to be used in differential equation problems,
  /// like in Poisson equation. The first point of input points is supposed to
  /// be the center of the stencil.
  /// \param points list of positions of sample points
  /// \return laplacian weights (weights are stored on the same order of points)
  std::vector<double>
  laplacianWeights(const std::vector<Point<double, D>> &points);

  /// Compute laplacian weights to be used in differential equation problems,
  /// like in Poisson equation. The weights are comuted for a regular stencil
  /// with center and four other points equidistants from the center point
  /// The target point must be in points[0]
  /// \param h euclidian distance between points regularly spaced:
  ///         2
  ///         |
  ///    3 -- 0 -- 1
  ///         |
  ///         4
  /// \return laplacian weights (weights are stored on the same order of points)
  std::vector<double> laplacianWeights(double h);

  double laplacianAt(Point<double, D> target,
                     const std::vector<Point<double, D>> &points,
                     const std::vector<double> &values);

private:
  Definitions::PolynomialType
      polynomialBaseType; //!< used to store current polynomial basis
  //!< augmentation
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qr;
  std::vector<Point<double, D>> _points;
  Point<double, D> _target;
  std::function<std::vector<double>(const double *, size_t)>
      base; //!< stores the current polynomial basis values
  std::function<std::vector<double>(int, const double *, size_t)>
      dbase; //!< stores the current differential valus for the polynoial basis
  std::function<std::vector<double>(int, const double *, size_t)>
      d2base; //!< stores the current differential valus for the polynoial basis
  RadialKernel<Point<double, D>> *_kernel =
      nullptr; //!< radial kernel used to compute differential weights

  static RadialKernel<Point<double, D>>* defaultKernel()
  {
    return new QuinticKernel<Point<double, D>>{};
  }
};
#include <blas/differential_rbf.inl>

/************************************************************************/
/************************************************************************/
// TODO: deprecated
// Template para Vector2, ou Vector3
template <class T> class RBF {
private:
  RadialKernel<T> *_kernel;

  // TODO Implementar a trie para armazenar pesos de stencils diferentes
  // Trie distanceDictionary;
  LinearVector _laplacian; // Hold laplacian weights for stenil
  LinearVector _gradientX,
      _gradientY; // Hold gradient weights for diffenret stencils
  LinearVector _divergent;

  LinearVector gradient(std::vector<T>, int);
  LinearVector gradientAt(T, std::vector<T>, int);

public:
  enum RBF_TYPE { GAUSS, IMQ, IQ };
  RBF();
  RBF(double);
  RBF(RadialKernel<T> *);
  RBF(RadialKernel<T> *, double);

  LinearVector divergent(std::vector<T>);
  LinearVector gradientX(std::vector<T>);
  LinearVector gradientY(std::vector<T>);
  LinearVector gradientZ(std::vector<T>);
  LinearVector divergentAt(T, std::vector<T>);
  LinearVector gradientXAt(T, std::vector<T>);
  LinearVector gradientYAt(T, std::vector<T>);
  LinearVector gradientZAt(T, std::vector<T>);
  LinearVector gradientXYAt(T, std::vector<T>);
  LinearVector gradientXZAt(T, std::vector<T>);
  LinearVector gradientYZAt(T, std::vector<T>);
  LinearVector RBFWeights(std::vector<T>);
  LinearVector laplacian(std::vector<T>);

  void setKernel(RadialKernel<T> *);
};

#include <blas/rbf.inl>
} // namespace furoo

#endif
