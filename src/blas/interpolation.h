#ifndef FUROO_BLAS_INTERPOLATION_H
#define FUROO_BLAS_INTERPOLATION_H

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <blas/linear_matrix.h>
#include <blas/linear_solver.h>
#include <blas/linear_vector.h>
#include <blas/rbf_kernels.h>
#include <cmath>
#include <common/code_profiler.h>
#include <common/definitions.h>
#include <geometry/point.h>
#include <numeric>
#include <vector>

namespace furoo {

// TODO: PLS CHANGE TO <typename T, int D> as ALL the other classes we have
template <class T> class Interpolant {
public:
  Interpolant() {}

  virtual ~Interpolant() {}

  virtual double interpolateAt(T, const std::vector<T> &,
                               const std::vector<double> &) = 0;

  virtual T interpolateAt(T target, const std::vector<T> &points,
                          const std::vector<T> &values) {
    NOT_IMPLEMENTED();
    return T();
  }

  virtual std::vector<double> interpolateAt(std::vector<T> target,
                                            const std::vector<T> &points,
                                            const std::vector<double> &values) {
    NOT_IMPLEMENTED();
    return std::vector<double>();
  }

  virtual void weights(const std::vector<T> &points,
                       const std::vector<double> &f,
                       std::vector<double> &w) = 0;

  virtual void computeWeights(const std::vector<T> &points,
                              const std::vector<double> &f) {
    NOT_IMPLEMENTED();
  }

  virtual double evaluate(T target) const {
    NOT_IMPLEMENTED();
    return 0.;
  }
};

template <class T> class KernelInterpolant : public Interpolant<T> {
private:
  RadialKernel<T> *_kernel;

public:
  KernelInterpolant(T d = T(1.))
      : KernelInterpolant(new TrilinearHatKernel<T>(d)) {}
  KernelInterpolant(RadialKernel<T> *k) { _kernel = k; }
  ~KernelInterpolant() {}
  double interpolateAt(T, const std::vector<T> &,
                       const std::vector<double> &) override;
  void weights(const std::vector<T> &, const std::vector<double> &f,
               std::vector<double> &) override;
};

template <class T> class RBFInterpolant : public Interpolant<T> {
private:
  RadialKernel<T> *_kernel;
  Definitions::PolynomialType _polynomialBaseType;
  std::function<std::vector<double>(const double *, size_t)> _base;
  std::function<std::vector<double>(int, const double *, size_t)> _dbase;
  std::vector<double> _weights;
  std::vector<T> _points;

public:
  RBFInterpolant(RadialKernel<T> &k); // TODO: k must be const

  RBFInterpolant();

  ~RBFInterpolant();

  RadialKernel<T> *kernel();

  std::function<std::vector<double>(const double *, size_t)> base() const;

  std::function<std::vector<double>(int, const double *, size_t)> dBase() const;

  /// Change the polynomial base augmentation used in rbf to compute
  /// differential weights.
  /// Currently supported polynomials: LINEAR, QUADRATIC
  /// \param type new polynomial type to be used
  void setBasePolynomial(Definitions::PolynomialType type);

  double interpolateAt(T, const std::vector<T> &,
                       const std::vector<double> &) override;

  T interpolateAt(T target, const std::vector<T> &points,
                  const std::vector<T> &values) override;

  std::vector<double> interpolateAt(std::vector<T> target,
                                    const std::vector<T> &points,
                                    const std::vector<double> &values) override;

  void weights(const std::vector<T> &points, const std::vector<double> &f,
               std::vector<double> &w) override;

  void computeWeights(const std::vector<T> &points,
                      const std::vector<double> &f) override;

  double evaluate(T target) const override;
};

class BilinearInterpolant : public Interpolant<Point2d> {
public:
  BilinearInterpolant() {}
  ~BilinearInterpolant() {}

  /// points are expected in the order:
  /// 2 ---- 3
  /// |      |
  /// 0 ---- 1
  /// \return
  double interpolateAt(Point2d, const std::vector<Point2d> &,
                       const std::vector<double> &) override;
  /** Computes weights for bilinear interpolation.
   *  Assumes points are ordered counter-clockwise starting from lower left.
   * \param target point
   * \param points (not used)
   * \param computed weights
   */
  void weights(const std::vector<Point2d> &, const std::vector<double> &,
               std::vector<double> &) override;
};

class BicubicInterpolant : public Interpolant<Point2d> {
public:
  BicubicInterpolant();
  ~BicubicInterpolant() {}
  double interpolateAt(Point2d, const std::vector<Point2d> &,
                       const std::vector<double> &) override;
  /** Computes the bicubic interpolation at the targets points. Assumes
   * everything is inside a unit box and the 4 corners are orderd in
   * counterclockwise starting from the lower left.
   * \param interpolation point [0,0]x[1,1]
   * \param function values
   * \param function gradient in x
   * \param function gradient in y
   * \param function gradient in xy
   * \param box size
   * \param interpolated gradient in y
   * \returns the interpolated function value
   */
  double interpolateAt(Point2d, const std::vector<double> &,
                       const std::vector<double> &, const std::vector<double> &,
                       const std::vector<double> &, Vector2d h,
                       double * = nullptr, double * = nullptr) const;
  void weights(const std::vector<Point2d> &, const std::vector<double> &,
               std::vector<double> &) override;

private:
  LinearMatrix m;
};

class MonotonicCatmullRom {

public:
  static void getFraction(double x, int low, int high, int *index, double *f);
  static double interpolate(double v0, double v1, double v2, double v3,
                            double f);
};

template <int D>
class ShepardInterpolant : public Interpolant<Point<double, D>> {
private:
  DistanceKernel<Point<double, D>> kernel;

public:
  ShepardInterpolant() = default;
  virtual ~ShepardInterpolant() = default;

  double interpolateAt(Point<double, D> target,
                       const std::vector<Point<double, D>> &points,
                       const std::vector<double> &f) override {
    if (points.size() == 0)
      return 0.;
    double sum = 0.;
    double wsum = 0.;
    for (size_t i = 0; i < points.size(); ++i) {
      if ((target - points[i]).length() < 1e-8) {
        return f[i];
      }
      double w = kernel.phi(points[i], target);
      wsum += w;
      sum += w * f[i];
    }
    ASSERT_FATAL(wsum > 0.);
    return sum / wsum;
  }
  void weights(const std::vector<Point<double, D>> &,
               const std::vector<double> &, std::vector<double> &) override{};
};

// TODO: MOVE CLASS IMPLEMENTATION TO INL FILE!
template <int D> class MLS : public Interpolant<Point<double, D>> {
public:
  explicit MLS(Definitions::PolynomialType baseType =
                   Definitions::PolynomialType::QUADRATIC) {
    base = Definitions::polynomial(baseType);
  }

  virtual ~MLS() = default;
  double interpolateAt(Point<double, D> target,
                       const std::vector<Point<double, D>> &points,
                       const std::vector<double> &f) override {
    // THROW(false, "YOU SHALL NOT MLS");
    THROW(points.size(), "MLS: empty list of points!");
    THROW(base != nullptr, "MLS: base not defined!");
    setBaseFromPointCount(points.size());
    auto targetBase = base(target.asConstArray(), D);
    // bases
    Eigen::MatrixXd X(points.size(), targetBase.size());
    for (size_t i = 0; i < points.size(); i++) {
      if (target == points[i])
        return f[i];
      auto pointBase = base(points[i].asConstArray(), D);
      for (size_t j = 0; j < pointBase.size(); j++)
        X(i, j) = pointBase[j];
    }
    Eigen::VectorXd w(points.size()), y(points.size());
    double maxF = -INFINITY, minF = INFINITY;
    for (size_t i = 0; i < points.size(); i++) {
      w(i) = kernel.phi(target, points[i]);
      y(i) = f[i];
      maxF = std::max(maxF, f[i]);
      minF = std::min(minF, f[i]);
    }
    Eigen::MatrixXd W = X.transpose() * w.asDiagonal();
    Eigen::MatrixXd A = W * X;
    Eigen::VectorXd b = W * y;

    Eigen::VectorXd alpha = A.fullPivLu().solve(b);
    double value = 0.;
    for (size_t i = 0; i < targetBase.size(); i++) {
      value += alpha[i] * targetBase[i];
    }

    if (std::isnan(value)) {
      std::cerr << "========================" << '\n';
      std::cerr << "w: " << w << std::endl;
      std::cerr << "w: " << W << std::endl;
      std::cerr << "X: " << X << std::endl;
      std::cerr << "y: " << y << std::endl;
      std::cerr << "A: " << A << std::endl << std::endl;
      std::cerr << "b: " << b << std::endl << std::endl;

      std::cerr << alpha << std::endl;
      for (size_t i = 0; i < targetBase.size(); i++)
        std::cerr << targetBase[i] << " ";
      std::cerr << std::endl;

      THROW(false, "MLS: NAN Value!");
    }
    //    THROW(std::isnan(value), "MLS: Nan value!");
    return clamp(value, minF, maxF);
  }

  void weights(const std::vector<Point<double, D>> &,
               const std::vector<double> &, std::vector<double> &) override{};

private:
  void setBasePolynomial(Definitions::PolynomialType type) {
    base = Definitions::polynomial(type);
  }
  void setBaseFromPointCount(size_t n) {
    if (D == 2) {
      if (n < 2)
        setBasePolynomial(Definitions::PolynomialType::CONSTANT);
      else if (n < 6)
        setBasePolynomial(Definitions::PolynomialType::LINEAR);
      else
        setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    } else {
      if (n < 3)
        setBasePolynomial(Definitions::PolynomialType::CONSTANT);
      else if (n < 9)
        setBasePolynomial(Definitions::PolynomialType::LINEAR);
      else
        setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
    }
    setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
  }
  GaussianKernel<Point<double, D>> kernel;
  std::function<std::vector<double>(const double *, size_t)> base;
};

#include <blas/interpolant.inl>

typedef ShepardInterpolant<2> ShepardInterpolant2;

} // namespace furoo

#endif // FUROO_BLAS_INTERPOLATION_H
