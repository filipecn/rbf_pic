#ifndef FUROO_BLAS_RBF_KERNELS_H
#define FUROO_BLAS_RBF_KERNELS_H

#include <algorithm>
#include <cmath>
#include <common/debug.h>
#include <geometry/numeric.h>
#include <iostream>

namespace furoo {
// T: Dimensional point
template<class T>
class RadialKernel {
 public:
  explicit RadialKernel(double h = 1.0) : _h(h) {};
  virtual double phi(T, T) = 0;
  virtual double dphi(T, T) = 0;
  virtual double d2phi(T, T) = 0;
  virtual double phi(double d) = 0;
  virtual double dphi(double d) = 0;
  virtual double d2phi(double d) = 0;
  virtual double lphi(double d) { return d2phi(d) + dphi(d); };
  virtual void setH(double h) { _h = h; };
  virtual double getH() { return _h; };
 protected:
  double _h = 1.0;
};

template<class T>
class GaussianKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  explicit GaussianKernel(double = 1.);
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
  //void setH(double) override;
//  double getH() override;
};

template<class T>
class GaussianInverseKernel : public RadialKernel<T> {
 public:
  explicit GaussianInverseKernel(double = 1.);
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
};

template<class T>
class IMQKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
//  void setH(double) override;
//  double getH() override;
};

template<class T>
class CubicKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
//  void setH(double) override {}
//  double getH() override { return 0.0; }
};

template<class T>
class QuinticKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
//  void setH(double) override {}
//  double getH() override { return 0.0; }
};

template<class T>
class HepticKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double d) override;
  double dphi(double d) override;
  double d2phi(double d) override;
//  void setH(double) override {}
//  double getH() override { return 0.0; }
};

template<class T>
class WendlandKernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double) override;
  double dphi(double) override;
  double d2phi(double) override;
//  void setH(double) override;
//  double getH() override;
};

template<class T>
class Wendland32Kernel : public RadialKernel<T> {
// private:
//  double h = 1.0;

 public:
  double phi(T, T) override;
  double dphi(T, T) override;
  double d2phi(T, T) override;
  double phi(double) override;
  double dphi(double) override;
  double d2phi(double) override;
//  void setH(double) override;
//  double getH() override;
};

template<typename T>
class TrilinearHatKernel : public RadialKernel<T> {
 public:
  /**
   * \param d spacing / cubic root of volume
   */
  TrilinearHatKernel(T d = T(1.)) : _d(d) {}
  virtual ~TrilinearHatKernel() {}
  double phi(T xp, T xi) override {
    return h((xp[0] - xi[0]) / _d[0]) *
        ((xp.dimension() > 1) ? h((xp[1] - xi[1]) / _d[1]) : 1.) *
        ((xp.dimension() > 2) ? h(xp[2] - xi[2]) / _d[2] : 1.);
  }
  double dphi(T xp, T xi) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(xp);
    UNUSED_VARIABLE(xi);
    return 0.;
  }
  double d2phi(T xp, T xi) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(xp);
    UNUSED_VARIABLE(xi);
    return 0.;
  }
  double phi(double d) override {
    UNUSED_VARIABLE(d);
    return 0.;
  }
  double dphi(double d) override {
    UNUSED_VARIABLE(d);
    return 0.;
  }
  double d2phi(double d) override {
    UNUSED_VARIABLE(d);
    return 0.;
  }
//  void setH(double d) override { UNUSED_VARIABLE(d); }
//  double getH() override { return 0.; }

 private:
  T _d;
  double h(double r) const {
    if (r >= 0. && r <= 1.)
      return 1. - r;
    if (r >= -1. && r < 0.)
      return 1. + r;
    return 0.;
  }
};

template<typename T>
class QuadraticBSplineKernel : public RadialKernel<T> {
 public:
  QuadraticBSplineKernel(T d = T(1.)) : _d(d) {}
  virtual ~QuadraticBSplineKernel() {}
  double phi(T xp, T xi) override {
    return h((xp[0] - xi[0]) / _d[0]) *
        ((xp.dimension() > 1) ? h((xp[1] - xi[1]) / _d[1]) : 1.) *
        ((xp.dimension() > 2) ? h(xp[2] - xi[2]) / _d[2] : 1.);
  }
  double dphi(T xp, T xi) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(xp);
    UNUSED_VARIABLE(xi);
    return 0.;
  }
  double d2phi(T xp, T xi) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(xp);
    UNUSED_VARIABLE(xi);
    return 0.;
  }
  double phi(double d) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(d);
    return 0.;
  }
  double dphi(double d) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(d);
    return 0.;
  }
  double d2phi(double d) override {
    NOT_IMPLEMENTED();
    UNUSED_VARIABLE(d);
    return 0.;
  }
//  void setH(double d) override { UNUSED_VARIABLE(d); }
//  double getH() override { return 0.; }

 private:
  T _d;
  double h(double r) const {
    if (r >= -1.5 && r < -.5)
      return .5 * SQR(r + 1.5);
    if (r >= -.5 && r < .5)
      return .75 - SQR(r);
    if (r >= .5 && r < 1.5)
      return .5 * SQR(1.5 - r);
    return 0.;
  }
};

template<class T>
class DistanceKernel : public RadialKernel<T> {
// private:
//  double h = 3;

 public:
  DistanceKernel(double h = 3.) : RadialKernel<T>(h) {}
  virtual double phi(T, T) override;
  virtual double dphi(T, T) override { return 0.; }
  virtual double d2phi(T, T) override { return 0.; }
  double phi(double) override { return 0.; }
  double dphi(double) override { return 0.; };
  double d2phi(double) override { return 0.; };
//  virtual void setH(double) override {}
//  virtual double getH() override { return 0.; }
};

#include <blas/rbf_kernels.inl>
} // namespace furoo
#endif
