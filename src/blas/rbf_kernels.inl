/*****************************************************************************/
// template<class T>
// void GaussianKernel<T>::setH(double newH) { h = newH; }
// template<class T>
// double GaussianKernel<T>::getH() { return h; }
template <class T> double GaussianKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return std::exp(-d2 * h2);
}
template <class T>
GaussianKernel<T>::GaussianKernel(double h_) : RadialKernel<T>(h_) {}
template <class T> double GaussianKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return (-2 * h2) * std::exp(-d2 * h2);
}
template <class T> double GaussianKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return (2 * h2) * std::exp(-d2 * h2) * (2 * d2 * h2 - 1);
}
template <class T> double GaussianKernel<T>::phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return std::exp(-d2 * h2);
}

template <class T> double GaussianKernel<T>::dphi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return (-2 * h2) * std::exp(-d2 * h2);
}

template <class T> double GaussianKernel<T>::d2phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return (2 * h2) * std::exp(-d2 * h2) * (2 * d2 * h2 - 1);
}

/*****************************************************************************/
// template<class T>
// void GaussianInverseKernel<T>::setH(double newH) { h = newH; }
// template<class T>
// double GaussianInverseKernel<T>::getH() { return h; }
template <class T> double GaussianInverseKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return std::exp(-d2 / h2);
}
template <class T>
GaussianInverseKernel<T>::GaussianInverseKernel(double h_)
    : RadialKernel<T>(h_) {}
template <class T> double GaussianInverseKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return (-2 / h2) * std::exp(-d2 / h2);
}
template <class T> double GaussianInverseKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return (2 / h2) * std::exp(-d2 / h2) * (2 * d2 / h2 - 1);
}
template <class T> double GaussianInverseKernel<T>::phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return std::exp(-d2 / h2);
}

template <class T> double GaussianInverseKernel<T>::dphi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return (-2 / h2) * std::exp(-d2 / h2);
}

template <class T> double GaussianInverseKernel<T>::d2phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return (2 / h2) * std::exp(-d2 / h2) * (2 * d2 / h2 - 1);
}

/*****************************************************************************/
// template<class T>
// void IMQKernel<T>::setH(double newH) { h = newH; }
// template<class T>
// double IMQKernel<T>::getH() { return h; }
template <class T> double IMQKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return 1 / std::sqrt(1 + h2 * d2);
}
template <class T> double IMQKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;
  double den = (1 + d2 * h2);

  return -(h2) / (den * std::sqrt(den));
}
template <class T> double IMQKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;
  double den = (1 + d2 * h2);

  return h2 * (2 * h2 * d2 - 1) / (den * den * std::sqrt(den));
}
template <class T> double IMQKernel<T>::phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  return 1 / std::sqrt(1 + h2 * d2);
}
template <class T> double IMQKernel<T>::dphi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  double den = (1 + d2 * h2);

  return -(h2) / (den * std::sqrt(den));
}
template <class T> double IMQKernel<T>::d2phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;
  double den = (1 + d2 * h2);
  return h2 * (2 * h2 * d2 - 1) / (den * den * std::sqrt(den));
}

/*****************************************************************************/
template <class T> double CubicKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;
  return d2 * d;
}
template <class T> double CubicKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  return 3 * d;
}
template <class T> double CubicKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  return 6 * d;
}
template <class T> double CubicKernel<T>::phi(double d) {
  double d2 = d * d;
  return d2 * d;
}
template <class T> double CubicKernel<T>::dphi(double d) { return 3 * d; }
template <class T> double CubicKernel<T>::d2phi(double d) { return 6 * d; }

/*****************************************************************************/
template <class T> double QuinticKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return d2 * d2 * d;
}
template <class T> double QuinticKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return 5 * d2 * d;
}
template <class T> double QuinticKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return 20 * d2 * d;
}
template <class T> double QuinticKernel<T>::phi(double d) {
  double d2 = d * d;
  return d2 * d2 * d;
}
template <class T> double QuinticKernel<T>::dphi(double d) {
  double d2 = d * d;
  return 5 * d2 * d;
}
template <class T> double QuinticKernel<T>::d2phi(double d) {
  double d2 = d * d;
  return 20 * d2 * d;
}

/*****************************************************************************/
template <class T> double HepticKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return d2 * d2 * d2 * d;
}
template <class T> double HepticKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return 7 * d2 * d2 * d;
}
template <class T> double HepticKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d;

  return 42 * d2 * d2 * d;
}
template <class T> double HepticKernel<T>::phi(double d) {
  double d2 = d * d;
  return d2 * d2 * d2 * d;
}
template <class T> double HepticKernel<T>::dphi(double d) {
  double d2 = d * d;
  return 7 * d2 * d2 * d;
}
template <class T> double HepticKernel<T>::d2phi(double d) {
  double d2 = d * d;
  return 42 * d2 * d2 * d;
}

/*****************************************************************************/
// template<class T>
// void WendlandKernel<T>::setH(double newH) { h = newH; }
// template<class T>
// double WendlandKernel<T>::getH() { return h; }
template <class T> double WendlandKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  return std::pow(std::max(0.0, 1.0 - this->_h * d), 8) *
         (32.0 * CUBE(this->_h * d) + 25.0 * SQR(this->_h * d) +
          8.0 * this->_h * d + 1.0);
}
template <class T> double WendlandKernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  return -22.0 * SQR(this->_h) *
             std::pow(std::max(0.0, 1.0 - this->_h * d), 7) +
         (16.0 * SQR(this->_h * d) + 7.0 * this->_h * d + 1.0);
}
template <class T> double WendlandKernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  return 528.0 * std::pow(this->_h, 4) *
             std::pow(std::max(0.0, 1.0 - this->_h * d), 6) +
         (6.0 * this->_h * d + 1.0);
}
template <class T> double WendlandKernel<T>::phi(double d) {
  return std::pow(std::max(0.0, 1.0 - this->_h * d), 8) *
         (32.0 * CUBE(this->_h * d) + 25.0 * SQR(this->_h * d) +
          8.0 * this->_h * d + 1.0);
}
template <class T> double WendlandKernel<T>::dphi(double d) {
  return -22.0 * SQR(this->_h) *
             std::pow(std::max(0.0, 1.0 - this->_h * d), 7) +
         (16.0 * SQR(this->_h * d) + 7.0 * this->_h * d + 1.0);
}
template <class T> double WendlandKernel<T>::d2phi(double d) {
  return 528.0 * std::pow(this->_h, 4) *
             std::pow(std::max(0.0, 1.0 - this->_h * d), 6) +
         (6.0 * this->_h * d + 1.0);
}

/*****************************************************************************/
// template<class T>
// void Wendland32Kernel<T>::setH(double newH) { h = newH; }
// template<class T>
// double Wendland32Kernel<T>::getH() { return h; }
template <class T> double Wendland32Kernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return std::pow(1 - d * this->_h, 6) + (35 * d2 * h2 + 18 * d * this->_h + 3);
}
template <class T> double Wendland32Kernel<T>::dphi(T v, T w) {
  double d = (v - w).length();
  double h2 = this->_h * this->_h;

  return -56 * h2 * (5 * this->_h * d + 1) * std::pow(1 - d * this->_h, 5);
}
template <class T> double Wendland32Kernel<T>::d2phi(T v, T w) {
  double d = (v - w).length();
  double d2 = d * d, h2 = this->_h * this->_h;

  return 56 * h2 * (35 * d2 * h2 - 4 * this->_h * d - 1) *
         std::pow(1 - d * this->_h, 4);
}
template <class T> double Wendland32Kernel<T>::phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;

  return std::pow(1 - d * this->_h, 6) + (35 * d2 * h2 + 18 * d * this->_h + 3);
}
template <class T> double Wendland32Kernel<T>::dphi(double d) {
  double h2 = this->_h * this->_h;

  return -56 * h2 * (5 * this->_h * d + 1) * std::pow(1 - d * this->_h, 5);
}
template <class T> double Wendland32Kernel<T>::d2phi(double d) {
  double d2 = d * d, h2 = this->_h * this->_h;

  return 56 * h2 * (35 * d2 * h2 - 4 * this->_h * d - 1) *
         std::pow(1 - d * this->_h, 4);
}

/*****************************************************************************/

/*****************************************************************************/
template <class T> double DistanceKernel<T>::phi(T v, T w) {
  double d = (v - w).length();
  return std::pow(1 / d, this->_h);
}
