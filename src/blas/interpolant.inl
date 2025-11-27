template <class T> RBFInterpolant<T>::RBFInterpolant(RadialKernel<T> &k) {
  _kernel = &k;
}

template <class T> RBFInterpolant<T>::RBFInterpolant() {
  static QuinticKernel<T> defaultKernel;
  _kernel = &defaultKernel;
}

template <class T> RBFInterpolant<T>::~RBFInterpolant() {}

template <class T>
void RBFInterpolant<T>::setBasePolynomial(Definitions::PolynomialType type) {
  _base = Definitions::polynomial(type);
  _dbase = Definitions::polynomialDerivative(type);
  _polynomialBaseType = type;
}

template <class T>
double RBFInterpolant<T>::interpolateAt(T target, const std::vector<T> &points,
                                        const std::vector<double> &f) {
  int D = target.dimension();
  size_t n = points.size();
  THROW(n > 0, "RBFInterpolant::interpolateAt not enough points");
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
  auto targetBase = _base(target.asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + targetBase.size(),
                                            points.size() + targetBase.size());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + targetBase.size());
  Eigen::VectorXd _f = Eigen::VectorXd::Zero(
      f.size() + targetBase.size()); // Rbf values for target - stencil points

  double maxF = -INFINITY, minF = INFINITY;
  for (size_t i = 0; i < points.size(); i++) {
    b[i] = f[i];
    maxF = std::max(maxF, f[i]);
    minF = std::min(minF, f[i]);
    auto pointBase = _base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));
    _f[i] = _kernel->phi(target, points.at(i));
  }
  // Setting rhs polynomial augmentation base
  for (size_t d = 0; d < targetBase.size(); d++) {
    b[n + d] = 0;
    _f[n + d] = targetBase[d];
  }

  // Eigen::VectorXd w = M.fullPivLu().solve(b);
  Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
  // Eigen::PartialPivLU<Eigen::MatrixXd> solver(M);
  // Eigen::VectorXd w = solver.solve(b);
#ifdef FUROO_PROFILE
  CodeProfiler::count("RBFinterpolateAt");
  CodeProfiler::count(concat("interpolateBase", targetBase.size()));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif

  // Evaluating the point over the computed weights
  double value = w.transpose() * _f;
  return clamp(value, minF, maxF);
}

template <class T>
T RBFInterpolant<T>::interpolateAt(T target, const std::vector<T> &points,
                                   const std::vector<T> &values) {
  size_t n = points.size();
  int D = target.dimension();
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
  auto targetBase = _base(target.asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + targetBase.size(),
                                            points.size() + targetBase.size());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + targetBase.size());
  Eigen::VectorXd _f = Eigen::VectorXd::Zero(
      values.size() +
      targetBase.size()); // Rbf values for target - stencil points

  for (size_t i = 0; i < points.size(); i++) {
    auto pointBase = _base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));
    _f[i] = _kernel->phi(target, points.at(i));
  }
  // Setting rhs polynomial augmentation base
  for (size_t d = 0; d < targetBase.size(); d++) {
    b[n + d] = 0;
    _f[n + d] = targetBase[d];
  }
  T ans;
  for (int d = 0; d < D; d++) {
    double maxF = -INFINITY, minF = INFINITY;
    for (size_t i = 0; i < points.size(); i++) {
      b[i] = values[i][d];
      maxF = std::max(maxF, values[i][d]);
      minF = std::min(minF, values[i][d]);
    }
    Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
    if ((M * w - b).norm() / b.norm() > 1e-3) {
      std::cerr << "colPiv error too high, trying FullPiv LU\n";
      auto MS = M.fullPivLu();
      Eigen::VectorXd w = MS.solve(b);
    }
#ifdef FUROO_PROFILE
    CodeProfiler::count("RBFinterpolateAt");
    CodeProfiler::count(concat("interpolateBase", targetBase.size()));
    CodeProfiler::count(concat("EigenSolve", b.size()));
#endif
    double value = w.transpose() * _f;
    ans[d] = clamp(value, minF, maxF);
  }
  return ans;
}

template <class T>
std::vector<double>
RBFInterpolant<T>::interpolateAt(std::vector<T> target,
                                 const std::vector<T> &points,
                                 const std::vector<double> &values) {
  if (target.size() == 0)
    return std::vector<double>();

  size_t n = points.size();
  int D = target[0].dimension();
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
  // setBasePolynomial(Definitions::PolynomialType::QUADRATIC);

  auto targetBaseSize = _base(target[0].asConstArray(), D).size();
  std::vector<std::vector<double>> targetBases;
  for (T &ts : target) {
    targetBases.push_back(_base(ts.asConstArray(), D));
  }

  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + targetBaseSize,
                                            points.size() + targetBaseSize);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + targetBaseSize);
  std::vector<Eigen::VectorXd> pointF(
      target.size(), Eigen::VectorXd::Zero(values.size() + targetBaseSize));

  double maxF = -INFINITY, minF = INFINITY;
  for (size_t i = 0; i < points.size(); i++) {
    auto pointBase = _base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));

    for (int tid = 0; tid < target.size(); tid++)
      pointF[tid][i] = _kernel->phi(target[tid], points.at(i));

    b[i] = values[i];
    maxF = std::max(maxF, values[i]);
    minF = std::min(minF, values[i]);
  }
  // Setting rhs polynomial augmentation base
  for (size_t d = 0; d < targetBaseSize; d++) {
    b[n + d] = 0;
    for (int ti = 0; ti < target.size(); ti++) {
      pointF[ti][n + d] = targetBases[ti][d];
    }
  }
  std::vector<double> ans;
  Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
  double sysError = (M * w - b).norm() / b.norm();
  if (sysError > 5e-3) {
    std::cerr << "colPiv error too high (" << sysError
              << "), trying FullPiv LU (";
    auto MS = M.fullPivLu();
    Eigen::VectorXd w = MS.solve(b);
    sysError = (M * w - b).norm() / b.norm();
    std::cerr << sysError << ")\n";
  }
#ifdef FUROO_PROFILE
  CodeProfiler::count("RBFinterpolateAt");
  CodeProfiler::count(concat("interpolateBase", targetBase.size()));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif

  for (int tid = 0; tid < target.size(); tid++) {
    double value = w.transpose() * pointF[tid];
    ans.emplace_back(clamp(value, minF, maxF));
  }

  return ans;
}

template <class T>
void RBFInterpolant<T>::weights(const std::vector<T> &points,
                                const std::vector<double> &f,
                                std::vector<double> &w) {
  THROW(points.size() && points.size() == f.size(),
        "RBFInterpolant<>::weights empty stencil");
  int D = points[0].dimension();
  size_t n = points.size();
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
  _points = points;
  size_t baseSize = _base(_points[0].asConstArray(), D).size();
  Eigen::MatrixXd M =
      Eigen::MatrixXd::Zero(points.size() + baseSize, points.size() + baseSize);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + baseSize);
  Eigen::VectorXd _f = Eigen::VectorXd::Zero(f.size() + baseSize);
  for (size_t i = 0; i < points.size(); i++) {
    b[i] = f[i];
    auto pointBase = _base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));
  }
  // Setting rhs polynomial augmentation base
  for (size_t d = 0; d < baseSize; d++) {
    b[n + d] = 0;
  }
  Eigen::VectorXd _w = M.fullPivLu().solve(b);
#ifdef FUROO_PROFILE
  CodeProfiler::count("RBFweights");
  CodeProfiler::count(concat("interpolateBase", baseSize));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif
  w.clear();
  for (size_t i = 0; i < _w.size(); i++)
    w.emplace_back(_w[i]);
}

template <class T>
void RBFInterpolant<T>::computeWeights(const std::vector<T> &points,
                                       const std::vector<double> &f) {
  weights(points, f, _weights);
}

template <class T> double RBFInterpolant<T>::evaluate(T target) const {
  int D = target.dimension();
  size_t n = _points.size();
  auto targetBase = _base(target.asConstArray(), D);
  THROW(_points.size() && _points.size() == _weights.size() - targetBase.size(),
        "RBFInterpolant<>::evaluate invalid stencil");
  Eigen::VectorXd _f =
      Eigen::VectorXd::Zero(_points.size() + targetBase.size());
  for (size_t i = 0; i < _points.size(); i++)
    _f[i] = _kernel->phi(target, _points.at(i));
  for (size_t d = 0; d < targetBase.size(); d++)
    _f[n + d] = targetBase[d];
  double sum = 0;
  for (size_t i = 0; i < _weights.size(); i++)
    sum += _weights[i] * _f[i];
#ifdef FUROO_PROFILE
  CodeProfiler::count("RBFevaluateAt");
#endif
  return sum;
}

template <class T> RadialKernel<T> *RBFInterpolant<T>::kernel() {
  return _kernel;
}

template <class T>
std::function<std::vector<double>(const double *, size_t)>
RBFInterpolant<T>::base() const {
  return _base;
}

template <class T>
std::function<std::vector<double>(int, const double *, size_t)>
RBFInterpolant<T>::dBase() const {
  return _dbase;
}

template <class T>
double KernelInterpolant<T>::interpolateAt(T target,
                                           const std::vector<T> &points,
                                           const std::vector<double> &f) {
  double sum = 0.;
  double W = 0;
  for (size_t i = 0; i < points.size(); i++) {
    double kVal = _kernel->phi(points[i], target);
    sum += f[i] * kVal;
    if (kVal != 0.)
      W++;
  }
  std::cout << "  W: " << W << '\n';
  return sum / W;
}

template <class T>
void KernelInterpolant<T>::weights(const std::vector<T> &points,
                                   const std::vector<double> &f,
                                   std::vector<double> &w) {
  UNUSED_VARIABLE(points);
  UNUSED_VARIABLE(w);
}
