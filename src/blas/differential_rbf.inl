template <int D>
DifferentialRBF<D>::DifferentialRBF(RadialKernel<Point<double, D>> &kernel)
    : _kernel(&kernel) {
  setBasePolynomial(Definitions::PolynomialType::LINEAR);
  polynomialBaseType = Definitions::PolynomialType::LINEAR;
}

template <int D>
void DifferentialRBF<D>::setKernel(RadialKernel<Point<double, D>> &kernel) {
  _kernel = &kernel;
}

template <int D>
void DifferentialRBF<D>::setBasePolynomial(Definitions::PolynomialType type) {
  base = Definitions::polynomial(type);
  dbase = Definitions::polynomialDerivative(type);
  d2base = Definitions::polynomialSecondDerivative(type);
  polynomialBaseType = type;
}

template <int D>
double
DifferentialRBF<D>::divergentAt(Point<double, D> target,
                                const std::vector<Point<double, D>> &points,
                                const std::vector<Vector<double, D>> &f) {
  UNUSED_VARIABLE(target);
  UNUSED_VARIABLE(points);
  UNUSED_VARIABLE(f);
  return 0.;
  //    return gradientAt(target, 0, points, f) + gradientAt(target, 1,
  //    points, f);
}

template <int D>
double DifferentialRBF<D>::gradientAt(
    Point<double, D> target, int dimension,
    const std::vector<Point<double, D>> &originalPoints,
    const std::vector<double> &f) {

  size_t n = originalPoints.size();
  std::vector<Point<double, D>> points;
  for (auto &point : originalPoints)
    points.emplace_back(point - (target - Point<double, D>()));
  target = Point<double, D>(0.);
  if (D == 2) {
    if (n < 2)
      setBasePolynomial(Definitions::PolynomialType::CONSTANT);
    else if (n < 8)
      setBasePolynomial(Definitions::PolynomialType::LINEAR);
    else
      setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
  } else {
    if (n < 2)
      setBasePolynomial(Definitions::PolynomialType::CONSTANT);
    else if (n < 8)
      setBasePolynomial(Definitions::PolynomialType::LINEAR);
    else
      setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
  }
  // setBasePolynomial(Definitions::PolynomialType::LINEAR);
  auto targetBase = base(target.asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + targetBase.size(),
                                            points.size() + targetBase.size());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + targetBase.size());
  for (size_t i = 0; i < points.size(); i++) {
    b[i] = -(points[i] - target)[dimension] * _kernel->dphi(points[i], target);
    auto pointBase = base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points[i], points[j]);
  }
  // Setting rhs polynomial augmentation base
  auto diffBase = dbase(dimension, target.asConstArray(), D);
  for (size_t d = 0; d < diffBase.size(); d++)
    b[n + d] = diffBase[d];
  double grad = 0.;
  //  Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
  Eigen::VectorXd w = M.fullPivLu().solve(b);

#ifdef FUROO_PROFILE
  CodeProfiler::count("gradientAt");
  CodeProfiler::count(concat("gradientBase", targetBase.size()));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif

  for (size_t i = 0; i < n; i++)
    grad += w[i] * f[i];
  //    auto dtarget = base(target.asConstArray(), D);
  //    for (size_t i = 0; i < targetBase.size(); i++)
  //      grad += target[i] * w[n + i];
  return grad;
}

template <int D>
std::vector<double> DifferentialRBF<D>::laplacianWeights(
    const std::vector<Point<double, D>> &originalPoints) {
  std::vector<Point<double, D>> points;
  auto center = originalPoints[0];
  for (size_t i = 0; i < originalPoints.size(); i++)
    points.emplace_back(originalPoints[i] - (center - Point<double, D>()));
  center = Point<double, D>(0.);
  size_t n = points.size();
  /*if (D == 2) {
    if (n < 2)
      setBasePolynomial(Definitions::PolynomialType::CONSTANT);
    else if (n < 8)
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
  }*/
  setBasePolynomial(Definitions::PolynomialType::QUADRATIC);
  auto laplBaseX = d2base(0, points[0].asConstArray(), D);
  auto laplBaseY = d2base(1, points[0].asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + laplBaseX.size(),
                                            points.size() + laplBaseX.size());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + laplBaseX.size());

  for (size_t i = 0; i < points.size(); i++) {
    //    b[i] = _kernel->d2phi(points[0], points[i]) -
    //           _kernel->dphi(points[0], points[i]);
    b[i] = _kernel->d2phi(Point<double, D>(0.), points[i]) +
           _kernel->dphi(Point<double, D>(0.), points[i]);
    if (D == 3)
      b[i] += _kernel->dphi(Point<double, D>(0.), points[i]);
    auto pointBase = base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));

    std::vector<std::vector<double>> laplBases;
    for (int d = 0; d < D; d++)
      laplBases.push_back(d2base(d, points[i].asConstArray(), D));
    for (size_t j = 0; j < laplBaseX.size(); j++) {
      b[n + j] = 0.;
      for (int d = 0; d < D; d++)
        b[n + j] += laplBases[d][j];
    }
  }

  // Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
  Eigen::VectorXd w = M.fullPivLu().solve(b);

#ifdef FUROO_PROFILE
  CodeProfiler::count("laplacianWeights");
  CodeProfiler::count(concat("laplacianBase", laplBaseX.size()));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif

  if (w.hasNaN()) {
    THROW(false,
          "DifferentialRBF::laplacianWeights final weights vector has NaN.");
  }
  std::vector<double> weights;
  for (int i = 0; i < w.rows(); i++)
    weights.emplace_back(w[i]);
  return weights;
}

template <int D>
double DifferentialRBF<D>::laplacianAt(
    Point<double, D> target,
    const std::vector<Point<double, D>> &originalPoints,
    const std::vector<double> &values) {
  std::vector<Point<double, D>> points;
  auto center = target;
  for (size_t i = 0; i < originalPoints.size(); i++)
    points.emplace_back(originalPoints[i] - (center - Point<double, D>()));
  center = Point<double, D>(0.);
  size_t n = points.size();
  if (D == 2) {
    if (n < 2)
      setBasePolynomial(Definitions::PolynomialType::CONSTANT);
    else if (n < 8)
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
  auto laplBaseX = d2base(0, center.asConstArray(), D);
  auto laplBaseY = d2base(1, center.asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(points.size() + laplBaseX.size(),
                                            points.size() + laplBaseX.size());
  Eigen::VectorXd b = Eigen::VectorXd::Zero(points.size() + laplBaseX.size());

  for (size_t i = 0; i < points.size(); i++) {
    //    b[i] = _kernel->d2phi(points[0], points[i]) -
    //           _kernel->dphi(points[0], points[i]);
    b[i] = _kernel->d2phi(Point<double, D>(0.), points[i]) +
           _kernel->dphi(Point<double, D>(0.), points[i]);
    if (D == 3)
      b[i] += _kernel->dphi(Point<double, D>(0.), points[i]);
    auto pointBase = base(points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(points.at(i), points.at(j));

    auto laplBaseX = d2base(0, points[i].asConstArray(), D);
    auto laplBaseY = d2base(1, points[i].asConstArray(), D);
    for (size_t d = 0; d < laplBaseX.size(); d++) {
      b[n + d] = laplBaseX[d] + laplBaseY[d];
    }
  }
#ifdef FUROO_PROFILE
  CodeProfiler::count("laplacianAt");
  CodeProfiler::count(concat("laplacianBase", laplBaseX.size()));
  CodeProfiler::count(concat("EigenSolve", b.size()));
#endif
  //  Eigen::VectorXd w = M.colPivHouseholderQr().solve(b);
  Eigen::VectorXd w = M.fullPivLu().solve(b);
  if (w.hasNaN()) {
    THROW(false,
          "DifferentialRBF::laplacianWeights final weights vector has NaN.");
  }
  double laplacianValue = 0.;
  for (int i = 0; i < values.size(); i++)
    laplacianValue += w[i] * values[i];
  return laplacianValue;
}

template <int D>
std::vector<double> DifferentialRBF<D>::laplacianWeights(double h) {
  if (D == 2) {
    double lambda =
        (_kernel->lphi(h) - _kernel->lphi(0.)) /
        (5. * _kernel->phi(0.) - 8. * _kernel->phi(h) + _kernel->phi(2. * h) +
         2. * _kernel->phi(std::sqrt(2.) * h));
    std::vector<double> weights(5, lambda);
    weights[0] = -4. * lambda;
    return weights;
  } else {
    double p = 1 / (h * h);
    return {-6.0 / (h * h), p, p, p, p, p, p};
  }
}

template <int D>
void DifferentialRBF<D>::setTarget(const Point<double, D> target) {
  _target = target;
}

template <int D>
void DifferentialRBF<D>::addStencilPoint(Point<double, D> point) {
  _points.emplace_back(point);
}

template <int D> void DifferentialRBF<D>::setStencil() {
  THROW(_points.size() > 0,
        "DifferentailRBF::setStencil points of the stencil not set");
  auto targetBase = base(_target.asConstArray(), D);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(_points.size() + targetBase.size(),
                                            _points.size() + targetBase.size());
  size_t n = _points.size();
  for (size_t i = 0; i < _points.size(); i++) {
    auto pointBase = base(_points[i].asConstArray(), D);
    for (size_t d = 0; d < pointBase.size(); d++)
      M(n + d, i) = M(i, n + d) = pointBase[d];
    for (size_t j = i; j < _points.size(); j++)
      M(i, j) = M(j, i) = _kernel->phi(_points.at(i), _points.at(j));
  }
  _qr = M.colPivHouseholderQr();
}

template <int D>
double DifferentialRBF<D>::gradientAt(int dimension,
                                      const std::vector<double> &f) {
  auto targetBase = base(_target.asConstArray(), D);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(_points.size() + targetBase.size());
  size_t n = _points.size();
  for (size_t i = 0; i < _points.size(); i++)
    b[i] = (_points.at(i) - _target)[dimension] *
           _kernel->dphi(_points.at(i), _target);
  // Setting rhs polynomial augmentation base
  auto diffBase = dbase(dimension, _target.asConstArray(), D);
  for (size_t d = 0; d < diffBase.size(); d++)
    b[n + d] = diffBase[d];
}

template <int D> RadialKernel<Point<double, D>> *DifferentialRBF<D>::kernel() {
  return _kernel;
}
