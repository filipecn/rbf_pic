template<class T>
RBF<T>::RBF() {
  _kernel = new GaussianKernel<T>();
  _gradientX.resize(3);
  _gradientY.resize(3);
}

template<class T>
RBF<T>::RBF(double h) : RBF() { _kernel->setH(h); }

template<class T>
RBF<T>::RBF(RadialKernel <T> *kernel) : RBF() {
  _kernel = kernel;
}

template<class T>
RBF<T>::RBF(RadialKernel <T> *newKernel, double h) {
  _kernel = newKernel;
  _kernel->setH(h);
}
template<class T>
void RBF<T>::setKernel(RadialKernel <T> *newKernel) {
  _kernel = newKernel;
}

template<class T>
LinearVector RBF<T>::laplacian(std::vector<T> neighbors) {
  if (_laplacian.size() == 0) {
    // std::cerr << "Laplacian" << std::endl;
    LinearMatrix M(neighbors.size() + 3);
    LinearVector b(neighbors.size() + 3);
    _laplacian = LinearVector(neighbors.size());
    size_t n = neighbors.size();

    for (size_t i = 0; i < n; i++) {
      b[i] = _kernel->d2phi(neighbors.at(0), neighbors.at(i)) +
          _kernel->dphi(neighbors.at(0), neighbors.at(i));

      T p = neighbors[i];
      M(n, i) = 1;
      M(i, n) = 1;
      for (size_t d = 1; d <= p.dimension(); d++) {
        M(n + d, i) = p[d - 1];
        M(i, n + d) = p[d - 1];
      }
      for (size_t j = i; j < neighbors.size(); j++) {
        M(i, j) = _kernel->phi(neighbors.at(i), neighbors.at(j));
        M(j, i) = _kernel->phi(neighbors.at(i), neighbors.at(j));
      }
    }
    // Solve RBF system and convert to vector
    GaussJordanLinearSolver solver;
    // std::cerr << M << std::endl;
    // std::cerr << b << std::endl;
    solver.solve(&M, &b, &_laplacian);
    _laplacian = b;
    for (unsigned i = 0; i < _laplacian.size(); i++)
      _laplacian[i] = -_laplacian[i];
  }
  return _laplacian;
}

template<class T>
LinearVector RBF<T>::gradientAt(T point, std::vector<T> neighbors, int dim) {
  // std::cerr << "Computeing gradient at " << std::endl;
  LinearMatrix M(neighbors.size());
  LinearVector b(neighbors.size());
  for (size_t i = 0; i < neighbors.size(); i++) {
    b[i] =
        (point - neighbors.at(i))[dim] * _kernel->dphi(point, neighbors.at(i));
    for (size_t j = i; j < neighbors.size(); j++) {
      M(i, j) = _kernel->phi(neighbors.at(i), neighbors.at(j));
      M(j, i) = _kernel->phi(neighbors.at(i), neighbors.at(j));
    }
  }

  GaussJordanLinearSolver solver;
  LinearVector _g(neighbors.size());
  solver.solve(&M, &b, &_g);
  _g = b;
  return _g;
}

template<class T>
LinearVector RBF<T>::gradient(std::vector<T> neighbors, int dim) {
  // std::cerr << "Computig gradient coord-" << dim << std::endl;
  LinearMatrix M(neighbors.size() + 3);
  LinearVector b(neighbors.size() + 3);
  LinearVector _gradient = LinearVector(neighbors.size() + 1);
  size_t n = neighbors.size();

  for (size_t i = 0; i < neighbors.size(); i++) {
    b[i] = (neighbors.at(i) - neighbors.at(0))[dim] *
        _kernel->dphi(neighbors.at(i), neighbors.at(0));

    T p = neighbors[i];
    for (size_t d = 1; d <= p.dimension(); d++) {
      M(n + d, i) = p[d - 1];
      M(i, n + d) = p[d - 1];
    }
    M(n, i) = 1;
    M(i, n) = 1;
    for (size_t j = i; j < neighbors.size(); j++) {
      M(i, j) = _kernel->phi(neighbors.at(i), neighbors.at(j));
      M(j, i) = _kernel->phi(neighbors.at(i), neighbors.at(j));
    }
  }
  b[n + dim + 1] = 1;

  GaussJordanLinearSolver solver;
  // std::cerr << M << std::endl;
  // std::cerr << b << std::endl;
  solver.solve(&M, &b, &_gradient);
  _gradient = b;
  // std::cerr << _gradient << std::endl;
  // std::cerr << b << std::endl;
  return _gradient;
}

template<class T>
LinearVector RBF<T>::gradientX(std::vector<T> neighbors) {
  // int n = neighbors.size() - 3;

  // if (_gradientX.size() == 0) {
  _gradientX = gradient(neighbors, 0);
  // }
  return _gradientX;
}

template<class T>
LinearVector RBF<T>::gradientXAt(T point, std::vector<T> neighbors) {
  // int n = neighbors.size() - 3;

  // if (_gradientX.size() == 0) {
  _gradientX = gradientAt(point, neighbors, 0);
  // }
  return _gradientX;
}

template<class T>
LinearVector RBF<T>::gradientY(std::vector<T> neighbors) {
  // int n = neighbors.size() - 3;

  // if (_gradientY.size() == 0) {
  _gradientY = gradient(neighbors, 1);
  // }
  return _gradientY;
}

template<class T>
LinearVector RBF<T>::gradientYAt(T point, std::vector<T> neighbors) {
  // int n = neighbors.size() - 3;

  // if (_gradientY.size() == 0) {
  _gradientY = gradientAt(point, neighbors, 1);
  // }
  return _gradientY;
}

template<class T>
LinearVector RBF<T>::gradientXYAt(T point, std::vector<T> neighbors) {
  // std::cerr << "Computeing gradient at " << std::endl;
  LinearMatrix M(neighbors.size() + 3);
  LinearVector b(neighbors.size() + 3);
  size_t n = neighbors.size();

  for (size_t i = 0; i < neighbors.size(); i++) {
    auto vdiff = point - neighbors.at(i);
    b[i] = (vdiff[0] * vdiff[1] / SQR(vdiff.length())) *
        (_kernel->d2phi(point, neighbors.at(i)) -
            _kernel->dphi(point, neighbors.at(i)));

    T p = neighbors[i];
    for (size_t d = 1; d <= p.dimension(); d++) {
      M(n + d, i) = p[d - 1];
      M(i, n + d) = p[d - 1];
    }
    M(n, i) = 1;
    M(i, n) = 1;
    for (size_t j = i; j < neighbors.size(); j++) {
      M(i, j) = _kernel->phi(neighbors.at(i), neighbors.at(j));
      M(j, i) = _kernel->phi(neighbors.at(i), neighbors.at(j));
    }
  }

  GaussJordanLinearSolver solver;
  LinearVector _g(neighbors.size() + 3);
  solver.solve(&M, &b, &_g);
  _g = b;
  return _g;
}

template<class T>
LinearVector RBF<T>::gradientZ(std::vector<T> neighbors) {
  // int n = neighbors.size() - 3;
  //
  // if( gradientZ.empty() ) {
  //     gradientZ = gradient( neighbors, 2 );
  // }
  // return gradientZ;
  UNUSED_VARIABLE(neighbors);
  return LinearVector();
}

template<class T>
LinearVector RBF<T>::divergent(std::vector<T> neighbors) {
  if (_divergent.size() == 0) {
    LinearVector gx = gradientX(neighbors);
    LinearVector gy = gradienty(neighbors);
    _divergent = gx + gy;
  }
  return _divergent;
}

template<class T>
LinearVector RBF<T>::divergentAt(T point, std::vector<T> neighbors) {
  // if (divergent.size() == 0) {
  // std::cerr << "Computing divergent " << std::endl;

  LinearVector gx = gradientAt(point, neighbors, 0);
  LinearVector gy = gradientAt(point, neighbors, 1);

  _divergent = gx + gy;
  // }
  return _divergent;
}

template<class T>
LinearVector RBF<T>::RBFWeights(std::vector<T> neighbors) {
  gradientX(neighbors);
  gradientY(neighbors);
  return laplacian(neighbors);
}
