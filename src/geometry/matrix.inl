template <typename T, size_t N> Matrix<T, N>::Matrix() {
  for (size_t i = 0; i < N; i++) {
    memset(m[i], 0, N * sizeof(T));
  }
}

template <typename T, size_t N>
Matrix<T, N>::Matrix(std::initializer_list<T> vs, bool cm) {
  size_t l = 0, c = 0;
  for (auto v = vs.begin(); v != vs.end(); v++) {
    m[l][c] = *v;
    if (!cm) {
      c++;
      if (c >= N) {
        c = 0, l++;
      }
    } else {
      l++;
      if (l >= N) {
        l = 0, c++;
      }
    }
  }
}

template <typename T, size_t N>
Matrix<T, N>::Matrix(std::initializer_list<Vector<T, N>> vs, bool columns) {
  size_t lc = 0;
  for (auto v : vs) {
    for (size_t i = 0; i < N; i++)
      if (columns)
        m[i][lc] = v[i];
      else
        m[lc][i] = v[i];
    lc++;
  }
}

template <typename T, size_t N>
Vector<T, N> Matrix<T, N>::operator*(const Vector<T, N> &v) const {
  Vector<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = 0.0;
    for (size_t j = 0; j < N; j++) {
      x[i] += v[j] * m[i][j];
    }
  }
  return x;
}

template <typename T, size_t N> T Matrix<T, N>::determinant() const {
  if (N == 3)
    return m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] +
           m[0][2] * m[1][0] * m[2][1] - m[2][0] * m[1][1] * m[0][2] -
           m[2][1] * m[1][2] * m[0][0] - m[2][2] * m[1][0] * m[0][1];
  NOT_IMPLEMENTED();
}

template <typename T, size_t N>
Matrix<T, N> Matrix<T, N>::operator*(const Matrix<T, N> &b) const {
  Matrix<T, N> c;
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      c.m[i][j] = 0.0;
      for (size_t k = 0; k < N; k++) {
        c.m[i][j] += m[i][k] * b.m[k][j];
      }
    }
  }
  return c;
}

template <typename T, size_t N>
T Matrix<T, N>::operator()(size_t i, size_t j) const {
  return m[i][j];
}

template <typename T, size_t N>
T &Matrix<T, N>::operator()(size_t i, size_t j) {
  return m[i][j];
}

template <typename T, size_t N> Matrix<T, N> Matrix<T, N>::transpose() const {
  Matrix<T, N> r;
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      r(j, i) = m[i][j];
    }
  }
  return r;
}

template <typename T, size_t N> void Matrix<T, N>::setIdentity() {
  for (size_t i = 0; i < N; i++) {
    memset(m[i], 0, N * sizeof(T));
    m[i][i] = 1.0;
  }
}

template <typename T, size_t N> bool Matrix<T, N>::isIdentity() const {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      if ((i == j) && !IS_EQUAL(m[i][j], 1.0)) {
        return false;
      } else if ((i != j) && IS_NOT_ZERO(m[i][j])) {
        return false;
      }
    }
  }
  return true;
}

template <typename T, size_t N> Matrix<T, N> inverse(const Matrix<T, N> &m) {
  Vector<T, 2 * N> em[N];
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      em[i][j] = m.m[i][j];
      if (i == j) {
        em[i][j + N] = 1.f;
      } else {
        em[i][j + N] = 0.f;
      }
    }
  }
  // Gauss-Jordan Elimination
  for (size_t k = 0; k < N; k++) {
    em[k] /= em[k][k];
    for (size_t i = 0; i < N; i++) {
      if (i != k) {
        em[i] -= em[k] * em[i][k];
      }
    }
  }
  Matrix<T, N> r;
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      r.m[i][j] = em[i][j + N];
    }
  }
  return r;
}
