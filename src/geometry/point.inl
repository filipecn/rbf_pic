template<class T, int D>
Point<T, D>::Point(std::initializer_list<T> p) {
  int k = 0;
  for (auto it = p.begin(); it != p.end(); ++it) {
    if (k >= D) {
      break;
    }
    _v[k++] = *it;
  }
}

template<class T, int D>
Point<T, D>::Point() {
  for (int i = 0; i < D; i++) {
    _v[i] = static_cast<T>(0);
  }
}

template<class T, int D>
Point<T, D>::Point(T V) {
  for (int i = 0; i < D; i++) {
    _v[i] = static_cast<T>(V);
  }
}

template<class T, int D>
Point<T, D>::Point(T x_, T y_) {
  _v[0] = x_;
  _v[1] = y_;
}

template<class T, int D>
Point<T, D>::Point(T x_, T y_, T z_) {
  _v[0] = x_;
  _v[1] = y_;
  _v[2] = z_;
}

template<class T, int D>
bool operator==(const Point <T, D> &lhs, const Point <T, D> &rhs) {
  for (size_t i = 0; i < lhs.size; ++i) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

template<class T, int D>
bool Point<T, D>::operator==(const Point <T, D> &p) const {
  bool eq = true;

  for (int i = 0; i < D; i++) {
    eq = (eq && IS_EQUAL(_v[i], p[i]));
  }
  return eq;
}

template<class T, int D>
bool Point<T, D>::operator>=(const Point <T, D> &p) const {
  for (int i = 0; i < D; i++) {
    if (_v[i] < p[i]) {
      return false;
    }
  }
  return true;
}

template<class T, int D>
Point<T, D>::Point(const std::vector<T> &v) {
  for (size_t d = 0; d < D; d++)
    _v[d] = v[d];
}

template<class T, int D>
Point<T, D>::Point(const T *v) {
  for (size_t d = 0; d < D; d++)
    _v[d] = v[d];
}

// INLINE FUNCTIONS

template<class T, int D>
bool operator!=(const Point <T, D> &lhs, const Point <T, D> &rhs) {
  return !(lhs == rhs);
}

template<class T, int D>
double distance(const Point <T, D> &a, const Point <T, D> &b) {
  return (a - b).length();
}

template<class T, int D>
double distance2(const Point <T, D> &a, const Point <T, D> &b) {
  return (a - b).length2();
}
