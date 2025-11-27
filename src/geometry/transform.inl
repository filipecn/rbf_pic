template <int D> Transform<D>::Transform() {
  _m.setIdentity();
  _m_inv.setIdentity();
}

template <int D>
Transform<D>::Transform(Matrix<double, D + 1> m, Matrix<double, D + 1> m_inv)
    : _m(m), _m_inv(m_inv) {}

template <int D>
Point<double, D> Transform<D>::operator()(const Point<double, D> &p) const {
  Point<double, D + 1> tp(0.);
  for (int i = 0; i <= D; i++) {
    for (int j = 0; j < D; j++)
      tp[i] += p[j] * _m(i, j);
    tp[i] += _m(i, D);
  }
  Point<double, D> r;
  for (int i = 0; i < D; i++)
    r[i] = tp[i];
  if (tp[D] == 1.f)
    return r;
  return r / tp[D];
}

template <int D>
Vector<double, D> Transform<D>::operator()(const Vector<double, D> &v) const {
  Vector<double, D> tv(0.);
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      tv[i] += v[j] * _m(i, j);
  return tv;
}

template <int D>
Vector<double, D> Transform<D>::applyReflection(const Vector<double, D> &v,
                                                const Vector<double, D> &n) {
  return v - n * (2. * v.dot(n));
}

template <int D> Transform<D> Transform<D>::scale(Vector<double, D> s) {
  Matrix<double, D + 1> m;
  m.setIdentity();
  for (int i = 0; i < D; i++)
    m(i, i) = s[i];
  return Transform<D>(m, furoo::inverse(m));
}

template <int D> Transform<D> Transform<D>::translate(Vector<double, D> t) {
  Matrix<double, D + 1> m;
  m.setIdentity();
  for (int i = 0; i < D; i++)
    m(i, D) = t[i];
  return Transform<D>(m, furoo::inverse(m));
}

template <> inline Transform<2> Transform<2>::rotate(size_t axis, double angle) {
  Matrix<double, 3> m;
  m.setIdentity();
  double cos = std::cos(angle);
  double sin = std::sin(angle);
  m(0, 0) = cos;
  m(0, 1) = -sin;
  m(1, 0) = sin;
  m(1, 1) = cos;
  return Transform<2>(m, furoo::inverse(m));
}

template <> Transform<3> inline Transform<3>::rotate(size_t axis, double angle) {
  Matrix<double, 4> m;
  m.setIdentity();
  double cos = std::cos(angle);
  double sin = std::sin(angle);
  switch (axis) {
  case 0:
    m(1, 1) = cos;
    m(1, 2) = -sin;
    m(2, 1) = sin;
    m(2, 2) = cos;
    break;
  case 1:
    m(0, 0) = cos;
    m(2, 0) = -sin;
    m(0, 2) = sin;
    m(2, 2) = cos;
    break;
  case 2:
    m(0, 0) = cos;
    m(0, 1) = -sin;
    m(1, 0) = sin;
    m(1, 1) = cos;
    break;
  default:
    THROW(false, "Transform<>::rotate invalid axis");
  }
  return Transform<3>(m, furoo::inverse(m));
}

template <int D> Transform<D> Transform<D>::toUnitBox(BBox<double, D> b) {
  Vector<double, D> t;
  Vector<double, D> s;
  for (int d = 0; d < D; d++) {
    t[d] = -b.lower()[d];
    s[d] = 1. / b.size()[d];
  }
  auto m = scale(s).matrix() * translate(t).matrix();
  return Transform<D>(m, furoo::inverse(m));
}

template <int D> Matrix<double, D + 1> Transform<D>::inverseMatrix() const {
  return _m_inv;
}

template <int D> Matrix<double, D + 1> Transform<D>::matrix() const {
  return _m;
}

template <int D> Transform<D> Transform<D>::inverse() const {
  return Transform<D>(_m_inv, _m);
}

template <int D>
BBox<double, D> Transform<D>::operator()(const BBox<double, D> &b) const {
  return BBox<double, D>((*this)(b.lower()), (*this)(b.upper()));
}
