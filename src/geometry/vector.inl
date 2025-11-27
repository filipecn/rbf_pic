template <typename T, size_t D> Vector<T, D>::Vector() {
  memset(v, 0, D * sizeof(T));
}

template <typename T, size_t D>
Vector<T, D>::Vector(std::initializer_list<T> values) : Vector() {
  int k = 0;
  for (auto value = values.begin(); value != values.end(); value++)
    v[k++] = *value;
}

template <typename T, size_t D>
Vector<T, D>::Vector(size_t n, const T *t) : Vector() {
  for (size_t i = 0; i < D && i < n; i++)
    v[i] = t[i];
}

template <typename T, size_t D> Vector<T, D>::Vector(const T &t) : Vector() {
  for (size_t i = 0; i < D; i++)
    v[i] = t;
}

template <typename T, size_t D>
Vector<T, D>::Vector(const T &x, const T &y) : Vector() {
  v[0] = x;
  v[1] = y;
}

template <typename T, size_t D>
Vector<T, D>::Vector(const T &x, const T &y, const T &z) : Vector() {
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

template <typename T, size_t D>
Vector<T, D>::Vector(const Vector<T, D> &other) {
  memcpy(v, other.v, D * sizeof(T));
}

template <typename T, size_t D> T Vector<T, D>::operator[](int i) const {
  ASSERT_FATAL(i >= 0 && i <= static_cast<int>(D));
  return v[i];
}

template <typename T, size_t D> T &Vector<T, D>::operator[](int i) {
  ASSERT_FATAL(i >= 0 && i <= static_cast<int>(D));
  return v[i];
}

template <typename T, size_t D>
bool Vector<T, D>::operator==(const Vector<T, D> &_v) const {
  for (size_t i = 0; i < D; i++)
    if (!IS_EQUAL(v[i], _v[i]))
      return false;
  return true;
}

template <typename T, size_t D>
bool Vector<T, D>::operator!=(const Vector<T, D> &_v) const {
  bool dif = false;
  for (size_t i = 0; i < D; i++)
    if (!IS_EQUAL(v[i], _v[i])) {
      dif = true;
      break;
    }
  return dif;
}

template <typename T, size_t D>
bool Vector<T, D>::operator<=(const Vector<T, D> &_v) const {
  for (size_t i = 0; i < D; i++)
    if (v[i] > _v[i])
      return false;
  return true;
}

template <typename T, size_t D>
bool Vector<T, D>::operator<(const Vector<T, D> &_v) const {
  for (size_t i = 0; i < D; i++)
    if (v[i] >= _v[i])
      return false;
  return true;
}

template <typename T, size_t D>
bool Vector<T, D>::operator>=(const Vector<T, D> &_v) const {
  for (size_t i = 0; i < D; i++)
    if (v[i] < _v[i])
      return false;
  return true;
}

template <typename T, size_t D>
bool Vector<T, D>::operator>(const Vector<T, D> &_v) const {
  for (size_t i = 0; i < D; i++)
    if (v[i] <= _v[i])
      return false;
  return true;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator-(const Vector<T, D> &_v) const {
  Vector<T, D> v_;
  for (size_t i = 0; i < D; i++)
    v_[i] = v[i] - _v[i];
  return v_;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator+(const Vector<T, D> &_v) const {
  Vector<T, D> v_;
  for (size_t i = 0; i < D; i++)
    v_[i] = v[i] + _v[i];
  return v_;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator*(const Vector<T, D> &_v) const {
  Vector<T, D> v_;
  for (size_t i = 0; i < D; i++)
    v_[i] = v[i] * _v[i];
  return v_;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator*(const T &f) const {
  Vector<T, D> v_;
  for (size_t i = 0; i < D; i++)
    v_[i] = v[i] * f;
  return v_;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator/(const Vector<T, D> &_v) const {
  Vector<T, D> v_;
  for (size_t i = 0; i < D; i++)
    v_[i] = v[i] / _v[i];
  return v_;
}

template <typename T, size_t D> Vector<T, D> Vector<T, D>::operator/=(T f) {
  for (size_t i = 0; i < D; i++)
    v[i] /= f;
  return *this;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator+=(const Vector<T, D> &_v) {
  for (size_t i = 0; i < D; i++)
    v[i] += _v[i];
  return *this;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator-=(const Vector<T, D> &_v) {
  for (size_t i = 0; i < D; i++)
    v[i] -= _v[i];
  return *this;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::operator/(T f) const {
  ASSERT_FATAL(IS_NOT_ZERO(f));
  T inv = static_cast<T>(1) / f;
  Vector<T, D> a;
  for (size_t i = 0; i < D; i++)
    a.v[i] = v[i] * inv;
  return a;
}

template <typename T, size_t D>
Vector<T, 2> Vector<T, D>::xy(size_t x, size_t y) const {
  return Vector<T, 2>(v[x], v[y]);
}

template <typename T, size_t D>
Vector<double, 2> Vector<T, D>::doubleXY(size_t x, size_t y) const {
  return Vector<double, 2>(static_cast<double>(v[x]),
                           static_cast<double>(v[y]));
}

template <typename T, size_t D>
Vector<double, 3> Vector<T, D>::doubleXYZ(size_t x, size_t y, size_t z) const {
  return Vector<double, 3>(static_cast<double>(v[x]), static_cast<double>(v[y]),
                           static_cast<double>(v[z]));
}

template <typename T, size_t D> T Vector<T, D>::max() const {
  T m = v[0];
  for (size_t i = 1; i < D; i++)
    m = std::max(m, v[i]);
  return m;
}

template <typename T, size_t D> double Vector<T, D>::length2() const {
  double sum = 0.f;
  for (size_t i = 0; i < D; i++)
    sum += SQR(v[i]);
  return sum;
}

template <typename T, size_t D> double Vector<T, D>::length() const {
  return std::sqrt(length2());
}

template <typename T, size_t D>
Vector<double, D> Vector<T, D>::normalized() const {
  double d = length();
  ASSERT_FATAL(IS_NOT_ZERO(d));
  Vector<double, D> r;
  for (size_t i = 0; i < D; i++)
    r[i] = v[i] / d;
  return r;
}

template <typename T, size_t D> Vector<T, 2> Vector<T, D>::right() const {
  return Vector<T, 2>(v[1], -v[0]);
}

template <typename T, size_t D> Vector<T, 2> Vector<T, D>::left() const {
  return Vector<T, 2>(-v[1], v[0]);
}

template <typename T, size_t D>
Vector<T, D> operator*(T f, const Vector<T, D> &v) {
  return v * f;
}

template <typename T> Vector<int, 3> round(const Vector<T, 3> &v) {
  return Vector<int, 3>(static_cast<int>(v[0] + (v[0] < 0 ? -1 : 1) * 0.5),
                        static_cast<int>(v[1] + (v[1] < 0 ? -1 : 1) * 0.5),
                        static_cast<int>(v[2] + (v[2] < 0 ? -1 : 1) * 0.5));
}

template <typename T> Vector<int, 2> round(const Vector<T, 2> &v) {
  return Vector<int, 2>(static_cast<int>(v[0] + (v[0] < 0 ? -1 : 1) * 0.5),
                        static_cast<int>(v[1] + (v[1] < 0 ? -1 : 1) * 0.5));
}

template <typename T> Vector<int, 3> ceil(const Vector<T, 3> &v) {
  return Vector<int, 3>(static_cast<int>(v[0] + 1.0),
                        static_cast<int>(v[1] + 1.0),
                        static_cast<int>(v[2] + 1.0));
}

template <typename T> Vector<int, 2> ceil(const Vector<T, 2> &v) {
  return Vector<int, 2>(static_cast<int>(v[0] + 1.0),
                        static_cast<int>(v[1] + 1.0));
}

template <typename T> Vector<int, 3> floor(const Vector<T, 3> &v) {
  return Vector<int, 3>(static_cast<int>(v[0]), static_cast<int>(v[1]),
                        static_cast<int>(v[2]));
}

template <typename T> Vector<int, 2> floor(const Vector<T, 2> &v) {
  return Vector<int, 2>(static_cast<int>(v[0]), static_cast<int>(v[1]));
}

template <typename T> Vector<T, 3> min(Vector<T, 3> a, Vector<T, 3> b) {
  return Vector<T, 3>(std::min(a[0], b[0]), std::min(a[1], b[1]),
                      std::min(a[2], b[2]));
}

template <typename T> Vector<T, 3> max(Vector<T, 3> a, Vector<T, 3> b) {
  return Vector<T, 3>(std::max(a[0], b[0]), std::max(a[1], b[1]),
                      std::max(a[2], b[2]));
}

template <typename T> Vector<T, 2> min(Vector<T, 2> a, Vector<T, 2> b) {
  return Vector<T, 2>(std::min(a[0], b[0]), std::min(a[1], b[1]));
}

template <typename T> Vector<T, 2> max(Vector<T, 2> a, Vector<T, 2> b) {
  return Vector<T, 2>(std::max(a[0], b[0]), std::max(a[1], b[1]));
}

template <typename T, size_t D> void normalize(Vector<T, D> &v) {
  double d = v.length();
  ASSERT_FATAL(IS_NOT_ZERO(d));
  for (size_t i = 0; i < D; i++)
    v[i] = v[i] / d;
}

template <typename T, size_t D>
double Vector<T, D>::dot(const Vector<T, D> &w) const {
  double sum = 0.;
  for (size_t i = 0; i < D; i++)
    sum += v[i] * w[i];
  return sum;
}

template <typename T, size_t D>
Vector<T, D> Vector<T, D>::projectInto(Vector<T, D> target) {
  double length = target.length();
  return target * (dot(target) / (length * length));
}
