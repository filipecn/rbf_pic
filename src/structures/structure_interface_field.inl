template <int D>
void StructureInterface<D>::FieldInterface::setName(const std::string &name) {
  _name = name;
}

template <int D>
std::string StructureInterface<D>::FieldInterface::name() const {
  return _name;
}

template <int D>
template <typename T>
StructureInterface<D>::Field<T>::Field(size_t size, T value,
                                       const std::string fieldName)
    : _defaultValue(value) {
  _data.resize(size, value);
  this->_name = fieldName;
}

template <int D>
template <typename T>
StructureInterface<D>::Field<T>::Field(size_t size, T value)
    : _defaultValue(value) {
  _data.resize(size, value);
}

template <int D>
template <typename T>
void StructureInterface<D>::Field<T>::setAll(T value) {
  std::fill(_data.begin(), _data.end(), value);
}

template <int D>
template <typename T>
T StructureInterface<D>::Field<T>::operator[](size_t i) const {
  return _data[i];
}

template <int D>
template <typename T>
T &StructureInterface<D>::Field<T>::operator[](size_t i) {
  THROW(i < _data.size(), "StructureInterface<>:operator[] invalid id");
  return _data[i];
}

template <int D>
template <typename T>
void StructureInterface<D>::Field<T>::clear() {
  _data.clear();
}

template <int D>
template <typename T>
void StructureInterface<D>::Field<T>::increaseSize(size_t size) {
  if (size == 1)
    _data.emplace_back(_defaultValue);
  else
    _data.resize(size, _defaultValue);
}

template <int D>
template <typename T>
void StructureInterface<D>::Field<T>::reset(int i) {
  if (i < 0)
    setAll(_defaultValue);
  else
    _data[i] = _defaultValue;
}

template <int D>
template <typename T>
std::string StructureInterface<D>::Field<T>::fieldDataTypeName() const {
  if /*constexpr*/ (std::is_same<T, double>::value)
    return "double";
  if /*constexpr*/ (std::is_same<T, int>::value)
    return "int";
  if /*constexpr*/ (std::is_same<T, Definitions::Material>::value)
    return "Definitions::Material";
  if /*constexpr*/ (std::is_same<T, unsigned char>::value)
    return "unsigned_char";
  if /*constexpr*/ (std::is_same<T, Definitions::Boundary>::value)
    return "Definitions::Boundary";
  if /*constexpr*/ (std::is_same<T, Point<double, D>>::value)
    return "Point<double,D>";
  return "void";
}
