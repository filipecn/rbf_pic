template<int D>
StructureInterface<D>::StructureInterface() = default;

template<int D>
StructureInterface<D>::~StructureInterface() = default;

template<int D>
template<typename T>
void StructureInterface<D>::changeFieldType(size_t fieldId, T defaultValue) {
  switch (_fieldsLocations[fieldId]) {
    case Definitions::MeshLocation::CELL_CENTER:
      _fields[fieldId].reset(new Field<T>(allCellDataCount(), defaultValue));
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      _fields[fieldId].reset(new Field<T>(allVertexDataCount(), defaultValue));
      break;
    case Definitions::MeshLocation::FACE_CENTER:
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      _fields[fieldId].reset(new Field<T>(allFaceDataCount(), defaultValue));
      break;
    default:break;
  }
}

template<int D>
template<typename T>
size_t StructureInterface<D>::addCellField(T value, std::string fieldName) {
  _fieldsLocations.emplace_back(Definitions::MeshLocation::CELL_CENTER);
  _fields.emplace_back(new Field<T>(allCellDataCount() + 1, value, fieldName));
  _cellFieldIds.emplace_back(_fields.size() - 1u);
  return _fields.size() - 1u;
}

template<int D>
template<typename T>
size_t StructureInterface<D>::addFaceField(T value, std::string fieldName,
                                           Definitions::MeshLocation location) {
  _fieldsLocations.emplace_back(location);
  _fields.emplace_back(new Field<T>(allFaceDataCount(), value, fieldName));
  _faceFieldIds.emplace_back(_fields.size() - 1u);
  return _fields.size() - 1u;
}
template<int D>
template<typename T>
size_t StructureInterface<D>::addVertexField(T value, std::string fieldName) {
  _fieldsLocations.emplace_back(Definitions::MeshLocation::VERTEX_CENTER);
  _fields.emplace_back(new Field<T>(allVertexDataCount(), value, fieldName));
  _vertexFieldIds.emplace_back(_fields.size() - 1u);
  return _fields.size() - 1u;
}

template<int D>
template<typename T>
size_t StructureInterface<D>::addField(Definitions::MeshLocation location,
                                       T value, std::string fieldName) {
  switch (location) {
    case Definitions::MeshLocation::CELL_CENTER:
      _fields.emplace_back(new Field<T>(allCellDataCount(),
                                        value, fieldName));
      _cellFieldIds.emplace_back(_fields.size() - 1u);
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      _fields.emplace_back(new Field<T>(allVertexDataCount(),
                                        value,
                                        fieldName));
      _vertexFieldIds.emplace_back(_fields.size() - 1u);
      break;
    case Definitions::MeshLocation::FACE_CENTER:
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      _fields.emplace_back(new Field<T>(allFaceDataCount(), value, fieldName));
      _faceFieldIds.emplace_back(_fields.size() - 1u);
      break;
    default:break;
  }
  _fieldsLocations.emplace_back(location);
  return _fields.size() - 1u;
}

template<int D>
void StructureInterface<D>::removeFields() {
  _fieldsLocations.clear();
  _cellFieldIds.clear();
  _faceFieldIds.clear();
  _vertexFieldIds.clear();
  _fields.clear();
}

template<int D>
std::string StructureInterface<D>::fieldName(size_t fieldId) const {
  return _fields[fieldId]->name();
}

template<int D>
template<typename T>
typename StructureInterface<D>::template Field<T> *StructureInterface<D>::field(
    size_t fieldId) {
  return dynamic_cast<Field <T> *>(_fields[fieldId].get());
}

template<int D>
template<typename T>
const typename StructureInterface<D>::template Field<T> *StructureInterface<D>::field(
    size_t fieldId) const {
  return dynamic_cast<const Field <T> *>(_fields[fieldId].get());
}

template<int D>
Definitions::MeshLocation StructureInterface<D>::fieldLocation(size_t fieldId) const {
  THROW(fieldId < _fieldsLocations.size(),
        "StructureInterface::fieldLocation invalid field id.")
  return _fieldsLocations[fieldId];
}

template<int D>
size_t StructureInterface<D>::fieldCount() const {
  return _fields.size();
}

template<int D>
std::string StructureInterface<D>::fieldDataTypeName(size_t fieldId) const {
  return _fields[fieldId]->fieldDataTypeName();
}

template<int D>
double StructureInterface<D>::maxScalarFieldValue(size_t fieldId) const {
  THROW(fieldId < this->_fields.size(),
        "StructureInterfacae<>::maxScalarFieldValue invalid field id");
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  double value = -INFINITY;
  switch (this->_fieldsLocations[fieldId]) {
    case Definitions::MeshLocation::CELL_CENTER:
      iterateCells([&](size_t i) {
        value = std::max(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::FACE_CENTER:
      iterateFaces([&](size_t i) {
        value = std::max(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::HORIZONTAL)
          value = std::max(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::DEPTH)
          value = std::max(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::VERTICAL)
          value = std::max(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      iterateVertices([&](size_t i) {
        value = std::max(value, field[i]);
      });
      break;
    default: THROW(false,
                   "StructureInterface<>::maxScalarFieldValue invalid mesh location")
  }
  return value;
}

template<int D>
double StructureInterface<D>::minScalarFieldValue(size_t fieldId) const {
  THROW(fieldId < this->_fields.size(),
        "StructureInterface<>::minScalarFieldValue invalid field id");
  auto &field = *dynamic_cast<Field<double> *>(this->_fields[fieldId].get());
  double value = INFINITY;
  switch (this->_fieldsLocations[fieldId]) {
    case Definitions::MeshLocation::CELL_CENTER:
      iterateCells([&](size_t i) {
        value = std::min(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::FACE_CENTER:
      iterateFaces([&](size_t i) {
        value = std::min(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::HORIZONTAL)
          value = std::min(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::VERTICAL)
          value = std::min(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      iterateFaces([&](size_t i) {
        if (faceOrientation(i) == Definitions::Orientation::DEPTH)
          value = std::min(value, field[i]);
      });
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      iterateVertices([&](size_t i) {
        value = std::min(value, field[i]);
      });
      break;
    default: THROW(false,
                   "StructureInterface<>::minScalarFieldValue invalid mesh location")
  }
  return value;
}
