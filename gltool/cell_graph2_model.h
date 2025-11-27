#ifndef FUROO_CELL_GRAPH2_MODEL_H
#define FUROO_CELL_GRAPH2_MODEL_H

#include "graphic_element_set.h"
#include <aergia/aergia.h>
#include <furoo.h>

class FieldModelInterface {
public:
  FieldModelInterface() = default;
  virtual ~FieldModelInterface() = default;
  virtual aergia::Color colorAt(size_t id);
  std::string name;
  size_t structureFieldId = 0;
  float alpha = .5f;
  double minValue = 0.f, maxValue = 1.f;
  bool isCategorical = false;
  bool visible = true;
};

template <typename T, int D> class FieldModel : public FieldModelInterface {
public:
  FieldModel() = default;
  explicit FieldModel(
      typename furoo::StructureInterface<D>::template Field<T> *f)
      : field(f) {}
  aergia::Color colorAt(size_t id) override {
    if (!this->visible)
      return aergia::COLOR_TRANSPARENT;
    if (!field)
      return FieldModelInterface::colorAt(id);
    if (isCategorical)
      return colorMap[(*field)[id]];
    return palette(1.f - ponos::linearStep(static_cast<float>((*field)[id]),
                                           minValue, maxValue),
                   alpha);
  }
  T valueAt(size_t id) {
    THROW(field, "FieldModel::valueAt null field pointer");
    return (*field)[id];
  }
  typename furoo::StructureInterface<D>::template Field<T> *field = nullptr;
  aergia::ColorPalette palette = aergia::HEAT_MATLAB_PALETTE;
  std::map<T, aergia::Color> colorMap;
};

template <int D> class PointFieldModel : public FieldModelInterface {
public:
  PointFieldModel() = default;
  explicit PointFieldModel(typename furoo::StructureInterface<
                           D>::template Field<furoo::Point<double, D>> *f)
      : field(f) {}
  furoo::Point<double, D> positionAt(size_t id) {
    THROW(field, "PointFieldModel::valueAt null field pointer");
    return (*field)[id];
  }
  typename furoo::StructureInterface<D>::template Field<furoo::Point<double, D>>
      *field = nullptr;
};

template <int D> class CellGraphTopologyModel : public aergia::SceneObject {
public:
  CellGraphTopologyModel() {}
  void set(const furoo::StructureInterface<D> *graph) {
    static bool first = true;
    if (first) {
      _edgesSet.reset(new GraphicElementSet(createArrow<D>()));
      if (D == 2) {
        _nodesSet.reset(new GraphicElementSet(ponos::create_quad_wireframe_mesh(
            ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
            ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0))));
      } else {
        static auto cubeWired = ponos::RawMeshes::cubeWireframe();
        static auto cube = ponos::RawMeshes::cube();
        _nodesSet = std::make_shared<GraphicElementSet>(cubeWired.get());
      }
      first = false;
    }
    this->structure = graph;
    updateEdges();
    updateNodes();
  }

  void updateEdges() {
    if (!structure)
      return;
    _edgesSet->resize(structure->faceCount());
    // update selected field min and max values
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _edgesSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _edgesSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateFaces([&](size_t faceId) {
      auto cells = structure->faceCells(faceId);
      auto pos = structure->cellCenterPosition(cells[0].id);
      auto m = _edgesSet->transformBuffer(bufferElementIndex);
      float size[D];
      for (int d = 0; d < D; d++)
        if (cells.size() == 2)
          size[d] = (structure->cellCenterPosition(cells[1].id) - pos)[d];
        else
          size[d] = 0.f;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size[0], size[1], (D == 3) ? size[2] : 0.f))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];

      aergia::Color color = aergia::COLOR_RED;
      auto c = _edgesSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _edgesSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }

  void updateNodes() {
    if (!structure)
      return;
    _nodesSet->resize(structure->cellCount());
    // unactivate dead particles, unfortunatelly we have to kill them all first
    // and make active ones live again
    {
      auto a = _nodesSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _nodesSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateCells([&](size_t cellId) {
      auto pos = structure->cellCenterPosition(cellId) -
                 structure->cellRegion(cellId).size() / 6.f;
      auto m = _nodesSet->transformBuffer(bufferElementIndex);
      float size = structure->cellRegion(cellId).size()[0] / 3;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_RED;
      auto c = _nodesSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _nodesSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      // if (useCutPlanes) {
      //   auto cellRegion = structure->cellRegion(cellId);
      //   a[0] = 0;
      //   for (int d = 0; d < D; d++) {
      //     if (useCut[d] && cutPlanes[d] >= cellRegion.lower()[d] &&
      //         cutPlanes[d] <= cellRegion.upper()[d])
      //       a[0] = 1;
      //   }
      // }
      bufferElementIndex++;
    });
  }
  // INTERFACE
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    if (!structure)
      return;
    _nodesSet->draw(camera, t);
    _edgesSet->draw(camera, t);
  }

private:
  std::shared_ptr<GraphicElementSet> _nodesSet;
  std::shared_ptr<GraphicElementSet> _edgesSet;
  const furoo::StructureInterface<D> *structure = nullptr;
};

template <int D> class CellGraphModel : public aergia::SceneObject {
public:
  CellGraphModel() {
    text.reset(new aergia::TextRenderer(FONT_PATH));
    //        "/mnt/windows/Projects/furoo/gltool/external/ponos/aergia/examples/assets/arial.ttf"));
    for (int d = 0; d < D; d++) {
      _selectedCellVectorFields[d] = -1;
      useCut[d] = false;
      cutPlanes[d] = 0.;
    }
  }
  void set(const furoo::StructureInterface<D> *graph) {
    static bool first = true;
    if (first) {
      _cellVectorSet.reset(new GraphicElementSet(createArrow<D>()));
      if (D == 2) {
        _cellsEdgesSet.reset(
            new GraphicElementSet(ponos::create_quad_wireframe_mesh(
                ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
                ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0))));
        _cellRegionsSet.reset(new GraphicElementSet(ponos::create_quad_mesh(
            ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
            ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0), false, false)));
        _facesEdgesSet.reset(
            new GraphicElementSet(ponos::create_quad_wireframe_mesh(
                ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
                ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0))));
        _faceRegionsSet.reset(new GraphicElementSet(ponos::create_quad_mesh(
            ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
            ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0), false, false)));
        _verticesEdgesSet.reset(
            new GraphicElementSet(ponos::create_quad_wireframe_mesh(
                ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
                ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0))));
        _vertexRegionsSet.reset(new GraphicElementSet(ponos::create_quad_mesh(
            ponos::Point3(0, 0, 0), ponos::Point3(1, 0, 0),
            ponos::Point3(1, 1, 0), ponos::Point3(0, 1, 0), false, false)));
      } else {
        static auto cubeWired = ponos::RawMeshes::cubeWireframe();
        static auto cube = ponos::RawMeshes::cube();
        _cellsEdgesSet = std::make_shared<GraphicElementSet>(cubeWired.get());
        _cellRegionsSet = std::make_shared<GraphicElementSet>(cube.get());
        _facesEdgesSet = std::make_shared<GraphicElementSet>(cubeWired.get());
        _faceRegionsSet = std::make_shared<GraphicElementSet>(cube.get());
        _verticesEdgesSet =
            std::make_shared<GraphicElementSet>(cubeWired.get());
        _vertexRegionsSet = std::make_shared<GraphicElementSet>(cube.get());
      }
      first = false;
    }
    this->structure = graph;
    updateCellsVector();
    updateCellsEdges();
    updateCellRegions();
    updateFaceRegions();
    updateFacesEdges();
    updateVertexRegions();
    updateVerticesEdges();
  }
  void updateCellsVector() {
    if (!structure)
      return;
    _cellVectorSet->resize(structure->cellCount());
    // update selected field min and max values
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _cellVectorSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _cellVectorSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateCells([&](size_t cellId) {
      auto pos = structure->cellCenterPosition(cellId);
      if (_selectedCellPointField >= 0)
        pos = cellPointFields[_selectedCellPointField]->positionAt(cellId);
      auto m = _cellVectorSet->transformBuffer(bufferElementIndex);
      float size[D];
      for (int d = 0; d < D; d++)
        if (_selectedCellVectorFields[d] >= 0)
          size[d] =
              cellVectorScale * dynamic_cast<FieldModel<double, D> *>(
                                    cellFields[_selectedCellVectorFields[d]])
                                    ->valueAt(cellId);
        else
          size[d] = 0.f;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size[0], size[1], (D == 3) ? size[2] : 0.f))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];

      aergia::Color color = aergia::COLOR_BLUE;
      auto c = _cellVectorSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _cellVectorSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }
  void updateCellsEdges() {
    if (!structure)
      return;
    _cellsEdgesSet->resize(structure->cellCount());
    // unactivate dead particles, unfortunatelly we have to kill them all first
    // and make active ones live again
    {
      auto a = _cellsEdgesSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _cellsEdgesSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateCells([&](size_t cellId) {
      auto pos = structure->cellCenterPosition(cellId) -
                 structure->cellRegion(cellId).size() / 2.f;
      auto m = _cellsEdgesSet->transformBuffer(bufferElementIndex);
      float size = structure->cellRegion(cellId).size()[0];
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_BLACK;
      auto c = _cellsEdgesSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _cellsEdgesSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      if (useCutPlanes) {
        auto cellRegion = structure->cellRegion(cellId);
        a[0] = 0;
        for (int d = 0; d < D; d++) {
          if (useCut[d] && cutPlanes[d] >= cellRegion.lower()[d] &&
              cutPlanes[d] <= cellRegion.upper()[d])
            a[0] = 1;
        }
      }
      bufferElementIndex++;
    });
  }
  void updateCellRegions() {
    if (!structure)
      return;
    _cellRegionsSet->resize(structure->cellCount());
    // update selected field min and max values
    if (_selectedCellField >= 0 &&
        !cellFields[_selectedCellField]->isCategorical) {
      cellFields[_selectedCellField]->minValue = structure->minScalarFieldValue(
          cellFields[_selectedCellField]->structureFieldId);
      cellFields[_selectedCellField]->maxValue = structure->maxScalarFieldValue(
          cellFields[_selectedCellField]->structureFieldId);
    }
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _cellRegionsSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _cellRegionsSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateCells([&](size_t cellId) {
      auto pos = structure->cellCenterPosition(cellId) -
                 structure->cellRegion(cellId).size() * 0.15;
      if (_selectedCellPointField >= 0)
        pos = cellPointFields[_selectedCellPointField]->positionAt(cellId);
      auto m = _cellRegionsSet->transformBuffer(bufferElementIndex);
      float size = structure->cellRegion(cellId).size()[0] * 0.3;
      if (_selectedCellPointField >= 0)
        size = 0.001;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];

      aergia::Color color = aergia::COLOR_TRANSPARENT;
      if (_selectedCellField >= 0) {
        color = cellFields[_selectedCellField]->colorAt(cellId);
        color.a = 1.f;
      }
      auto c = _cellRegionsSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _cellRegionsSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      if (useCutPlanes) {
        auto cellRegion = structure->cellRegion(cellId);
        a[0] = 0;
        for (int d = 0; d < D; d++) {
          if (useCut[d] && cutPlanes[d] >= cellRegion.lower()[d] &&
              cutPlanes[d] <= cellRegion.upper()[d])
            a[0] = 1;
        }
      }
      bufferElementIndex++;
    });
  }
  void updateFacesEdges() {
    if (!structure)
      return;
    _facesEdgesSet->resize(structure->faceCount());
    // unactivate dead faces, unfortunatelly we have to kill them all first
    // and make active ones live again
    {
      auto a = _facesEdgesSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _facesEdgesSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateFaces([&](size_t faceId) {
      auto pos = structure->faceCenterPosition(faceId) - faceSize / 2.f;
      auto m = _facesEdgesSet->transformBuffer(bufferElementIndex);

      float size = faceSize;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];

      aergia::Color color = aergia::COLOR_BLACK;
      auto c = _facesEdgesSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = 0.4;
      auto a = _facesEdgesSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }
  void updateFaceRegions() {
    if (!structure)
      return;
    _faceRegionsSet->resize(structure->faceCount());
    // update selected field min and max values
    if (_selectedFaceField >= 0 &&
        !faceFields[_selectedFaceField]->isCategorical) {
      faceFields[_selectedFaceField]->minValue = structure->minScalarFieldValue(
          faceFields[_selectedFaceField]->structureFieldId);
      faceFields[_selectedFaceField]->maxValue = structure->maxScalarFieldValue(
          faceFields[_selectedFaceField]->structureFieldId);
    }
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _faceRegionsSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _faceRegionsSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateFaces([&](size_t faceId) {
      auto pos = structure->faceCenterPosition(faceId) - faceSize / 2.f;
      auto m = _faceRegionsSet->transformBuffer(bufferElementIndex);
      float size = faceSize;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_TRANSPARENT;
      if (_selectedFaceField >= 0)
        color = faceFields[_selectedFaceField]->colorAt(faceId);
      auto c = _faceRegionsSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _faceRegionsSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }
  void updateVerticesEdges() {
    if (!structure)
      return;
    _verticesEdgesSet->resize(structure->vertexCount());
    // unactivate dead vertices, unfortunatelly we have to kill them all first
    // and make active ones live again
    {
      auto a = _verticesEdgesSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _verticesEdgesSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateVertices([&](size_t vertexId) {
      auto pos = structure->vertexPosition(vertexId) - vertexSize / 2.f;
      auto m = _verticesEdgesSet->transformBuffer(bufferElementIndex);
      float size = vertexSize;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_BLACK;
      auto c = _verticesEdgesSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _verticesEdgesSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }
  void updateVertexRegions() {
    if (!structure)
      return;
    _vertexRegionsSet->resize(structure->vertexCount());
    // update selected field min and max values
    if (_selectedVertexField >= 0 &&
        !vertexFields[_selectedVertexField]->isCategorical) {
      vertexFields[_selectedVertexField]->minValue =
          structure->minScalarFieldValue(
              cellFields[_selectedVertexField]->structureFieldId);
      vertexFields[_selectedVertexField]->maxValue =
          structure->maxScalarFieldValue(
              cellFields[_selectedVertexField]->structureFieldId);
    }
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _vertexRegionsSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _vertexRegionsSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    structure->iterateVertices([&](size_t vertexId) {
      auto pos = structure->vertexPosition(vertexId) - vertexSize / 2.f;
      auto m = _vertexRegionsSet->transformBuffer(bufferElementIndex);
      float size = vertexSize;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_TRANSPARENT;
      if (_selectedVertexField >= 0)
        color = vertexFields[_selectedVertexField]->colorAt(vertexId);
      auto c = _vertexRegionsSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _vertexRegionsSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    });
  }
  size_t addField(FieldModelInterface *field,
                  furoo::Definitions::MeshLocation fieldLocation,
                  bool spatial = false) {
    if (!spatial)
      switch (fieldLocation) {
      case furoo::Definitions::MeshLocation::CELL_CENTER:
        cellFields.push_back(field);
        return cellFields.size() - 1;
      case furoo::Definitions::MeshLocation::FACE_CENTER:
      case furoo::Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      case furoo::Definitions::MeshLocation::VERTICAL_FACE_CENTER:
        faceFields.push_back(field);
        return faceFields.size() - 1;
      case furoo::Definitions::MeshLocation::VERTEX_CENTER:
        vertexFields.push_back(field);
        return vertexFields.size() - 1;
      default:
        return 0;
      }
    else
      switch (fieldLocation) {
      case furoo::Definitions::MeshLocation::CELL_CENTER:
        cellPointFields.push_back(dynamic_cast<PointFieldModel<D> *>(field));
        return cellPointFields.size() - 1;
      case furoo::Definitions::MeshLocation::FACE_CENTER:
      case furoo::Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      case furoo::Definitions::MeshLocation::VERTICAL_FACE_CENTER:
        facePointFields.push_back(dynamic_cast<PointFieldModel<D> *>(field));
        return facePointFields.size() - 1;
      case furoo::Definitions::MeshLocation::VERTEX_CENTER:
        vertexPointFields.push_back(dynamic_cast<PointFieldModel<D> *>(field));
        return vertexPointFields.size() - 1;
      default:
        return 0;
      }
    return 0;
  }
  size_t addCellField(FieldModelInterface *field) {
    cellFields.push_back(field);
    return cellFields.size() - 1;
  }
  size_t addFaceField(FieldModelInterface *field) {
    faceFields.push_back(field);
    return faceFields.size() - 1;
  }
  size_t addVertexField(FieldModelInterface *field) {
    vertexFields.push_back(field);
    return vertexFields.size() - 1;
  }
  void selectCellVectorField(int i, int d) { _selectedCellVectorFields[d] = i; }
  void selectCellField(int i) { _selectedCellField = i; }
  void selectFaceField(int i) { _selectedFaceField = i; }
  void selectVertexField(int i) { _selectedVertexField = i; }
  void selectCellPointField(int i) { _selectedCellPointField = i; }
  double selectedCellFieldMinValue() const {
    if (_selectedCellField >= 0)
      return cellFields[_selectedCellField]->minValue;
    return 0;
  }
  double selectedCellFieldMaxValue() const {
    if (_selectedCellField >= 0)
      return cellFields[_selectedCellField]->maxValue;
    return 0;
  }
  double selectedFaceFieldMinValue() const {
    if (_selectedFaceField >= 0)
      return faceFields[_selectedFaceField]->minValue;
    return 0;
  }
  double selectedFaceFieldMaxValue() const {
    if (_selectedFaceField >= 0)
      return faceFields[_selectedFaceField]->maxValue;
    return 0;
  }
  double selectedVertexFieldMinValue() const {
    if (_selectedVertexField >= 0)
      return vertexFields[_selectedVertexField]->minValue;
    return 0;
  }
  double selectedVertexFieldMaxValue() const {
    if (_selectedVertexField >= 0)
      return vertexFields[_selectedVertexField]->maxValue;
    return 0;
  }
  void clearFields() {
    _selectedVertexField = -1;
    _selectedFaceField = -1;
    _selectedCellField = -1;
    cellFields.clear();
    faceFields.clear();
    vertexFields.clear();
    cellPointFields.clear();
    facePointFields.clear();
    vertexPointFields.clear();
  }
  // INTERFACE
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    if (!structure)
      return;
    // glDisable(GL_DEPTH_TEST);
    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    _cellsEdgesSet->draw(camera, t);
    if (showCellRegions)
      _cellRegionsSet->draw(camera, t);
    //  _facesEdgesSet->draw(camera, t);
    if (showFaceRegions)
      _faceRegionsSet->draw(camera, t);
    //  _verticesEdgesSet->draw(camera, t);
    if (showVertexRegions)
      _vertexRegionsSet->draw(camera, t);
    if (showCellVectors) {
      glLineWidth(4.);
      _cellVectorSet->draw(camera, t);
      glLineWidth(1.);
    }
    // render cell ids
    text->setCamera(camera);
    text->textSize = textSize;
    if (showCellText)
      structure->iterateCells([&](size_t i) {
        std::ostringstream s;
        if (cellFields.size() && _selectedCellField >= 0 &&
            !cellFields[_selectedCellField]->isCategorical) {
          double v = (*dynamic_cast<FieldModel<double, D> *>(
                           cellFields[_selectedCellField])
                           ->field)[i];
          s << v;
        } else
          s << "c" << i;
        auto pos = structure->cellCenterPosition(i);
        text->at(ponos::Point3(pos.x(), pos.y(), (D == 2) ? 1.001 : pos.z()))
            << s.str();
      });
    if (showFaceText)
      structure->iterateFaces([&](size_t i) {
        std::vector<furoo::Definitions::Neighbor> edge =
            structure->faceCells(i);
        std::ostringstream s;
        if (faceFields.size() && _selectedFaceField >= 0 &&
            !faceFields[_selectedFaceField]->isCategorical) {
          double v = (*dynamic_cast<FieldModel<double, D> *>(
                           faceFields[_selectedFaceField])
                           ->field)[i];
          if (std::fabs(v) < 1e-15)
            s << "0";
          else
            s << v;
        } else
          s << "f" << i << "(" << edge[0].id << ","
            << ((edge.size() > 1) ? edge[1].id : -1) << ")";
        auto pos = structure->faceCenterPosition(i);
        text->at(ponos::Point3(pos.x(), pos.y(), (D == 2) ? 1.001 : pos.z()))
            << s.str();
      });
    if (showVertexText)
      structure->iterateVertices([&](size_t i) {
        std::ostringstream s;
        if (vertexFields.size() && _selectedVertexField >= 0 &&
            !vertexFields[_selectedVertexField]->isCategorical) {
          s << "NOT FOUND";
        } else
          s << "v" << i;
        auto pos = structure->vertexPosition(i);
        text->at(ponos::Point3(pos.x(), pos.y(), (D == 2) ? 1.001 : pos.z()))
            << s.str();
      });
    if (showText) {
      std::ostringstream s;
      s << structure->faceCount() << " faces " << structure->cellCount()
        << " cells ";
      s << structure->vertexCount() << " vertices.";
      text->render(s.str(), 0, 0);
    }
  }
  bool intersect(const ponos::Ray3 &r, float *t) override { return false; }

  bool useCutPlanes = false;
  double cutPlanes[D];
  bool useCut[D];
  float textSize = 0.0003;
  bool showText = false;
  bool showFaceText = false;
  bool showCellText = false;
  bool showVertexText = false;
  bool showFaceRegions = true;
  bool showCellRegions = true;
  bool showVertexRegions = true;
  bool showCellVectors = false;
  float vertexSize = 0.01;
  float faceSize = 0.01;
  float cellVectorScale = 1.f;
  aergia::Color nodeColor = aergia::Color(0, 0, 0, 0.4);
  aergia::Color edgeColor = aergia::Color(0, 0, 0, 0.4);
  aergia::Color regionColor = aergia::Color(1, 0, 0, 0.2);
  std::vector<size_t> selectedCells;
  std::vector<FieldModelInterface *> cellFields;
  std::vector<FieldModelInterface *> vertexFields;
  std::vector<FieldModelInterface *> faceFields;
  std::vector<PointFieldModel<D> *> cellPointFields;
  std::vector<PointFieldModel<D> *> vertexPointFields;
  std::vector<PointFieldModel<D> *> facePointFields;

protected:
  std::shared_ptr<aergia::TextRenderer> text;
  // instances fields
  std::shared_ptr<GraphicElementSet> _cellsEdgesSet;
  std::shared_ptr<GraphicElementSet> _cellRegionsSet;
  std::shared_ptr<GraphicElementSet> _cellVectorSet;
  std::shared_ptr<GraphicElementSet> _facesEdgesSet;
  std::shared_ptr<GraphicElementSet> _faceRegionsSet;
  std::shared_ptr<GraphicElementSet> _verticesEdgesSet;
  std::shared_ptr<GraphicElementSet> _vertexRegionsSet;
  int _selectedCellVectorFields[D];
  int _selectedCellField = -1;
  int _selectedFaceField = -1;
  int _selectedVertexField = -1;
  int _selectedCellPointField = -1;

  const furoo::StructureInterface<D> *structure = nullptr;
};

// TODO: deprecated
class CellGraph2Model : public aergia::SceneObject {
public:
  CellGraph2Model();
  explicit CellGraph2Model(const furoo::CellGraph2 *graph);
  void set(const furoo::CellGraph2 *graph);
  void updateCellsEdges();
  void updateCellRegions();
  void updateFacesEdges();
  void updateFaceRegions();
  void updateVerticesEdges();
  void updateVertexRegions();
  void selectCellField(int id);
  void selectFaceField(int id);
  void selectVertexField(int id);
  size_t addCellField(FieldModelInterface *field);
  size_t addFaceField(FieldModelInterface *field);
  size_t addVertexField(FieldModelInterface *field);
  double selectedCellFieldMinValue() const;
  double selectedCellFieldMaxValue() const;
  double selectedFaceFieldMinValue() const;
  double selectedFaceFieldMaxValue() const;
  double selectedVertexFieldMinValue() const;
  double selectedVertexFieldMaxValue() const;
  // INTERFACE
  void draw(const aergia::CameraInterface *camera, ponos::Transform t) override;
  bool intersect(const ponos::Ray3 &r, float *t) override;

  float textSize = 0.0003;
  bool showText = false;
  bool showFaceText = false;
  bool showCellText = false;
  bool showVertexText = false;
  float vertexSize = 0.01;
  float faceSize = 0.01;
  aergia::Color nodeColor = aergia::Color(0, 0, 0, 0.4);
  aergia::Color edgeColor = aergia::Color(0, 0, 0, 0.4);
  aergia::Color regionColor = aergia::Color(1, 0, 0, 0.2);
  const furoo::CellGraph2 *graph = nullptr;
  std::vector<size_t> selectedCells;
  std::vector<FieldModelInterface *> cellFields;
  std::vector<FieldModelInterface *> vertexFields;
  std::vector<FieldModelInterface *> faceFields;

private:
  std::shared_ptr<aergia::TextRenderer> text;
  // instances fields
  std::shared_ptr<GraphicElementSet> _cellsEdgesSet;
  std::shared_ptr<GraphicElementSet> _cellRegionsSet;
  std::shared_ptr<GraphicElementSet> _facesEdgesSet;
  std::shared_ptr<GraphicElementSet> _faceRegionsSet;
  std::shared_ptr<GraphicElementSet> _verticesEdgesSet;
  std::shared_ptr<GraphicElementSet> _vertexRegionsSet;
  int _selectedCellField = -1;
  int _selectedFaceField = -1;
  int _selectedVertexField = -1;
};

#endif // FUROO_CELL_GRAPH2_MODEL_H
