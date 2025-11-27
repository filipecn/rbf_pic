#include "cell_graph2_model.h"
#include "ponos_wrapper.h"

CellGraph2Model::CellGraph2Model() {
  text.reset(new aergia::TextRenderer(
      "/mnt/windows/Projects/furoo/gltool/external/ponos/aergia/examples/assets/arial.ttf"));
  _cellsEdgesSet.reset(new GraphicElementSet(ponos::create_quad_wireframe_mesh(
      ponos::Point3(0, 0, 0),
      ponos::Point3(1, 0, 0),
      ponos::Point3(1, 1, 0),
      ponos::Point3(0, 1, 0))));
  _cellRegionsSet.reset(new GraphicElementSet(
      ponos::create_quad_mesh(
          ponos::Point3(0, 0, 0),
          ponos::Point3(1, 0, 0),
          ponos::Point3(1, 1, 0),
          ponos::Point3(0, 1, 0),
          false,
          false)));
  _facesEdgesSet.reset(new GraphicElementSet(ponos::create_quad_wireframe_mesh(
      ponos::Point3(0, 0, 0),
      ponos::Point3(1, 0, 0),
      ponos::Point3(1, 1, 0),
      ponos::Point3(0, 1, 0))));
  _faceRegionsSet.reset(new GraphicElementSet(
      ponos::create_quad_mesh(
          ponos::Point3(0, 0, 0),
          ponos::Point3(1, 0, 0),
          ponos::Point3(1, 1, 0),
          ponos::Point3(0, 1, 0),
          false,
          false)));
  _verticesEdgesSet.reset(new GraphicElementSet(ponos::create_quad_wireframe_mesh(
      ponos::Point3(0, 0, 0),
      ponos::Point3(1, 0, 0),
      ponos::Point3(1, 1, 0),
      ponos::Point3(0, 1, 0))));
  _vertexRegionsSet.reset(new GraphicElementSet(
      ponos::create_quad_mesh(
          ponos::Point3(0, 0, 0),
          ponos::Point3(1, 0, 0),
          ponos::Point3(1, 1, 0),
          ponos::Point3(0, 1, 0),
          false,
          false)));
}

CellGraph2Model::CellGraph2Model(const furoo::CellGraph2 *graph)
    : CellGraph2Model() {
  set(graph);
}

void CellGraph2Model::set(const furoo::CellGraph2 *graph) {
  this->graph = graph;
  updateCellsEdges();
  updateCellRegions();
  updateFaceRegions();
  updateFacesEdges();
  updateVertexRegions();
  updateVerticesEdges();
}

void CellGraph2Model::draw(const aergia::CameraInterface *camera,
                           ponos::Transform t) {
  if (!graph)
    return;
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  _cellsEdgesSet->draw(camera, t);
  _cellRegionsSet->draw(camera, t);
//  _facesEdgesSet->draw(camera, t);
  _faceRegionsSet->draw(camera, t);
//  _verticesEdgesSet->draw(camera, t);
  _vertexRegionsSet->draw(camera, t);
  // render cell ids
  text->setCamera(camera);
  text->textSize = textSize;
  if (showCellText && _selectedCellField >= 0
      && !cellFields[_selectedCellField]->isCategorical) {
    graph->iterateCells([&](size_t i) {
      double v = (*dynamic_cast<FieldModel<double, 2> *>(
          cellFields[_selectedCellField])->field)[i];
      std::ostringstream s;
      s << v;
      auto pos = graph->cellCenterPosition(i);
      text->at(ponos::Point3(pos.x(), pos.y(), graph->node(i).level()))
          << s.str();
    });
  }
  if (showFaceText && _selectedFaceField >= 0
      && !faceFields[_selectedFaceField]->isCategorical) {
    graph->iterateFaces([&](size_t i) {
      double v = (*dynamic_cast<FieldModel<double, 2> *>(
          faceFields[_selectedFaceField])->field)[i];
      std::ostringstream s;
      if (std::fabs(v) < 1e-15)
        s << "0";
      else s << v;
      auto pos = graph->faceCenterPosition(i);
      text->at(ponos::Point3(pos.x(), pos.y(), 1)) << s.str();
    });
  }
  if (showVertexText && _selectedVertexField >= 0
      && !vertexFields[_selectedVertexField]->isCategorical) {
    graph->iterateVertices([&](size_t i) {
      std::ostringstream s;
      s << "v" << i;
      auto pos = graph->vertexPosition(i);
      text->at(ponos::Point3(pos.x(), pos.y(), 1)) << s.str();
    });
  }
  if (showText) {
    std::ostringstream s;
    s << graph->faceCount() << " faces " << graph->cellCount() << " cells ";
    s << graph->vertexCount() << " vertices.";
    text->render(s.str(), 0, 0);
  }
}

bool CellGraph2Model::intersect(const ponos::Ray3 &r, float *t) {
  return SceneObject::intersect(r, t);
}

void CellGraph2Model::updateCellsEdges() {
  if (!graph)
    return;
  _cellsEdgesSet->resize(graph->cellCount());
  // unactivate dead particles, unfortunatelly we have to kill them all first
  // and make active ones live again
  {
    auto a = _cellsEdgesSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _cellsEdgesSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateCells([&](size_t cellId) {
    auto pos =
        graph->cellCenterPosition(cellId)
            - graph->cellRegion(cellId).size() / 2.f;
    auto p = _cellsEdgesSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0];
    p[1] = pos[1];
    p[2] = 0.f;//pos[2];
    p[3] = graph->cellRegion(cellId).size()[0];
    aergia::Color color = aergia::COLOR_BLACK;
    auto c = _cellsEdgesSet->colorBuffer(bufferElementIndex);
    c[0] = color.r;
    c[1] = color.g;
    c[2] = color.b;
    c[3] = color.a;
    auto a = _cellsEdgesSet->activeBuffer(bufferElementIndex);
    a[0] = 1;
    bufferElementIndex++;
  });
}

void CellGraph2Model::updateCellRegions() {
  if (!graph)
    return;
  _cellRegionsSet->resize(graph->cellCount());
  // update selected field min and max values
  if (_selectedCellField >= 0
      && !cellFields[_selectedCellField]->isCategorical) {
    cellFields[_selectedCellField]->minValue =
        graph->minScalarFieldValue(cellFields[_selectedCellField]->structureFieldId);
    cellFields[_selectedCellField]->maxValue =
        graph->maxScalarFieldValue(cellFields[_selectedCellField]->structureFieldId);
  }
  {
    // unactivate dead cells, unfortunatelly we have to kill them all first
    // and make active ones live again
    auto a = _cellRegionsSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _cellRegionsSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateCells([&](size_t cellId) {
    auto pos =
        graph->cellCenterPosition(cellId)
            - graph->cellRegion(cellId).size() / 2.f;
    auto p = _cellRegionsSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0]; // x
    p[1] = pos[1]; // y
    p[2] = 0.f;    // z
    p[3] = graph->cellRegion(cellId).size()[0];
    aergia::Color color = aergia::COLOR_TRANSPARENT;
    if (_selectedCellField >= 0)
      color = cellFields[_selectedCellField]->colorAt(cellId);
    auto c = _cellRegionsSet->colorBuffer(bufferElementIndex);
    c[0] = color.r;
    c[1] = color.g;
    c[2] = color.b;
    c[3] = color.a;
    auto a = _cellRegionsSet->activeBuffer(bufferElementIndex);
    a[0] = 1;
    bufferElementIndex++;
  });
}

void CellGraph2Model::updateFacesEdges() {
  if (!graph)
    return;
  _facesEdgesSet->resize(graph->faceCount());
  // unactivate dead faces, unfortunatelly we have to kill them all first
  // and make active ones live again
  {
    auto a = _facesEdgesSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _facesEdgesSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateFaces([&](size_t faceId) {
    auto pos =
        graph->faceCenterPosition(faceId) - faceSize / 2.f;
    auto p = _facesEdgesSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0];
    p[1] = pos[1];
    p[2] = 0.f;//pos[2];
    p[3] = faceSize;
    aergia::Color color = aergia::COLOR_BLACK;
    auto c = _facesEdgesSet->colorBuffer(bufferElementIndex);
    c[0] = color.r;
    c[1] = color.g;
    c[2] = color.b;
    c[3] = color.a;
    auto a = _facesEdgesSet->activeBuffer(bufferElementIndex);
    a[0] = 1;
    bufferElementIndex++;
  });
}

void CellGraph2Model::updateFaceRegions() {
  if (!graph)
    return;
  _faceRegionsSet->resize(graph->faceCount());
  // update selected field min and max values
  if (_selectedFaceField >= 0
      && !faceFields[_selectedFaceField]->isCategorical) {
    faceFields[_selectedFaceField]->minValue =
        graph->minScalarFieldValue(faceFields[_selectedFaceField]->structureFieldId);
    faceFields[_selectedFaceField]->maxValue =
        graph->maxScalarFieldValue(faceFields[_selectedFaceField]->structureFieldId);
  }
  {
    // unactivate dead cells, unfortunatelly we have to kill them all first
    // and make active ones live again
    auto a = _faceRegionsSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _faceRegionsSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateFaces([&](size_t faceId) {
    auto pos =
        graph->faceCenterPosition(faceId) - faceSize / 2.f;
    auto p = _faceRegionsSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0]; // x
    p[1] = pos[1]; // y
    p[2] = 0.f;    // z
    p[3] = faceSize;
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

void CellGraph2Model::updateVerticesEdges() {
  if (!graph)
    return;
  _verticesEdgesSet->resize(graph->vertexCount());
  // unactivate dead vertices, unfortunatelly we have to kill them all first
  // and make active ones live again
  {
    auto a = _verticesEdgesSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _verticesEdgesSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateVertices([&](size_t vertexId) {
    auto pos =
        graph->vertexPosition(vertexId) - vertexSize / 2.f;
    auto p = _verticesEdgesSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0];
    p[1] = pos[1];
    p[2] = 0.f;//pos[2];
    p[3] = vertexSize;
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

void CellGraph2Model::updateVertexRegions() {
  if (!graph)
    return;
  _vertexRegionsSet->resize(graph->vertexCount());
  // update selected field min and max values
  if (_selectedVertexField >= 0
      && !vertexFields[_selectedVertexField]->isCategorical) {
    vertexFields[_selectedVertexField]->minValue =
        graph->minScalarFieldValue(cellFields[_selectedVertexField]->structureFieldId);
    vertexFields[_selectedVertexField]->maxValue =
        graph->maxScalarFieldValue(cellFields[_selectedVertexField]->structureFieldId);
  }
  {
    // unactivate dead cells, unfortunatelly we have to kill them all first
    // and make active ones live again
    auto a = _vertexRegionsSet->activeBuffer();
    memset(a, 0, sizeof(unsigned int) * _vertexRegionsSet->count());
  }
  // update positions
  size_t bufferElementIndex = 0;
  graph->iterateVertices([&](size_t vertexId) {
    auto pos =
        graph->vertexPosition(vertexId) - vertexSize / 2.f;
    auto p = _vertexRegionsSet->transformBuffer(bufferElementIndex);
    p[0] = pos[0]; // x
    p[1] = pos[1]; // y
    p[2] = 0.f;    // z
    p[3] = vertexSize;
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

size_t CellGraph2Model::addCellField(FieldModelInterface *field) {
  cellFields.push_back(field);
  return cellFields.size() - 1;
}

size_t CellGraph2Model::addFaceField(FieldModelInterface *field) {
  faceFields.push_back(field);
  return faceFields.size() - 1;
}

size_t CellGraph2Model::addVertexField(FieldModelInterface *field) {
  vertexFields.push_back(field);
  return vertexFields.size() - 1;
}

void CellGraph2Model::selectCellField(int i) {
  _selectedCellField = i;
}

void CellGraph2Model::selectFaceField(int i) {
  _selectedFaceField = i;
}

void CellGraph2Model::selectVertexField(int i) {
  _selectedVertexField = i;
}

double CellGraph2Model::selectedCellFieldMinValue() const {
  if (_selectedCellField >= 0)
    return cellFields[_selectedCellField]->minValue;
  return 0;
}

double CellGraph2Model::selectedCellFieldMaxValue() const {
  if (_selectedCellField >= 0)
    return cellFields[_selectedCellField]->maxValue;
  return 0;
}

double CellGraph2Model::selectedFaceFieldMinValue() const {
  if (_selectedFaceField >= 0)
    return faceFields[_selectedFaceField]->minValue;
  return 0;
}

double CellGraph2Model::selectedFaceFieldMaxValue() const {
  if (_selectedFaceField >= 0)
    return faceFields[_selectedFaceField]->maxValue;
  return 0;
}

double CellGraph2Model::selectedVertexFieldMinValue() const {
  if (_selectedVertexField >= 0)
    return vertexFields[_selectedVertexField]->minValue;
  return 0;
}

double CellGraph2Model::selectedVertexFieldMaxValue() const {
  if (_selectedVertexField >= 0)
    return vertexFields[_selectedVertexField]->maxValue;
  return 0;
}

aergia::Color FieldModelInterface::colorAt(size_t id) {
  return aergia::COLOR_TRANSPARENT;
}
