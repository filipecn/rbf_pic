#ifndef FUROO_GRAPHIC_ELEMENT_H
#define FUROO_GRAPHIC_ELEMENT_H

#include <aergia/aergia.h>
#include <furoo.h>

template <int D> inline ponos::RawMesh *createArrow() {
  ponos::RawMesh *mesh = new ponos::RawMesh();
  furoo::Point<double, D> a, b(1.);
  for (int d = 0; d < D; d++)
    mesh->addPosition({a[d]});
  for (int d = 0; d < D; d++)
    mesh->addPosition({b[d]});
  mesh->addFace({0, 1});
  // fix normal and uvs indices
  for (auto &i : mesh->indices)
    i.normalIndex = i.texcoordIndex = i.positionIndex;
  size_t vertexCount = mesh->positions.size() / 2;
  // describe mesh
  mesh->primitiveType = ponos::GeometricPrimitiveType::LINES;
  mesh->meshDescriptor.count = mesh->indices.size() / 2;
  mesh->meshDescriptor.elementSize = 2;
  mesh->positionDescriptor.count = vertexCount;
  mesh->positionDescriptor.elementSize = D;
  mesh->computeBBox();
  mesh->splitIndexData();
  mesh->buildInterleavedData();
  return mesh;
}

class GraphicElementSet : public aergia::SceneObject {
public:
  explicit GraphicElementSet(ponos::RawMesh *rawMesh);
  explicit GraphicElementSet(ponos::RawMeshSPtr &rawMesh);
  // INTERFACE
  void draw(const aergia::CameraInterface *camera, ponos::Transform t) override;
  bool intersect(const ponos::Ray3 &r, float *t) override;
  size_t count() const;
  void resize(size_t n);
  uint *activeBuffer(size_t i = 0);
  float *colorBuffer(size_t i = 0);
  float *transformBuffer(size_t i = 0);

private:
  void setShader();
  // instances fields
  size_t _posBuffer = 0, _colorBuffer = 0, _activeBuffer = 0;
  ponos::RawMeshSPtr _rawMesh;
  std::shared_ptr<aergia::SceneMesh> _sceneMesh;
  std::shared_ptr<aergia::InstanceSet> _instances;
  size_t _instancesCount = 1;
};

#endif // FUROO_GRAPHIC_ELEMENT_H
