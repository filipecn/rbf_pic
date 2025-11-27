#ifndef FUROO_MESH_MODEL_H
#define FUROO_MESH_MODEL_H

#include <aergia/aergia.h>
#include <furoo.h>
#include <graphic_element_set.h>

template <int D> class MeshModel : public aergia::SceneObject {
public:
  MeshModel(furoo::Mesh<D> &mesh) : _mesh(mesh) {
    _rawMesh.primitiveType = (D == 2)
                                 ? ponos::GeometricPrimitiveType::LINES
                                 : ponos::GeometricPrimitiveType::TRIANGLES;
    _rawMesh.meshDescriptor.count = mesh.elementCount();
    _rawMesh.meshDescriptor.elementSize = D;
    _rawMesh.positionDescriptor.count = mesh.vertexCount();
    _rawMesh.positionDescriptor.elementSize = D;
    for (size_t i = 0; i < mesh.vertexCount(); i++)
      for (int d = 0; d < D; d++)
        _rawMesh.addPosition({mesh.vertex(i)[d]});
    for (size_t i = 0; i < mesh.elementCount(); i++)
      for (int d = 0; d < D; d++)
        _rawMesh.addFace({mesh.element(i)[d]});
    for (auto &i : _rawMesh.indices)
      i.normalIndex = i.texcoordIndex = i.positionIndex;
    _rawMesh.apply(ponos::scale(2, 2, 2));
    _rawMesh.apply(ponos::translate(ponos::vec3(0.5,0,1)));
    _rawMesh.computeBBox();
    _rawMesh.splitIndexData();
    _rawMesh.buildInterleavedData();
    static auto ic =
        ponos::create_icosphere_mesh(ponos::Point3(), 1.f, 3, false, false);
    //    _sceneMesh.reset(new aergia::SceneMesh(*ic));
    _meshSet.reset(new GraphicElementSet(&_rawMesh));
    //    _meshSet.reset(new GraphicElementSet(ic));

    _meshSet->resize(1);
    // unactivate dead particles, unfortunatelly we have to kill them all first
    // and make active ones live again
    {
      auto a = _meshSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _meshSet->count());
    }
    // update positions
    for (size_t id = 0; id < 1; id++) {
      auto m = _meshSet->transformBuffer(id);
      furoo::Vector<double, D> size(1);
      float t[16];
      (ponos::translate(ponos::vec3(0, 0, 0)) *
       ponos::scale(size[0], size[1], (D == 3) ? size[2] : 0.f))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];
      aergia::Color color = aergia::COLOR_RED;
      auto c = _meshSet->colorBuffer(id);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = 0.2;
      auto a = _meshSet->activeBuffer(id);
      a[0] = 1;
    }
  }
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    _meshSet->draw(camera, t);
    //    aergia::glColor(aergia::COLOR_BLACK);
    //    _sceneMesh->bind();
    //    glDrawElements(_sceneMesh->indexBuffer()->bufferDescriptor.elementType,
    //                   _sceneMesh->indexBuffer()->bufferDescriptor.elementSize
    //                   *
    //                       _sceneMesh->indexBuffer()->bufferDescriptor.elementCount,
    //                   GL_UNSIGNED_INT, 0);
    //    _sceneMesh->unbind();
  }

private:
  std::shared_ptr<GraphicElementSet> _meshSet;
  //  std::shared_ptr<aergia::SceneMesh> _sceneMesh;
  ponos::RawMesh _rawMesh;
  furoo::Mesh<D> &_mesh;
};

#endif // FUROO_MESH_MODEL_H
