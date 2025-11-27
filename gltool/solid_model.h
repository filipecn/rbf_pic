#ifndef FUROO_GLTOOL_SOLID_MODEL_H
#define FUROO_GLTOOL_SOLID_MODEL_H

#include <furoo.h>
#include <ponos_wrapper.h>

template <int D> class SolidModel : public aergia::SceneObject {
public:
  SolidModel(furoo::SolidInterface<D> *solid) : _solid(solid) {
    _mesh = _solid->mesh();
    _sceneMeshObject.reset(new furoo::SolidObject<D>(_mesh));
  }

  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    if (_sceneMeshObject)
      _sceneMeshObject->draw(camera, t);
  }

  furoo::SolidInterface<D> *solid() { return _solid; }

private:
  std::shared_ptr<furoo::SolidObject<D>> _sceneMeshObject;
  furoo::Mesh<D> *_mesh = nullptr;
  furoo::SolidInterface<D> *_solid = nullptr;
};

#endif // FUROO_GLTOOL_SOLID_MODEL_H