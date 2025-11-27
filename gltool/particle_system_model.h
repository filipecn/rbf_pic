#ifndef FUROO_PARTICLE_SYSTEM_MODEL_H
#define FUROO_PARTICLE_SYSTEM_MODEL_H

#include "graphic_element_set.h"
#include "ponos_wrapper.h"
#include <furoo.h>
#include <sstream>

template <int D> class ParticleSystemModel : public aergia::SceneObject {
public:
  ParticleSystemModel();
  explicit ParticleSystemModel(furoo::ParticleSystem<double, D> *ps,
                               double radius);
  void set(furoo::ParticleSystem<double, D> *ps, double radius);
  // Interface
  void draw(const aergia::CameraInterface *camera, ponos::Transform t) override;
  bool intersect(const ponos::Ray3 &r, float *t) override;
  /// Defines a color palette for a property
  /// \tparam T property type
  /// \param i property id
  /// \param p color palette
  template <typename T>
  void addPropertyColor(size_t i, aergia::ColorPalette p) {
    if (std::is_same<T, double>::value)
      _paletteD[i] = std::move(p);
  }
  template <typename T> void selectProperty(size_t i) {
    if (std::is_same<T, double>::value && _paletteD.find(i) != _paletteD.end())
      _curProperty = static_cast<int>(i);
  }
  void addVectorGroup(furoo::Vector<size_t, D> g) {
    vectorGroups.emplace_back(g);
  }
  void selectPropertyGroup(size_t i) {
    _curPropertyGroup = i;
    updateVelocities();
  }
  /// update buffers
  void update();
  void updateParticles();
  void updateVelocities();
  bool showVelocities = false;
  double particleRadius = 0.07; ///< particle's sphere radius
  float velocityScale = 1.f;
  aergia::Color particleColor = aergia::COLOR_BLUE;
  furoo::ParticleSystem<double, D> *ps =
      nullptr; ///< particle system's reference
private:
  std::vector<furoo::Vector<size_t, D>> vectorGroups;
  int _curPropertyGroup = -1;
  int _curProperty = -1;    ///< current property to be drawn
  int _activeParticle = -1; ///< selected particle
  std::map<size_t, aergia::ColorPalette>
      _paletteD; ///< color pallets for double type properties
  // instances fields
  std::shared_ptr<GraphicElementSet> _particlesSet;
  std::shared_ptr<GraphicElementSet> _velocitiesSet;
};

#include "particle_system_model.inl"

typedef ParticleSystemModel<2> ParticleSystemModel2;
typedef ParticleSystemModel<3> ParticleSystemModel3;

#endif // FUROO_PARTICLE_SYSTEM_MODEL_H
