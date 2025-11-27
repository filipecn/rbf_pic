#ifndef FUROO_PARTICLE_SYSTEM2_MODEL_H
#define FUROO_PARTICLE_SYSTEM2_MODEL_H

#include <physics/particle_system2.h>
#include <aergia/scene/scene_object.h>
//#include <aergia/ui/text.h>

namespace furoo {

class ParticleSystem2Model : public aergia::SceneObject {
 public:
  explicit ParticleSystem2Model(ParticleSystem2 *, double);

  void draw() override;
  void update();

  bool intersect(const ponos::Ray3 &r, float *t) override;

  double dt;
  aergia::ColorPalette pallete;
  aergia::Color particleColor;
  double particleRadius;

 private:
  double _maxValue, _minValue;
  ParticleSystem2 *_ps;
  int _selectedParticle;
};

}  // namespace furoo

#endif  // FUROO_PARTICLE_SYSTEM2_MODEL_H
