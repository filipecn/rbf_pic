#ifndef FUROO_SEARCHRENDERER_H
#define FUROO_SEARCHRENDERER_H

#include <physics/particle_system2.h>
#include <aergia/scene/scene_object.h>
#include <physics/simulation_quad_tree.h>

namespace furoo {

class SearchRenderer : public aergia::SceneObject {
 public:
  explicit SearchRenderer(ParticleSystem2 *, SimQuadTree *);
  ~SearchRenderer() {}
  void draw() override;

  void set(ParticleSystem2 *, SimQuadTree *);
  bool intersect(const ponos::Ray3 &r, float *t) override;

 private:
  Point2d mouse;
  ParticleSystem2 *_ps;
  int _selectedParticle;
  std::unique_ptr<PointZGrid2> _faces;
  SimQuadTree *_domain;
  int _selectedCell;
  int _selectedFace;
};

}  // namespace furoo

#endif //FUROO_SEARCHRENDERER_H
