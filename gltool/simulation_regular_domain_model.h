#ifndef FUROO_SIMULATION_REGULAR_DOMAIN_MODEL_H
#define FUROO_SIMULATION_REGULAR_DOMAIN_MODEL_H

#include <physics/simulation_regular_domain.h>
#include <aergia/scene/scene_object.h>
#include <aergia/ui/text.h>
#include <furoo.h>

namespace furoo {

class SimulationRegularDomain2Model : public aergia::SceneObject {
 public:
  explicit SimulationRegularDomain2Model(SimRegularDomain2 *);

  void draw() const override;

  bool intersect(const ponos::Ray3 &r, float *t) override;

  void updatePressure();
  void updateDivergence();

  aergia::Color cellCenterColor;
  aergia::Color faceCenterColor;
  bool drawVelocities;
  bool drawPressure;
  bool drawDivergence;

 private:
  int _selectedCell;
  int _selectedFace;
  SimRegularDomain2 *_domain;
  aergia::ColorPalette pressurePalette;
  double minPressure, maxPressure;
  double minDivergence, maxDivergence;
  std::vector<double> divergence;
  std::shared_ptr<aergia::Text> _text;
  std::unique_ptr<PointZGrid2> _faces;
};

}  // namespace furoo

#endif  // FUROO_SIMULATION_REGULAR_DOMAIN_MODEL_H
