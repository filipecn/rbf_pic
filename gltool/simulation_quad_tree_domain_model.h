#ifndef FUROO_SIMULATION_QUAD_TREE_DOMAIN_H_H
#define FUROO_SIMULATION_QUAD_TREE_DOMAIN_H_H

#include <furoo.h>
#include <aergia/scene/scene_object.h>

namespace furoo {

class SimulationQuadTreeDomainModel : public aergia::SceneObject {
 public:
  explicit SimulationQuadTreeDomainModel(SimQuadTree *,
                                         DomainBoundaryHandler2 *);
  void set(SimQuadTree *, DomainBoundaryHandler2 *);
  void draw() override;
  bool intersect(const ponos::Ray3 &r, float *t) override;
  void updatePressure();
  void updateDivergence();
  void drawPressure();
  void drawDivergence();
  void drawMaterial();
  aergia::Color cellCenterColor;
  aergia::Color faceCenterColor;
 private:
  bool _drawMaterial;
  bool _drawVelocities;
  bool _drawPressure;
  bool _drawDivergence;
  std::vector<BBox2D> boxes;
  int _selectedCell;
  int _selectedFace;
  aergia::ColorPalette pressurePalette;
  double minPressure, maxPressure;
  double minDivergence, maxDivergence;
  std::vector<double> divergence;
  std::unique_ptr<PointZGrid2> _faces;
  SimQuadTree *_domain;
  DomainBoundaryHandler2 *_bh;
  double m, M;
};

} // furoo namespace

#endif //FUROO_SIMULATION_QUAD_TREE_DOMAIN_H_H
