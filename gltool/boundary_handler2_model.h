#ifndef FUROO_BOUNDARY_HANDLER2_MODEL_H
#define FUROO_BOUNDARY_HANDLER2_MODEL_H

#include <aergia/aergia.h>
#include <furoo.h>

namespace furoo {

class DomainBoundaryHandler2Model : public aergia::SceneObject {

 public:
  DomainBoundaryHandler2Model(SimRegularDomain2 *, DomainBoundaryHandler2 *);
  void draw() override;
  bool drawCells;
 private:
  SimRegularDomain2 *_dm;
  DomainBoundaryHandler2 *_bh;
};

} // furoo namespace

#endif //FUROO_BOUNDARY_HANDLER2_MODEL_H
