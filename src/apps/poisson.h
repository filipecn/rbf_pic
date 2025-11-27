#ifndef FURRO_APPS_POISSON_H
#define FURRO_APPS_POISSON_H

#include <structures/cell_graph.h>

namespace furoo {

class Poisson {
protected:
  CellGraph2 _structure;
  size_t _boundaryConditionsFieldId;
  size_t _cellValueFieldId;
  size_t _faceValueFieldId;
  size_t _gradientFieldId;

public:
  Poisson();
};
} // namespace furoo

#endif
