#ifndef FUROO_SOLVER_H
#define FUROO_SOLVER_H

namespace furoo {

template<int D>
class Solver {
 public:
  Solver();
  virtual ~Solver();

};

#include "physics/solver.inl"

}

#endif //FUROO_SOLVER_H
