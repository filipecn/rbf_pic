#ifndef FUROO_ALGORITHMS_FAST_MARCHING_METHOD_H
#define FUROO_ALGORITHMS_FAST_MARCHING_METHOD_H

#include <physics/simulation_regular_domain.h>

namespace furoo {

class FastMarchingMethod2 {
public:
  static void propagate(SimRegularDomain2 *domain, unsigned field,
                        unsigned frozen, const std::vector<unsigned> &sources);
};

} // furoo namespace

#endif // FUROO_ALGORITHMS_FAST_MARCHING_METHOD_H
