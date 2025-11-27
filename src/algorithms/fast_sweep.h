#ifndef FUROO_ALGORITHMS_FAST_SWEEP_H
#define FUROO_ALGORITHMS_FAST_SWEEP_H

#include <physics/simulation_regular_domain.h>

namespace furoo {

class FastSweep2 {
public:
  /**
   * \param
   * \param unsigned
   */
  static void propagate(SimRegularDomain2 *domain, unsigned field,
                        unsigned mask);
};

} // furoo namespace

#endif // FUROO_ALGORITHMS_FAST_SWEEP_H
