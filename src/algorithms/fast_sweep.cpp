#include "algorithms/fast_sweep.h"
#include <functional>
#include <geometry/numeric.h>

namespace furoo {

void FastSweep2::propagate(SimRegularDomain2 *domain, unsigned field,
                           unsigned mask) {
  auto solveDistance = [](double p, double q, double &r) {
    double d = std::min(p, q) + 1.;
    if (d > std::max(p, q))
      d = (p + q + std::sqrt(2 - SQR(p - q))) / 2.;
    if (d < r)
      r = d;
  };
  auto s = domain->grid().gridSize();
  auto cellId = [&s](unsigned i, unsigned j) -> unsigned {
    return j * s.x() + i;
  };
  int m = s.x();
  int n = s.y();
  // +1 +1
  for (int j = 1; j < n; j++)
    for (int i = 1; i < m; i++) {
      if (!domain->byteAtCell(mask, cellId(i, j)))
        solveDistance(domain->scalarAtCell(field, cellId(i - 1, j)),
                      domain->scalarAtCell(field, cellId(i, j - 1)),
                      domain->scalarAtCell(field, cellId(i, j)));
    }
  // -1 -1
  for (int j = n - 2; j >= 0; j--)
    for (int i = m - 2; i >= 0; i--) {
      if (!domain->byteAtCell(mask, cellId(i, j)))
        solveDistance(domain->scalarAtCell(field, cellId(i + 1, j)),
                      domain->scalarAtCell(field, cellId(i, j + 1)),
                      domain->scalarAtCell(field, cellId(i, j)));
    }
  // +1 -1
  for (int j = n - 2; j >= 0; j--)
    for (int i = 1; i < m; i++) {
      if (!domain->byteAtCell(mask, cellId(i, j)))
        solveDistance(domain->scalarAtCell(field, cellId(i - 1, j)),
                      domain->scalarAtCell(field, cellId(i, j + 1)),
                      domain->scalarAtCell(field, cellId(i, j)));
    }
  // -1 +1
  for (int j = 1; j < n; j++)
    for (int i = m - 2; i >= 0; i--) {
      if (!domain->byteAtCell(mask, cellId(i, j)))
        solveDistance(domain->scalarAtCell(field, cellId(i + 1, j)),
                      domain->scalarAtCell(field, cellId(i, j - 1)),
                      domain->scalarAtCell(field, cellId(i, j)));
    }
}

} // furoo namespace
