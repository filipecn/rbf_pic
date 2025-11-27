#include "algorithms/fast_marching_method.h"
#include <queue>

namespace furoo {

void FastMarchingMethod2::propagate(SimRegularDomain2 *domain, unsigned field,
                                    unsigned frozen,
                                    const std::vector<unsigned> &sources) {
  auto s = domain->grid().gridSize();
  auto cellId = [&s](unsigned i, unsigned j) -> unsigned {
    return j * s.x() + i;
  };
  auto cellvId = [&s](Vector2i ij) -> unsigned {
    return ij[1] * s.x() + ij[0];
  };
  int m = s.x();
  int n = s.y();
  auto safeScalar = [&](Vector2i ij) -> double {
    if (ij[0] < 0 || ij[0] >= m || ij[1] < 0 || ij[1] >= n)
      return INFINITY;
    return domain->scalarAtCell(field, cellvId(ij));
  };
  struct Point_ {
    Point_(unsigned I, unsigned J, double D) : i(I), j(J), t(D) {}
    unsigned i, j;
    double t;
  };
  auto cmp = [](Point_ a, Point_ b) { return a.t > b.t; };
  std::priority_queue<Point_, std::vector<Point_>, decltype(cmp)> q(cmp);
  for (auto p : sources) {
    unsigned j = p / m;
    unsigned i = p % m;
    q.push(Point_(i, j, domain->scalarAtCell(field, p)));
  }
  while (!q.empty()) {
    Point_ p = q.top();
    q.pop();
    domain->byteAtCell(frozen, cellId(p.i, p.j)) = 2;
    domain->scalarAtCell(field, cellId(p.i, p.j)) = p.t;
    Vector2i dir[4] = {Vector2i(-1, 0), Vector2i(1, 0), Vector2i(0, -1),
                       Vector2i(0, 1)};
    for (int i = 0; i < 4; i++) {
      Vector2i ij = Vector2i(p.i, p.j) + dir[i];
      if (ij[0] < 0 || ij[0] >= m || ij[1] < 0 || ij[1] >= n ||
          domain->byteAtCell(frozen, cellId(ij[0], ij[1])))
        continue;
      float T1 = std::min(safeScalar(ij + dir[0]), safeScalar(ij + dir[1]));
      float T2 = std::min(safeScalar(ij + dir[2]), safeScalar(ij + dir[3]));
      float Tmin = std::min(T1, T2);
      float d = fabs(T1 - T2);
      if (d < 1.f)
        domain->scalarAtCell(field, cellvId(ij)) =
            (T1 + T2 + sqrtf(2.f * SQR(1.f) - SQR(d))) / 2.f;
      else
        domain->scalarAtCell(field, cellvId(ij)) = Tmin + 1.f;
      q.push(Point_(ij[0], ij[1], domain->scalarAtCell(field, cellvId(ij))));
      domain->byteAtCell(frozen, cellvId(ij)) = 1;
    }
  }
}

} // furoo namespace
