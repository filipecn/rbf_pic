#ifndef FUROO_FUNCTIONS_H
#define FUROO_FUNCTIONS_H

#include <geometry/point.h>

namespace furoo {
class Functions {
 public:
  static Vector2d enright(Point2d p) {
    return {
        2. * SQR(sin(PI * p.x())) * sin(PI * p.y())
            * cos(PI * p.y()),
        -2. * SQR(sin(PI * p.y())) * sin(PI * p.x())
            * cos(PI * p.x())};
  }

  static Vector2d enrightGradient(Point2d p) {
    return {PI * sin(2. * PI * p.x()) * sin(2. * PI * p.y()),
            -PI * sin(2. * PI * p.x()) * sin(2. * PI * p.y())};
  }
};

}
#endif //FUROO_FUNCTIONS_H
