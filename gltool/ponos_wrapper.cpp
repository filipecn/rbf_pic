#include "ponos_wrapper.h"

namespace furoo {

ponos::BBox2D Ponos::bbox2D(const BBox2D &b) {
  return ponos::BBox2D(ponos::Point2(b.lower().x(), b.lower().y()),
                       ponos::Point2(b.upper().x(), b.upper().y()));
}

ponos::Point2 Ponos::point2(const Point2d &p) {
  return ponos::Point2(p.x(), p.y());
}
ponos::BBox2D Ponos::bbox2D(const BBox2d &b) {
  return ponos::BBox2D(ponos::Point2(b.lower().x(), b.lower().y()),
                       ponos::Point2(b.upper().x(), b.upper().y()));
}

} // furoo namespace
