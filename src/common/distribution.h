#ifndef FUROO_DISTRIBUTION_H
#define FUROO_DISTRIBUTION_H

#include <common/definitions.h>
#include <geometry/bbox.h>
#include <geometry/point.h>

namespace furoo {

template <int D> class Distribution {
public:
  /// Generates equally spaced points over a line segment.
  /// Including start and end points.
  /// \param a segment's start point
  /// \param b segment's end point
  /// \param n number of additional points
  /// \return list of generated positions
  static std::vector<Point<double, D>> linear(Point<double, D> a,
                                              Point<double, D> b, size_t n);

  /// Generates equally spaced points over a face of a box
  /// Including its vertices
  /// \param bbox box
  /// \param side face side
  /// \param n number of additional points
  /// \return list of generated positions
  static std::vector<Point<double, D>> linear(BBox<double, D> bbox,
                                              Definitions::Side side, size_t n);
};

#include "distribution.inl"

} // namespace furoo

#endif // FUROO_DISTRIBUTION_H
