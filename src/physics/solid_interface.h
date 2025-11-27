#ifndef FUROO_SOLID_INTERFACE_H
#define FUROO_SOLID_INTERFACE_H

#include <common/definitions.h>
#include <geometry/mesh.h>

namespace furoo {

/// Interface for physical solid objects
/// \tparam D dimensions
template <int D> class SolidInterface {
public:
  /// checks if point is inside solid
  /// \param point query point
  /// \return true if point is inside solid
  virtual bool contains(const Point<double, D> &point) const = 0;

  /// Given a initial and a final position of a particle, computes the
  /// intersection point with the solid. The final position MUST be inside the
  /// solid, otherwise, an exception is thrown
  /// \param origin origin point
  /// \param destination destination point
  /// \return the intersection information
  virtual Definitions::Intersection<D>
  intersectionPoint(const Point<double, D> &origin,
                    const Point<double, D> &destination) const = 0;

  /// It might be useful for solids to provide their mesh representation
  /// \return a pointer to its mesh object
  virtual Mesh<D> *mesh() = 0;
};

} // namespace furoo

#endif // FUROO_SOLID_INTERFACE_H
