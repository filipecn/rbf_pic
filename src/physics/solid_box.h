#ifndef FUROO_SOLID_BOX_H
#define FUROO_SOLID_BOX_H

#include <geometry/bbox.h>
#include <geometry/transform.h>
#include <physics/solid_interface.h>

namespace furoo {

template <int D> class SolidBox : public SolidInterface<D> {
public:
  explicit SolidBox(const BBox<double, D> &box);
  bool contains(const Point<double, D> &point) const override;
  Definitions::Intersection<D>
  intersectionPoint(const Point<double, D> &origin,
                    const Point<double, D> &destination) const override;
  Mesh<D> *mesh() override;
  BBox<double, D> box() const;

private:
  Mesh<D> _mesh;
  Transform<D> _toUnitBox;
  BBox<double, D> _box;
};

#include "solid_box.inl"

} // namespace furoo

#endif // FUROO_SOLID_BOX_H
