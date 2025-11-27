#ifndef FUROO_MESH_H
#define FUROO_MESH_H

#include <geometry/bbox.h>
#include <geometry/point.h>
#include <geometry/queries.h>

namespace furoo {

/// Simple mesh structure class for storing a set of vertices and elements.
/// Elements are faces in the mesh that represent line segments in 2D and
/// triangles in 3D.
/// \tparam D dimensions
template <int D> class Mesh {
public:
  /// \return number of elements
  size_t elementCount() const;

  /// \return number of vertices
  size_t vertexCount() const;

  /// \param i vertex index
  /// \return vertex position
  Point<double, D> vertex(size_t i) const;

  /// \param i element index
  /// \return element
  Vector<int, D> element(size_t i) const;

  /// The axis aligned bounding box containing the element
  /// \param i element index
  /// \return bounding box
  BBox<double, D> elementBBox(size_t i) const;

  /// Intersects a element and a ray
  /// \param elementId element index
  /// \param origin ray origin point
  /// \param direction ray direction
  /// \param hit [optional] parametric ray coordinate of the intersection point
  /// \return true if intersection exists
  bool intersectElement(size_t elementId, Point<double, D> origin,
                        Vector<double, D> direction,
                        double *hit = nullptr) const;

  /// Add a position to the mesh
  /// \param p vertex position
  void addVertex(Point<double, D> p);

  /// Add an element (composed by vertex indices)
  /// \param e vertices ids
  void addElement(Vector<int, D> e);

private:
  std::vector<Point<double, D>> _vertices; //!< list of vertices
  std::vector<Vector<int, D>> _elements;   //!< list of elements
};

#include "mesh.inl"

} // namespace furoo

#endif // FUROO_MESH_H
