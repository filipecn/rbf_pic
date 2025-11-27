// DEPRECATED
#ifndef FUROO_PHYSICS_SIMULATION_DOMAIN_H
#define FUROO_PHYSICS_SIMULATION_DOMAIN_H

#include <functional>
#include <geometry/point.h>
#include <memory>
#include <physics/domain_boundary_handler.h>

namespace furoo {

class SimulationDomain2 {
 public:
  enum class InterpolationMethod {
    BILINEAR, BICUBIC, CATMULL, RBF, OTHER
  };
  enum class Orientation {
    LEFT, RIGHT, BOTTOM, TOP, OTHER
  };

  struct NeighborOrientation {
    NeighborOrientation(int i, Orientation o) : id(i), orientation(o) {}

    int id;
    Orientation orientation;

    friend std::ostream &operator<<(std::ostream &o,
                                    const NeighborOrientation &ng) {
      switch (ng.orientation) {
        case Orientation::LEFT:o << "(" << ng.id << ", LEFT)" << std::endl;
          break;
        case Orientation::BOTTOM:o << "(" << ng.id << ", BOTTOM)" << std::endl;
          break;
        case Orientation::RIGHT:o << "(" << ng.id << ", RIGHT)" << std::endl;
          break;
        case Orientation::TOP:o << "(" << ng.id << ", TOP)" << std::endl;
          break;
        case Orientation::OTHER:o << "(" << ng.id << ", OTHER)" << std::endl;
          break;
      }
      return o;
    }
  };

  friend std::ostream &operator<<(std::ostream &o, const Orientation &ori) {
    switch (ori) {
      case Orientation::LEFT:o << "(LEFT)" << std::endl;
        break;
      case Orientation::BOTTOM:o << "(BOTTOM)" << std::endl;
        break;
      case Orientation::RIGHT:o << "(RIGHT)" << std::endl;
        break;
      case Orientation::TOP:o << "(TOP)" << std::endl;
        break;
      case Orientation::OTHER:o << "(OTHER)" << std::endl;
        break;
    }
    return o;
  }

  SimulationDomain2() = default;

  virtual ~SimulationDomain2() = default;

  virtual void setBoundaryHandler(DomainBoundaryHandler2 *bh) = 0;

  virtual void setInterpolationMethod(InterpolationMethod) = 0;

  virtual double &faceCenteredScalar(size_t, size_t) = 0;

  virtual Point2d cellCenterPosition(size_t cellId) const = 0;

  virtual Point2d faceCenterPosition(size_t faceId) const = 0;
  /// \param p target point
  /// \return id of containing cell, -1 if p is outside domain
  virtual int cellId(Point2d p) const = 0;

  virtual Vector2d smallestCellSize() const = 0;

  virtual BBox2D region() const = 0;
  /// \return cell domain region
  virtual BBox2D cellRegion(size_t) const = 0;

  virtual void iterateStencilAtCellCenter(
      size_t id,
      const std::function<void(int nid, double w, Point2d p, Orientation o)> &f)
  const = 0;

  virtual void iterateStencilAtCellCenter(
      size_t id,
      const std::function<void(int nid, double w,
                               DomainBoundaryHandler2::BoundaryType, double)>
      &f) const = 0;

  virtual std::vector<NeighborOrientation>
  faceNeighborhood(size_t faceId) const = 0;

  virtual std::vector<NeighborOrientation>
  cellNeighborhood(size_t cellId) const = 0;

  virtual std::vector<NeighborOrientation>
  cellFaces(unsigned cellId) const = 0;

  virtual void
  iterateDxScalar(size_t fieldId,
                  const std::function<void(double &, Point2d)> &f) = 0;

  virtual void
  iterateDxScalar(size_t fieldId,
                  const std::function<void(size_t, double &)> &f) = 0;

  virtual void
  iterateDyScalar(size_t fieldId,
                  const std::function<void(double &, Point2d)> &f) = 0;

  virtual void
  iterateDyScalar(size_t fieldId,
                  const std::function<void(size_t, double &)> &f) = 0;
  /// \param fieldId
  /// \param f callback for each cell
  virtual void
  iterateCellScalar(size_t fieldId,
                    const std::function<void(double &, Point2d)> &f) = 0;
  /// \param fieldId
  /// \param f callback for each cell
  virtual void
  iterateCellScalar(size_t fieldId,
                    const std::function<void(size_t, double &)> &f) = 0;

  virtual void
  iterateCellByte(size_t fieldId,
                  const std::function<void(unsigned char &)> &f) = 0;

  virtual void
  iterateVertexScalar(size_t fieldId,
                      const std::function<void(double &, Point2d)> &f) = 0;

  virtual size_t addFaceCenteredScalarField(double v) = 0;

  virtual size_t addCellCenteredScalarField(double v) = 0;

  virtual size_t addCellCenteredByteField(unsigned char v) = 0;

  virtual size_t addVertexCenteredScalarField(double v) = 0;
  /// \return active cell count
  virtual size_t cellCount() const = 0;
  /// \return active face count
  virtual size_t faceCount() const = 0;

  virtual size_t vertexCount() const = 0;

  virtual double sampleCellCenteredScalar(size_t, Point2d) const = 0;

  virtual double gradientXAtCellCenteredScalar(size_t, size_t) const = 0;

  virtual double gradientYAtCellCenteredScalar(size_t, size_t) const = 0;

  virtual double gradientXAtFaceCenteredScalar(size_t, size_t) const = 0;

  virtual double gradientYAtFaceCenteredScalar(size_t, size_t) const = 0;

  virtual double faceDivergentAtCellCenter(unsigned int,
                                           unsigned int) const = 0;

  virtual double cellGradientAtFaceCenter(unsigned int, unsigned int) const = 0;

  virtual Vector2d sampleFaceCenteredScalar(size_t, Point2d) const = 0;

  virtual void sampleDxScalar(size_t, const std::vector<Point2d> &,
                              LinearVector &) const = 0;

  virtual void sampleDyScalar(size_t, const std::vector<Point2d> &,
                              LinearVector &) const = 0;

  virtual void extrapolateToDomain(size_t fieldId) = 0;

  // random access
  virtual double scalarAtCell(unsigned, unsigned) const = 0;

  virtual double &scalarAtCell(unsigned, unsigned) = 0;

  virtual double scalarAtFace(unsigned, unsigned) const = 0;

  virtual double &scalarAtFace(unsigned, unsigned) = 0;

  virtual unsigned char byteAtCell(unsigned, unsigned) const = 0;

  virtual unsigned char &byteAtCell(unsigned, unsigned) = 0;

 protected:
  DomainBoundaryHandler2 *_boundaryHandler{};
  InterpolationMethod _interpolant;
}; // namespace furoo
} // namespace furoo

#endif // FUROO_PHYSICS_SIMULATION_DOMAIN_H
