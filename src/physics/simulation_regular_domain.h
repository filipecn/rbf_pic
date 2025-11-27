#ifndef FUROO_PHYSICS_SIMULATION_REGULAR_DOMAIN_H
#define FUROO_PHYSICS_SIMULATION_REGULAR_DOMAIN_H

#include <blas/linear_vector.h>
#include <common/debug.h>
#include <functional>
#include <geometry/point.h>
#include <geometry/transform.h>
#include <geometry/vector.h>
#include <physics/simulation_domain.h>
#include <structures/regular_grid.h>
#include <vector>

// Index map
// vn(m+1) ----- v ---------- v -------- ... ------- v(n+1)(m+1)
//  |           |            |            |          |
//  |  c(n-1)m  |            |           ...  cnm-1  |
//  |           |            |            |          |
// 2(m+1) ----- v ---------- v --------- ... ------- 3m+2
//  |           |            |            |          |
//  |    cm     |    cm+1    |    cm+2   ...  c2m-1  |
//  |           |            |            |          |
// m+1 --vm -- m+2 --vm+1---m+3 - vm+2-- ... v2m-1 - 2m+1
//  |           |            |            |          |
// u0    c0    u1      c1    u2     c2   ...  cm-1   um
//  |           |            |            |          |
//  0 ---v0---- 1 ----v1---  2 ----v2--- ...  vm-1 - m
namespace furoo {

class SimRegularDomain2 : public SimulationDomain2 {
 private:
  RegularGrid _grid;
  std::vector<LinearVector> _cellCenteredScalarFields;
  std::vector<std::vector<unsigned char>> _cellCenteredByteFields;
  std::vector<LinearVector> _uFaceFields;
  std::vector<LinearVector> _vFaceFields;
  Transform2D _cellToWorld, _cellToGrid;
  Transform2D _uFaceToWorld, _uFaceToGrid;
  Transform2D _vFaceToWorld, _vFaceToGrid;
  BicubicInterpolant _bicubicInterpolant;

 public:
  SimRegularDomain2();

  SimRegularDomain2(size_t);

  virtual ~SimRegularDomain2() {}

  unsigned uFacesCount() const;

  RegularGrid grid() const;

  unsigned cellIdFromIJ(Point2i) const;

  unsigned uFaceIdFromIJ(Point2i) const;

  unsigned vFaceIdFromIJ(Point2i) const;

  Point2i ijFromCellId(unsigned) const;

  Point2i ijFromFaceId(unsigned) const;

  int faceIdFromFaceId(unsigned, SimulationDomain2::Orientation) const;

  int faceIdFromCellId(unsigned, SimulationDomain2::Orientation) const;

  int cellIdFromFaceId(unsigned, SimulationDomain2::Orientation) const;

  int cellIdFromCellId(unsigned, SimulationDomain2::Orientation) const;

  int cellId(Point2d) const override;

  void setBoundaryHandler(DomainBoundaryHandler2 *bh) override;

  void setInterpolationMethod(InterpolationMethod) override;

  Vector2d smallestCellSize() const override;

  BBox2D region() const override;

  BBox2D cellRegion(size_t) const override;

  Point2d cellCenterPosition(size_t cellId) const override;

  Point2d faceCenterPosition(size_t faceId) const override;

  double &faceCenteredScalar(size_t, size_t) override;

  void iterateStencilAtCellCenter(
      size_t id, const std::function<void(int nid, double w, Point2d p,
                                          SimulationDomain2::Orientation o)> &f)
  const override;

  void iterateStencilAtCellCenter(
      size_t id,
      const std::function<void(int nid, double w,
                               DomainBoundaryHandler2::BoundaryType, double)>
      &f) const override;

  std::vector<SimulationDomain2::NeighborOrientation>
  faceNeighborhood(size_t faceId) const override;

  std::vector<SimulationDomain2::NeighborOrientation>
  cellNeighborhood(size_t cellId) const override;

  std::vector<SimulationDomain2::NeighborOrientation>
  cellFaces(unsigned cellId) const override;

  void
  iterateDxScalar(size_t fieldId,
                  const std::function<void(double &, Point2d)> &f) override;

  void iterateDxScalar(size_t fieldId,
                       const std::function<void(size_t, double &)> &f) override;

  void
  iterateDyScalar(size_t fieldId,
                  const std::function<void(double &, Point2d)> &f) override;

  void iterateDyScalar(size_t fieldId,
                       const std::function<void(size_t, double &)> &f) override;

  void
  iterateCellScalar(size_t fieldId,
                    const std::function<void(double &, Point2d)> &f) override;

  void
  iterateCellScalar(size_t fieldId,
                    const std::function<void(size_t, double &)> &f) override;

  void iterateCellByte(size_t fieldId,
                       const std::function<void(unsigned char &)> &f) override;

  void
  iterateVertexScalar(size_t fieldId,
                      const std::function<void(double &, Point2d)> &f) override;

  size_t addFaceCenteredScalarField(double) override;

  size_t addCellCenteredScalarField(double) override;

  size_t addCellCenteredByteField(unsigned char) override;

  size_t addVertexCenteredScalarField(double) override;

  size_t cellCount() const override;

  size_t faceCount() const override;

  size_t vertexCount() const override;

  double sampleCellCenteredScalar(size_t, Point2d) const override;

  double gradientXAtCellCenteredScalar(size_t fieldId,
                                       size_t cellId) const override;

  double gradientYAtCellCenteredScalar(size_t fieldId,
                                       size_t cellId) const override;

  double gradientXAtFaceCenteredScalar(size_t fieldId,
                                       size_t faceId) const override;

  double gradientYAtFaceCenteredScalar(size_t fieldId,
                                       size_t faceId) const override;

  double gradientXYAtCellCenteredScalar(size_t fieldId, size_t cellId) const;

  double gradientXYAtFaceCenteredScalar(size_t fieldId, size_t faceId) const;

  virtual double faceDivergentAtCellCenter(unsigned int fieldId,
                                           unsigned int cellId) const override;

  double cellGradientAtFaceCenter(unsigned int, unsigned int) const override;

  Vector2d sampleFaceCenteredScalar(size_t, Point2d) const override;

  double sampleVFaceField(size_t, Point2d) const;

  double sampleUFaceField(size_t, Point2d) const;

  Vector2d sampleFaceFieldAtCellCenter(size_t fieldId, size_t cellId) const;

  void sampleDxScalar(size_t, const std::vector<Point2d> &,
                      LinearVector &) const override;

  void sampleDyScalar(size_t, const std::vector<Point2d> &,
                      LinearVector &) const override;

  void extrapolateToDomain(size_t fieldId) override;

  /** Computes the cell of the grid a point **p** intersects.
   * \param target point
   * \param list of vertex ids (cellIds)
   * \param list of vertex positions
   * The lists of vertices follows the following order:
   *              2 ----- 3
   *              |       |
   *              0 ----- 1
   * \returns in cell coordinates of point **p**
   */
  Point2d cellCenteredGridCell(Point2d, unsigned *, Point2d * = nullptr) const;

  Point2d uFaceCenteredGridCell(Point2d,
                                unsigned *,
                                double * = nullptr,
                                unsigned = 0,
                                Point2d * = nullptr) const;

  Point2d vFaceCenteredGridCell(Point2d,
                                unsigned *,
                                double * = nullptr,
                                unsigned = 0,
                                Point2d * = nullptr) const;

  // random access
  double scalarAtCell(unsigned, unsigned) const override;

  double &scalarAtCell(unsigned, unsigned) override;

  double scalarAtFace(unsigned, unsigned) const override;

  double &scalarAtFace(unsigned, unsigned) override;

  unsigned char byteAtCell(unsigned, unsigned) const override;

  unsigned char &byteAtCell(unsigned, unsigned) override;

  void setScalarAtCell(unsigned fieldId, unsigned cellId, double value);

  void setScalarAtFace(unsigned fieldId, unsigned faceId, double value);
};
} // namespace furoo

#endif // FUROO_PHYSICS_SIMULATION_REGULAR_DOMAIN_H
