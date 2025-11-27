// DEPRECATED
#ifndef FUROO_PHYSICS_SIMULATION_QUAD_TREE_H
#define FUROO_PHYSICS_SIMULATION_QUAD_TREE_H

#include <blas/linear_vector.h>
#include <physics/simulation_domain.h>
#include <structures/point_z_grid.h>
#include <structures/quad_tree.h>

namespace furoo {

class SimQuadTree : public SimulationDomain2 {
 public:
  SimQuadTree();

  template<typename... Args>
  SimQuadTree(Args &&... args) : _qt(QuadTree(std::forward<Args>(args)...)) {
    if (_qt.leafCount()) {
      ASSERT_FATAL(_qt.edgeCount());
    }
    updateEdgeIndices();
    this->_interpolant = SimulationDomain2::InterpolationMethod::RBF;
  }

  virtual ~SimQuadTree();

  void setPointSetStructures(PointSet2Interface *cs, PointSet2Interface *fsx,
                             PointSet2Interface *fsy);

  void iterateCells(Point2d, double, const std::function<void(size_t)> &) const;

  void iterateUFaces(Point2d, double,
                     const std::function<void(size_t)> &) const;
  void iterateVFaces(Point2d, double,
                     const std::function<void(size_t)> &) const;
  void iterateVFaces(Point2d,
                     const std::function<void(size_t)> &,
                     size_t = 1) const;
  void iterateUFaces(Point2d,
                     const std::function<void(size_t)> &,
                     size_t = 1) const;

  // INTERFACE

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

  size_t addFaceCenteredScalarField(double v) override;

  size_t addCellCenteredScalarField(double v) override;

  size_t addCellCenteredByteField(unsigned char v) override;

  size_t addVertexCenteredScalarField(double v) override;

  size_t cellCount() const override { return _qt.leafCount(); }

  size_t faceCount() const override { return _qt.edgeCount(); }

  size_t vertexCount() const override {
    NOT_IMPLEMENTED();
    return 0;
  }

  double sampleCellCenteredScalar(size_t, Point2d) const override;

  double gradientXAtCellCenteredScalar(size_t fieldId,
                                       size_t cellId) const override;

  double gradientYAtCellCenteredScalar(size_t fieldId,
                                       size_t cellId) const override;

  double gradientXAtFaceCenteredScalar(size_t fieldId,
                                       size_t faceId) const override;

  double gradientYAtFaceCenteredScalar(size_t fieldId,
                                       size_t faceId) const override;

  double faceDivergentAtCellCenter(unsigned int fieldId,
                                   unsigned int faceId) const override;

  double cellGradientAtFaceCenter(unsigned int, unsigned int) const override;

  Vector2d sampleFaceCenteredScalar(size_t fieldId,
                                    Point2d target) const override;

  double sampleVFaceField(size_t, Point2d) const;

  double sampleUFaceField(size_t, Point2d) const;

  void sampleDxScalar(size_t, const std::vector<Point2d> &,
                      LinearVector &) const override;

  void sampleDyScalar(size_t, const std::vector<Point2d> &,
                      LinearVector &) const override;

  void extrapolateToDomain(size_t fieldId) override;

  // random access
  double scalarAtCell(unsigned, unsigned) const override;

  double &scalarAtCell(unsigned, unsigned) override;

  double scalarAtFace(unsigned, unsigned) const override;

  double &scalarAtFace(unsigned, unsigned) override;

  unsigned char byteAtCell(unsigned, unsigned) const override;

  unsigned char &byteAtCell(unsigned, unsigned) override;

  const QuadTree &tree() const;

  void setTree(const QuadTree &);

  // private:
  void updateEdgeIndices();

  Point2d neighborCellCenter(size_t i, size_t e) const;

  QuadTree _qt;
  // scalar quantities
  std::vector<LinearVector> _faceCenteredScalarFields;
  std::vector<LinearVector> _cellCenteredScalarFields;
  std::vector<std::vector<unsigned char>> _cellCenteredByteFields;
  std::vector<LinearVector> _vertexCenterd;
  std::vector<size_t> _vEdgeIndices;
  std::vector<size_t> _uEdgeIndices;
  std::unique_ptr<PointSet2Interface> _cellCenters;
  std::unique_ptr<PointSet2Interface> _vFaceCenters;
  std::unique_ptr<PointSet2Interface> _uFaceCenters;
};

} // namespace furoo

#endif // FUROO_PHYSICS_SIMULATION_QUAD_TREE_H
