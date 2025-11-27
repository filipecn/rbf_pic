#ifndef FUROO_PHYSICS_DOMAIN_BOUNDARY_HANDLER_H
#define FUROO_PHYSICS_DOMAIN_BOUNDARY_HANDLER_H

#include <blas/linear_vector.h>
#include <physics/particle_system2.h>
#include <structures/quad_tree.h>

namespace furoo {

class SimulationDomain2;
class DomainBoundaryHandler2 {
public:
  enum class BoundaryType { CUSTOM, DIRICHLET, NEUMANN, ROBIN, NONE };
  enum class MaterialType { FLUID, SOLID, AIR, CUSTOM };
  struct Intersection {
    Intersection(Point2d p, Vector2d n, bool iv = false)
        : point(p), normal(n), isValid(iv) {}
    Point2d point;
    Vector2d normal;
    bool isValid;
  };
  DomainBoundaryHandler2();
  DomainBoundaryHandler2(BoundaryType bt);
  DomainBoundaryHandler2(BoundaryType left, BoundaryType bottom,
                         BoundaryType right, BoundaryType top);
  virtual ~DomainBoundaryHandler2();
  virtual void markFaces(const SimulationDomain2 *domain);
  virtual void markCells(const SimulationDomain2 *, const ParticleSystem2 *);
  virtual BoundaryType faceBoundaryType(size_t faceId) const;
  virtual MaterialType faceType(size_t faceId) const;
  virtual double faceValue(size_t faceId) const;
  virtual void setBoundaryFaceValue(size_t faceId, double v);
  virtual MaterialType cellType(size_t cellId) const;
  virtual void resize(size_t cellCount, size_t faceCount);
  virtual size_t dirichletBoundaryCount() const;
  virtual unsigned fluidCellCount() const;
  friend std::ostream &operator<<(std::ostream &o, BoundaryType b) {
    switch (b) {
    case BoundaryType::NONE:
      o << "NONE";
      break;
    case BoundaryType::CUSTOM:
      o << "CUSTOM";
      break;
    case BoundaryType::DIRICHLET:
      o << "DIRICHLET";
      break;
    case BoundaryType::NEUMANN:
      o << "NEUMANN";
      break;
    case BoundaryType::ROBIN:
      o << "ROBIN";
      break;
    }
    return o;
  }
  friend std::ostream &operator<<(std::ostream &o, MaterialType b) {
    switch (b) {
    case MaterialType::SOLID:
      o << "SOLID";
      break;
    case MaterialType::FLUID:
      o << "FLUID";
      break;
    case MaterialType::AIR:
      o << "AIR";
      break;
    case MaterialType::CUSTOM:
      o << "CUSTOM";
      break;
    }
    return o;
  }

  void printFaces();
  void iterateBoundaryFaces(
      const std::function<void(size_t, BoundaryType)> &f) const;
  void iterateBoundaryCells(const SimulationDomain2 *, MaterialType,
                            const std::function<void(size_t)> &f) const;
  Intersection intersect(const SimulationDomain2 *domain,
                         Point2d initialPosition, Point2d finalPosition) const;
  void setFaceType(unsigned int, BoundaryType);
  void setCellType(unsigned int, MaterialType);

private:
  BoundaryType _boundaries[4]; //!< left | bottom | right | top
  std::vector<BoundaryType> _faceBoundaries;
  std::vector<MaterialType> _cells, _faces;
  std::vector<double> _boundaryFacesValues;
};
} // namespace furoo

#endif // FUROO_PHYSICS_DOMAIN_BOUNDARY_HANDLER_H
