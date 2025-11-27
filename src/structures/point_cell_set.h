#ifndef FUROO_POINT_CELL_SET_H
#define FUROO_POINT_CELL_SET_H

#include <structures/point_set.h>

namespace furoo {

/// Arranges points following morton code.
/// \tparam T
/// \tparam D
template<typename T, int D>
class PointCellSet : public PointSetInterface<T, D> {
 public:
  friend class Iterator;
  struct PointElement {
    PointElement() : id(0), cellId(0) { active = true; }
    size_t id;  ///< element id for reference from outside code
    size_t cellId; ///< point cell id
    bool active;  ///< element existence
  };
  /// Default constructor
  PointCellSet();
  /// Default destructor
  ~PointCellSet() override;
  /// necessary before any search after position modification
  void update();
  void setCell(size_t pointId, size_t cellId);
  void getCell(size_t pointId) const;
  // INTERFACE

  size_t size() override;
  void setDomainRegion(const BBox<T, D> &region) override;
  size_t add(Point<T, D> p) override;
  void setPosition(size_t i, Point<T, D> p) override;
  void remove(size_t i) override;
  Point<T, D> operator[](size_t i) const override;
  void search(const BBox<T, D> &b, const std::function<void(size_t)> &f)
  override;
  void iteratePoints(const std::function<void(size_t,
                                              Point<T,
                                                    D>)> &f) const
  override;
  void iterateClosestPoints(Point<T, D> center,
                            size_t n,
                            const std::function<void(size_t,
                                                     Point<T, D>)> &f,
                            std::function<bool(size_t)> isValid) override;
  void iterateCellPoints(size_t cellId,
                         const std::function<void(size_t, Point<T, D>)> &f);
 private:
  std::vector<PointElement> _points;  ///< array of point elements
  std::vector<Point<T, D>> _positions;  ///< points positions
  std::vector<size_t> _indices;        ///< map point id -> points_
  size_t _end;          ///< current number of active point element
  size_t _lastId;       ///< last point element id generated
  bool _needsUpdate;
};

#include "point_cell_set.inl"

}

#endif //FUROO_POINT_CELL_SET_H
