#ifndef FUROO_STRUCTURES_POINT_SET_H
#define FUROO_STRUCTURES_POINT_SET_H

#include <geometry/bbox.h>

namespace furoo {

template <typename T, int D> class PointSetInterface {
public:
  PointSetInterface() = default;

  virtual ~PointSetInterface() = default;

  /// \param i point id
  /// \return point i position
  virtual Point<T, D> operator[](size_t i) const = 0;

  /// \param r region where points points are allowed to exist
  virtual void setDomainRegion(const BBox<T, D> &r) = 0;

  /// remove all active points (doesn't free any memory)
  virtual void clear() = 0;

  /// \return number of active points
  virtual size_t size() = 0;

  /// \param p position of new point
  /// \return add new point
  virtual size_t add(Point<T, D> p) = 0;

  /// \param i point index
  /// \param position new position
  virtual void setPosition(size_t i, Point<T, D> position) = 0;

  /// \param i point index
  virtual void remove(size_t i) = 0;

  /// update this point set (for OpenMP, need to be fixed)
  virtual void update() = 0;

  /// \param b region of search
  /// \param f callback for each point found inside b
  virtual void search(const BBox<T, D> &b,
                      const std::function<void(size_t)> &f) = 0;

  virtual void searchZ(const BBox<T, D> &b, size_t level,
                       const std::function<void(size_t)> &f) = 0;

  /// Check the existence of a point inside a search region
  /// \param b search region
  /// \return true if at least one point is foun
  virtual bool containsPoint(const BBox<T, D> &b) = 0;
  virtual bool containsPointZ(const BBox<T, D> &b, size_t level) = 0;

  /// \param f callback for each point
  virtual void iteratePoints(
    const std::function<void(size_t, Point<T, D>)> &f) const = 0;

#ifdef _USE_OPENMP
  /// Iterate over points in parallel using OpenMP
  /// \param f callback (no race conditions) for each point
  virtual void iteratePoints_par(
    const std::function<void(size_t, Point<T, D>)> &f) const = 0;
#endif // _USE_OPENMP

  /// iterates n closest points to center filtered by isValid function
  /// \param center query point
  /// \param n number of nearest points
  /// \param f callback for processing each point
  /// \param isValid [optional] returns true if points is valid
  virtual void iterateClosestPoints(Point<T, D> center,
    size_t n,
    const std::function<void(size_t, Point<T, D>)> &f,
    std::function<bool(size_t)> isValid) = 0;
};

class PointSet2Interface {
public:
  PointSet2Interface() {}
  virtual ~PointSet2Interface() {}
  virtual Point2d operator[](size_t) const = 0;
  virtual void setDomainRegion(const BBox2D &) = 0;
  virtual size_t size() = 0;
  virtual size_t add(Point2d) = 0;
  virtual void setPosition(unsigned int, Point2d) = 0;
  virtual void remove(unsigned int) = 0;
  virtual void search(const BBox2D &, const std::function<void(size_t)> &) = 0;
  virtual void
  iteratePoints(const std::function<void(size_t, Point2d)> &f) const = 0;
  /// iterates n closest points to center filtered by isValid function
  /// \param center query point
  /// \param n number of nearest points
  /// \param f callback for processing each point
  /// \param isValid [optional] returns true if points is valid
  virtual void
  iterateClosestPoints(Point2d center, size_t n,
                       const std::function<void(size_t, Point2d)> &f,
                       std::function<bool(size_t)> isValid) = 0;
};

} // namespace furoo

#endif // FUROO_STRUCTURES_POINT_SET_H
