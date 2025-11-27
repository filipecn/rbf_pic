#ifndef FUROO_STRUCTURES_POINT_Z_GRID_H
#define FUROO_STRUCTURES_POINT_Z_GRID_H

#include <common/debug.h>
#include <common/search.h>
#include <common/timer.h>
#include <functional>
#include <geometry/numeric.h>
#include <geometry/point.h>
#include <geometry/queries.h>
#include <structures/point_set.h>
#include <structures/quad_tree.h>
#include <structures/regular_grid.h>
#include <vector>

namespace furoo {

/// Arranges points following morton code.
/// \tparam T
/// \tparam D
template <typename T, int D> class PointZGrid : public PointSetInterface<T, D> {
public:
  friend class Iterator;
  struct PointElement {
    PointElement() : id(0), zcode(0) { active = true; }
    size_t id;    ///< element id for reference from outside code
    size_t zcode; ///< computed morton code
    bool active;  ///< element existence
  };
  /// Helper class for iterating elements
  class Iterator {
  public:
    /// \param s z point set reference
    /// \param f first element id (morton code)
    /// \param depth tree level
    explicit Iterator(PointZGrid<T, D> &s, size_t f = 0, size_t depth = 0)
        : _lastIndex(-1), _first(f), _depth(depth), _cur(0), _upperBound(0),
          _zps(s) {
      if (_zps._end == 0)
        return;
      _firstIndex = find(s, 0, s._end, _first, depth, &_upperBound);
      _cur = _firstIndex;
      THROW(_lastIndex <= static_cast<int>(_zps._end),
            "PointZGrid<>::iterator invalid lastIndex");
    }
    /// Given a range of indices from points_ array, finds the first element
    /// obeying f and depth restrictions. It performs a lower_bound operation,
    /// so the element may not be valid.
    /// \param z z point set reference
    /// \param l first index of points_ array \param s size of range on
    /// points_array \param f first element id (morton code) \param depth tree
    /// level \param upperBound **[out | optional]** upper bound morton code for
    /// this search
    static size_t find(PointZGrid &z, size_t l, size_t s, size_t f,
                       size_t depth, size_t *upperBound = nullptr) {
      std::function<int(const PointElement &p, const size_t &v)> comp;
      comp = [](const PointElement &p, const size_t &v) {
        if (p.zcode < v)
          return -1;
        if (p.zcode > v)
          return 1;
        return 0;
      };
      THROW(z._end >= l + s, "PointZGrid<>::Iterator invalid element index");
      if (upperBound)
        *upperBound = f + (1 << ((z._nbits - depth) * D));
      THROW(depth <= z._maxDepth, "PointZGrid<>::Iterator depth > maxDepth");
      return lower_bound<PointElement, size_t>(&z._points[l], s, f, comp) + l +
             1;
    }
    /// \return true if there is more elements to iterate
    bool next() const {
      return _cur < static_cast<int>(_zps._end) &&
             _zps._points[_cur].zcode < _upperBound;
    }
    /// \return position of the current element in world coordinates
    Point<T, D> getWorldPosition() {
      THROW(_cur < static_cast<int>(_zps._points.size()), "");
      THROW(_zps._points[_cur].id < _zps._positions.size(), "");
      return _zps._positions[_zps._points[_cur].id];
    }
    /// sets a new position to current element
    /// \param p new position value
    void setPosition(const Point<T, D> &p) {
      _zps._positions[_zps._points[_cur].id] = p;
      Point<T, D> gp = _zps._toSet(p);
      _zps._points[_cur].zcode = computeIndex(gp);
    }
    /// \return id of current element
    size_t getId() {
      THROW(_cur < static_cast<int>(_zps._points.size()),
            "PointZGrid<>::iterator invalid current element");
      return _zps._points[_cur].id;
    }
    /// \return pointer to current point element struct
    PointElement *pointElement() {
      THROW(_cur < static_cast<int>(_zps._points.size()), "");
      return &_zps._points[_cur];
    }
    /// advance on iteration
    void operator++() { _cur++; }
    /// \return number of elements to iterate (first -> last)
    size_t count() {
      if (_lastIndex < 1) {
        _lastIndex = find(_zps, static_cast<size_t>(_firstIndex),
                          _zps._end - _firstIndex, _upperBound, _depth);
        THROW(_lastIndex >= _firstIndex,
              "PointZGrid<>::iterator invalid index for count()");
      }
      return static_cast<size_t>(_lastIndex - _firstIndex);
    }

  private:
    int _firstIndex;    ///< first element index
    int _lastIndex;     ///< last element index (lazily computed on count())
    size_t _first;      ///< morton code of first element
    size_t _depth;      ///< depth of node where points must reside
    int _cur;           ///< current element index
    size_t _upperBound; ///< upper bound for search (morton code)
    PointZGrid &_zps;
  };
  /// Default constructor
  PointZGrid();
  /// Default destructor
  ~PointZGrid() override;
  /// Constructor
  /// \param maxCoordinates maximum coordinates value for an input point
  explicit PointZGrid(size_t maxCoordinates);
  /// necessary before any search after position modification
  void update() override;
  // INTERFACE
  void clear() override;
  size_t size() override;
  void setDomainRegion(const BBox<T, D> &region) override;
  size_t add(Point<T, D> p) override;
  void setPosition(size_t i, Point<T, D> p) override;
  void remove(size_t i) override;
  Point<T, D> operator[](size_t i) const override;
  void search(const BBox<T, D> &b,
              const std::function<void(size_t)> &f) override;
  void searchZ(const BBox<T, D> &b, size_t level,
               const std::function<void(size_t)> &f) override;
  void iteratePoints(
    const std::function<void(size_t, Point<T, D>)> &f) const override;
#ifdef _USE_OPENMP
  void iteratePoints_par(
    const std::function<void(size_t, Point<T, D>)> &f) const override;
#endif // _USE_OPENMP
  void iterateClosestPoints(Point<T, D> center, size_t n,
                            const std::function<void(size_t, Point<T, D>)> &f,
                            std::function<bool(size_t)> isValid) override;
  bool containsPoint(const BBox<T, D> &b) override;
  bool containsPointZ(const BBox<T, D> &b, size_t level) override;

private:
  static size_t computeIndex(Point<T, D> p);
  Point<T, D> _resolution;             ///< maximum coordinates
  Transform<D> _toSet;                 ///< map to underling grid
  std::vector<PointElement> _points;   ///< array of point elements
  std::vector<Point<T, D>> _positions; ///< points positions
  std::vector<size_t> _indices;        ///< map point id -> points_
  size_t _end;      ///< current number of active point element
  size_t _lastId;   ///< last point element id generated
  size_t _nbits;    ///< number of bits used by the maximum coordinate
  size_t _maxDepth; ///< maximum depth of search tree, bounded by nbits
  size_t _maxZCode; ///< maximum allowed morton code
  bool _needUpdate;
  bool _needZUpdate;
  bool _sizeChanged;
};

#include "structures/point_z_grid.inl"

typedef PointZGrid<double, 2> PointZGrid2d;
typedef PointZGrid<double, 3> PointZGrid3d;

// TODO: deprecated
class PointZGrid2 : public PointSet2Interface {
public:
  friend class PointIterator;
  friend class SearchTree;
  struct PointElement {
    PointElement() { active = true; }
    size_t id;
    size_t zcode;
    bool active;
  };

  class PointIterator {
  public:
    PointIterator(PointZGrid2 &g, size_t f = 0, size_t depth = 0)
        : first(f), cur(0), last(0), grid(g) {
      comp = [](const PointElement &p, const size_t &v) {
        if (p.zcode < v)
          return -1;
        if (p.zcode > v)
          return 1;
        return 0;
      };
      if (grid._end == 0)
        return;

      int fi = lower_bound<PointElement, size_t>(&grid._points[0], grid._end,
                                                 first, comp);
      last = grid._end;
      if (depth) {
        last = first + (1 << ((grid._nbits - depth) * 2));
        last = lower_bound<PointElement, size_t>(&grid._points[0], grid._end,
                                                 last, comp) +
               1;
      }
      first = fi + 1;
      cur = first;
    }
    bool next() const { return cur < last; }
    Point2d getWorldPosition() {
      ASSERT_FATAL(cur < last);
      return grid._positions[grid._points[cur].id];
    }
    void setPosition(const Point2d &p) {
      grid._positions[grid._points[cur].id] = p;
      Point2d gp = grid._grid.transformToGrid(p);
      grid._points[cur].zcode = mortonCode(gp.x(), gp.y());
    }
    size_t getId() {
      ASSERT_FATAL(cur < last);
      return grid._points[cur].id;
    }
    PointElement *pointElement() { return &grid._points[cur]; }

    void operator++() { cur++; }
    int count() const { return last - first; }

  private:
    size_t first, cur, last;
    PointZGrid2 &grid;
    std::function<int(const PointElement &p, const size_t &v)> comp;
  };

  class SearchTree : public QuadTree {
  public:
    SearchTree(PointZGrid2 &g,
               std::function<bool(size_t id, size_t depth)> f =
                   [](size_t id, size_t depth) -> bool {
                 UNUSED_VARIABLE(id);
                 UNUSED_VARIABLE(depth);
                 return true;
               })
        : grid(g) {
      this->_root = new QuadTree::Node(
          BBox2D(Point2d(), Point2d(grid._grid.gridSize().x(),
                                    grid._grid.gridSize().y())));
      this->_count++;
      ids.emplace_back(0);
      this->refine(
          this->_root,
          [&](QuadTree::Node &node) -> bool {
            if (node.level() >= grid._maxDepth)
              return false;
            if (!f(ids[node.id], node.level()))
              return false;
            return true;
          },
          [&](QuadTree::Node &node) {
            ids.emplace_back(ids[node.id] |
                             (2 << ((grid._nbits - node.level() - 1) * 2)));
            ids.emplace_back(ids[node.id] |
                             (3 << ((grid._nbits - node.level() - 1) * 2)));
            ids.emplace_back(ids[node.id] |
                             (0 << ((grid._nbits - node.level() - 1) * 2)));
            ids.emplace_back(ids[node.id] |
                             (1 << ((grid._nbits - node.level() - 1) * 2)));
            return true;
          });
    }
    ~SearchTree() {}
    void iteratePoints(const BBox2D &bbox,
                       std::function<void(PointElement *o)> f) {
      BBox2D gbbox = grid._grid.getGridTransform()(bbox);
      traverse([&](QuadTree::Node &node) -> bool {
        if (bbox_bbox_intersection(gbbox, node.region())) {
          if (!node.children[0] || !node.children[1] || !node.children[2] ||
              !node.children[3]) {
            PointIterator it(grid, ids[node.id], node.level());
            while (it.next()) {
              if (bbox.inside(it.getWorldPosition()))
                f(it.pointElement());
              ++it;
            }
          } else
            return true;
        }
        return false;
      });
    }

    std::vector<size_t> ids;

  private:
    PointZGrid2 &grid;
  };

  PointZGrid2();
  template <typename... Args>
  PointZGrid2(Args &&... args)
      : _grid(RegularGrid(std::forward<Args>(args)...)) {
    _end = _lastId = 0;
    _tree = nullptr;
    _needUpdate = false;
    _sizeChanged = false;
    set();
  }
  virtual ~PointZGrid2();

  void setDomainRegion(const BBox2D &r) override {
    _grid = RegularGrid(64, 64, r);
    set();
  }
  size_t size() override;
  size_t add(Point2d p) override;
  void search(const BBox2D &b,
              const std::function<void(size_t id)> &f) override;
  Point2d operator[](size_t i) const override;
  void update();
  size_t nbits() const { return _nbits; }
  size_t maxDepth() const { return _maxDepth; }
  SearchTree &searchTree() const { return *_tree; }
  void setPosition(unsigned int i, Point2d p) override;
  void remove(unsigned int i) override;
  void
  iteratePoints(const std::function<void(size_t, Point2d)> &f) const override;
  void iterateClosestPoints(Point2d, size_t,
                            const std::function<void(size_t, Point2d)> &,
                            std::function<bool(size_t)>) override;

private:
  void set();

  RegularGrid _grid;
  SearchTree *_tree;
  size_t _end;
  size_t _lastId;
  std::vector<PointElement> _points;
  std::vector<Point2d> _positions;
  std::vector<unsigned int> _indices;

  size_t _nbits, _maxDepth, _maxZCode;
  bool _needUpdate;
  bool _needZUpdate;
  bool _sizeChanged;

  // position in index space
  // position in grid coordinates
  static size_t _computeIndex(const Point2d &p) {
    return mortonCode(p.x(), p.y());
  }
};

} // namespace furoo

#endif // FUROO_STRUCTURES_POINT_Z_GRID_H
