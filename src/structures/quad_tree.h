// DEPRECATED
#ifndef FUROO_STRUCTURES_QUAD_TREE_H
#define FUROO_STRUCTURES_QUAD_TREE_H

#include <functional>
#include <geometry/bbox.h>
#include <queue>

namespace furoo {

class QuadTree {
public:
  struct Node {
    friend class QuadTree;
    Node(const BBox2D &r, Node *p = nullptr) : _bbox(r), _level(0), _parent(p) {
      for (int i = 0; i < 4; i++)
        children[i] = nullptr;
      if (_parent)
        _level = _parent->_level + 1;
      id = 0;
    }
    bool isLeaf() const {
      return children[0] == nullptr && children[1] == nullptr &&
             children[2] == nullptr && children[3] == nullptr;
    }
    size_t level() const { return _level; }
    BBox2D region() const { return _bbox; }
    size_t leafIndex() const { return _leafIndex; }
    size_t id; //!< Node number
    Node *children[4];

  private:
    BBox2D _bbox;  //!< Node region in space
    size_t _level; //!< Tree level
    Node *_parent; //!< Pointer to parent
    //  | 0 | 1 |
    //  | 2 | 3 |
    size_t _childNumber; //!< Child number
    size_t _leafIndex;   //!< Leaf index
  };
  struct Leaf {
    Leaf(Node *n) : node(n) {}
    std::vector<size_t> edges;
    Node *node;
  };
  struct Edge {
    enum class NeighborOrientation { TOP, BOTTOM, LEFT, RIGHT };
    Edge(int a, NeighborOrientation ap, int b, NeighborOrientation bp,
         Point2d p) {
      v[0] = a;
      np[0] = ap;
      v[1] = b;
      np[1] = bp;
      faceCenterPosition = p;
      // ws[0] = ws[1] = 0.0;
      // isDirichilet = false;
    }
    /*double leafWeight(int i) const { return v[0] == i ? ws[0] : ws[1]; }
    double &leafWeight(int i) { return v[0] == i ? ws[0] : ws[1]; }
    double neighborWeight(int i) const { return v[0] == i ? ws[1] : ws[0]; }
    double &neighborWeight(int i) { return v[0] == i ? ws[1] : ws[0]; }*/
    int neighborIndex(int i) const { return v[0] == i ? v[1] : v[0]; }
    NeighborOrientation neighborOrientation(int i) const {
      return v[0] == i ? np[1] : np[0];
    }
    int v[2]; //!< leaf index
    // double ws[2]; //!< weights
    NeighborOrientation np[2];
    Point2d faceCenterPosition;
  };

  QuadTree();
  QuadTree(const BBox2D &region, const std::function<bool(Node &node)> &f,
           bool buildNeighborhood = true);
  QuadTree(const BBox2D &region, size_t maxLevel,
           bool buildNeighborhood = true);
  virtual ~QuadTree();
  BBox2D domainRegion() const { return _root->region(); }
  size_t height() const { return _height; }
  size_t nodeCount() const { return _count; }
  size_t leafCount() const { return _leafs.size(); }
  size_t edgeCount() const { return _leafEdges.size(); }
  const std::vector<Leaf> &leafs() const { return _leafs; }
  const std::vector<Edge> &leafEdges() const { return _leafEdges; }
  void traverse(const std::function<bool(Node &node)> &f);
  void traverse(const std::function<bool(const Node &node)> &f) const;
  void iterateLeafs(const std::function<void(Leaf &)> &f);
  int neighborIndex(size_t i, size_t e) const;
  void refine(const std::function<bool(Node &)> &f,
              const std::function<bool(Node &)> &sf);
  int intersect(Point2d) const;

protected:
  typename Edge::NeighborOrientation neighborOrientation(size_t i,
                                                         size_t e) const;
  Point2d neighborCenter(size_t i, size_t e) const;
  void traverse(Node *node, const std::function<bool(Node &node)> &f);
  void traverse(Node *node,
                const std::function<bool(const Node &node)> &f) const;
  void clear(Node *node);
  void refine(Node *node, const std::function<bool(Node &node)> &f,
              std::function<void(Node &)> sf = nullptr);
  void buildLeafNeighborhood();
  void buildLeafPhantomNeighborhood();

  size_t _height;
  size_t _count;
  Node *_root;
  std::vector<Leaf> _leafs;
  std::vector<Edge> _leafEdges;
};

} // namespace furoo

#endif
