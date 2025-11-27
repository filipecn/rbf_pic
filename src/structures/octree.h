#ifndef FUROO_STRUCTURES_OCTREE_H
#define FUROO_STRUCTURES_OCTREE_H

#include <functional>
#include <geometry/bbox.h>
#include <queue>

namespace furoo {

class Octree {
public:
  struct Node {
    friend class Octree;
    Node(const BBox3D &r, Node *p = nullptr) : _bbox(r), _level(0), _parent(p) {
      for (int i = 0; i < 8; i++)
        children[i] = nullptr;
      if (_parent)
        _level = _parent->_level + 1;
      id = 0;
    }
    bool isLeaf() const {
      return children[0] == nullptr && children[1] == nullptr &&
             children[2] == nullptr && children[3] == nullptr;
    }
    unsigned level() const { return _level; }
    BBox3D region() const { return _bbox; }
    unsigned leafIndex() const { return _leafIndex; }
    unsigned id; //!< Node number
    Node *children[8];

  private:
    BBox3D _bbox;    //!< Node region in space
    unsigned _level; //!< Tree level
    Node *_parent;   //!< Pointer to parent
    // NEIGHBORS
    //  y            3t     5ba
    //  | _x      0l    1r
    // z/            2b
    //           4f
    // CHILDREN
    //  y
    //  | _x
    // z/
    unsigned _childNumber; //!< Child number
    unsigned _leafIndex;   //!< Leaf index
  };
  struct Leaf {
    Leaf(Node *n) : node(n) {}
    std::vector<unsigned> edges;
    Node *node;
  };
  struct Edge {
    enum class NeighborOrientation { TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK };
    Edge(int a, NeighborOrientation ap, int b, NeighborOrientation bp,
         Point3d p) {
      v[0] = a;
      np[0] = ap;
      v[1] = b;
      np[1] = bp;
      faceCenterPosition = p;
    }
    int neighborIndex(int i) const { return v[0] == i ? v[1] : v[0]; }
    NeighborOrientation neighborOrientation(int i) const {
      return v[0] == i ? np[1] : np[0];
    }
    int v[2]; //!< leaf index
    NeighborOrientation np[2];
    Point3d faceCenterPosition;
  };

  Octree();
  Octree(const BBox3D &region, const std::function<bool(Node &node)> &f,
         bool buildNeighborhood = true);
  Octree(const BBox3D &region, unsigned maxLevel,
         bool buildNeighborhood = true);
  virtual ~Octree();
  BBox3D domainRegion() const { return _root->region(); }
  unsigned height() const { return _height; }
  unsigned nodeCount() const { return _count; }
  unsigned leafCount() const { return _leafs.size(); }
  unsigned edgeCount() const { return _leafEdges.size(); }
  const std::vector<Leaf> &leafs() const { return _leafs; }
  const std::vector<Edge> &leafEdges() const { return _leafEdges; }
  void traverse(const std::function<bool(Node &node)> &f);
  void traverse(const std::function<bool(const Node &node)> &f) const;
  void iterateLeafs(const std::function<void(Leaf &)> &f);
  int neighborIndex(unsigned i, unsigned e) const;
  void refine(const std::function<bool(Node &)> &f,
              const std::function<bool(Node &)> &sf);
  int intersect(Point3d) const;

protected:
  typename Edge::NeighborOrientation neighborOrientation(unsigned i,
                                                         unsigned e) const;
  Point3d neighborCenter(unsigned i, unsigned e) const;
  void traverse(Node *node, const std::function<bool(Node &node)> &f);
  void traverse(Node *node,
                const std::function<bool(const Node &node)> &f) const;
  void clear(Node *node);
  void refine(Node *node, const std::function<bool(Node &node)> &f,
              std::function<void(Node &)> sf = nullptr);
  void buildLeafNeighborhood();
  void buildLeafPhantomNeighborhood();

  unsigned _height;
  unsigned _count;
  Node *_root;
  std::vector<Leaf> _leafs;
  std::vector<Edge> _leafEdges;
};

} // namespace furoo

#endif
