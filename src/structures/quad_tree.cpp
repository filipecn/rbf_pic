#include <algorithm>
#include <common/debug.h>
#include <geometry/point.h>
#include <structures/quad_tree.h>

namespace furoo {

QuadTree::QuadTree() : _root(nullptr) {
  _height = 0;
  _count = 0;
}

QuadTree::QuadTree(const BBox2D &region,
                   const std::function<bool(Node &node)> &f,
                   bool buildNeighborhood)
    : _root(new Node(region, nullptr)) {
  _height = 0;
  _count = 1;
  refine(_root, f, nullptr);
  if (buildNeighborhood) {
    buildLeafNeighborhood();
    buildLeafPhantomNeighborhood();
  }
}

QuadTree::QuadTree(const BBox2D &region, size_t maxLevel,
                   bool buildNeighborhood)
    : QuadTree(region,
               [maxLevel](Node &node) -> bool {
                 if (node._level < maxLevel)
                   return true;
                 return false;
               },
               buildNeighborhood) {}

QuadTree::~QuadTree() {
  clear(_root);
  _height = _count = 0;
}

void QuadTree::traverse(const std::function<bool(Node &node)> &f) {
  traverse(_root, f);
}

void QuadTree::traverse(const std::function<bool(const Node &node)> &f) const {
  traverse(_root, f);
}

void QuadTree::traverse(Node *node, const std::function<bool(Node &node)> &f) {
  if (!node)
    return;
  if (!f(*node))
    return;
  for (auto c : node->children)
    traverse(c, f);
}

void QuadTree::traverse(Node *node,
                        const std::function<bool(const Node &node)> &f) const {
  if (!node)
    return;
  if (!f(*node))
    return;
  for (auto c : node->children)
    traverse(c, f);
}

void QuadTree::iterateLeafs(const std::function<void(Leaf &)> &f) {
  for (auto &leaf : _leafs)
    f(leaf);
}
/*
QuadTree::Node *QuadTree::neighborNode(size_t i, size_t e) {
  return _leafs[_leafEdges[_leafs[i].edges[e]].neighborIndex(i)].node;
}*/

int QuadTree::neighborIndex(size_t i, size_t e) const {
  return _leafEdges[_leafs[i].edges[e]].neighborIndex(i);
}

typename QuadTree::Edge::NeighborOrientation
QuadTree::neighborOrientation(size_t i, size_t e) const {
  return _leafEdges[_leafs[i].edges[e]].neighborOrientation(i);
}

Point2d QuadTree::neighborCenter(size_t i, size_t e) const {
  ASSERT_FATAL(neighborIndex(i, e) >= 0);
  return _leafs[_leafEdges[_leafs[i].edges[e]].neighborIndex(i)]
      .node->region()
      .center();
}

void QuadTree::refine(const std::function<bool(Node &)> &f,
                      const std::function<bool(Node &)> &sf) {
  refine(_root, f, sf);
}

int QuadTree::intersect(Point2d p) const {
  int r = -1;
  traverse(_root, [&](const Node &n) {
    if (n.region().inside(p)) {
      r = n.leafIndex();
      return true;
    }
    return false;
  });
  return r;
}

void QuadTree::clear(Node *node) {
  if (!node)
    return;
  for (int i = 0; i < 4; i++)
    clear(node->children[i]);
  delete node;
}

void QuadTree::refine(Node *node, const std::function<bool(Node &node)> &f,
                      std::function<void(Node &)> sf) {
  if (!node)
    return;
  _height = std::max(_height, node->_level);
  if (f(*node)) {
    Point2d pmin = node->_bbox.lower();
    Point2d pmax = node->_bbox.upper();
    Point2d mid = node->_bbox.center();
    node->children[0] = new Node(
        BBox2D(Point2d(pmin.x(), mid.y()), Point2d(mid.x(), pmax.y())), node);
    node->children[0]->id = _count++;
    node->children[1] = new Node(BBox2D(mid, pmax), node);
    node->children[1]->id = _count++;
    node->children[2] = new Node(BBox2D(pmin, mid), node);
    node->children[2]->id = _count++;
    node->children[3] = new Node(
        BBox2D(Point2d(mid.x(), pmin.y()), Point2d(pmax.x(), mid.y())), node);
    node->children[3]->id = _count++;
    if (sf)
      sf(*node);
    for (int i = 0; i < 4; i++) {
      node->children[i]->_childNumber = i;
      refine(node->children[i], f, sf);
    }
  } else {
    node->_leafIndex = _leafs.size();
    _leafs.emplace_back(node);
  }
}

void QuadTree::buildLeafNeighborhood() {
  std::function<void(Node *, Node *)> horizontalProcess = [&](Node *left,
                                                              Node *right) {
    if (!left || !right)
      return;
    if (left->isLeaf() && right->isLeaf()) {
      _leafs[left->leafIndex()].edges.emplace_back(_leafEdges.size());
      _leafs[right->leafIndex()].edges.emplace_back(_leafEdges.size());
      Point2d faceCenterPosition;
      if (left->level() >= right->level())
        faceCenterPosition =
            left->region().center() + Vector2d(left->region().size(0), 0.) * .5;
      else
        faceCenterPosition = right->region().center() -
                             Vector2d(right->region().size(0), 0.) * .5;
      _leafEdges.emplace_back(
          left->leafIndex(), Edge::NeighborOrientation::LEFT,
          right->leafIndex(), Edge::NeighborOrientation::RIGHT,
          faceCenterPosition);
    } else if (left->isLeaf()) {
      horizontalProcess(left, right->children[0]);
      horizontalProcess(left, right->children[2]);
    } else if (right->isLeaf()) {
      horizontalProcess(left->children[1], right);
      horizontalProcess(left->children[3], right);
    } else {
      horizontalProcess(left->children[1], right->children[0]);
      horizontalProcess(left->children[3], right->children[2]);
    }
  };
  std::function<void(Node *, Node *)> verticalProcess = [&](Node *top,
                                                            Node *bottom) {
    if (!top || !bottom)
      return;
    if (top->isLeaf() && bottom->isLeaf()) {
      _leafs[top->leafIndex()].edges.emplace_back(_leafEdges.size());
      _leafs[bottom->leafIndex()].edges.emplace_back(_leafEdges.size());
      Point2d faceCenterPosition;
      if (top->level() >= bottom->level())
        faceCenterPosition =
            top->region().center() - Vector2d(0., top->region().size(1)) * .5;
      else
        faceCenterPosition = bottom->region().center() +
                             Vector2d(0., bottom->region().size(1)) * .5;

      _leafEdges.emplace_back(
          top->leafIndex(), Edge::NeighborOrientation::TOP, bottom->leafIndex(),
          Edge::NeighborOrientation::BOTTOM, faceCenterPosition);
    } else if (top->isLeaf()) {
      verticalProcess(top, bottom->children[0]);
      verticalProcess(top, bottom->children[1]);
    } else if (bottom->isLeaf()) {
      verticalProcess(top->children[2], bottom);
      verticalProcess(top->children[3], bottom);
    } else {
      verticalProcess(top->children[2], bottom->children[0]);
      verticalProcess(top->children[3], bottom->children[1]);
    }
  };
  std::function<void(Node *)> faceProcess = [&](Node *node) {
    if (!node)
      return;
    for (auto child : node->children)
      faceProcess(child);
    horizontalProcess(node->children[0], node->children[1]);
    horizontalProcess(node->children[2], node->children[3]);
    verticalProcess(node->children[0], node->children[2]);
    verticalProcess(node->children[1], node->children[3]);
  };
  faceProcess(_root);
}

void QuadTree::buildLeafPhantomNeighborhood() {
  std::function<void(Node *, Node *)> verticalProcess = [&](Node *top,
                                                            Node *bottom) {
    if (!top && !bottom)
      return;
    if (top) {
      if (top->isLeaf()) {
        _leafs[top->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point2d faceCenterPosition =
            top->region().center() - Vector2d(0., top->region().size(1)) * .5;
        _leafEdges.emplace_back(
            top->leafIndex(), Edge::NeighborOrientation::TOP, -1,
            Edge::NeighborOrientation::BOTTOM, faceCenterPosition);
      } else {
        verticalProcess(top->children[2], bottom);
        verticalProcess(top->children[3], bottom);
      }
    } else {
      if (bottom->isLeaf()) {
        _leafs[bottom->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point2d faceCenterPosition =
            bottom->region().center() +
            Vector2d(0., bottom->region().size(1)) * .5;
        _leafEdges.emplace_back(
            bottom->leafIndex(), Edge::NeighborOrientation::BOTTOM, -1,
            Edge::NeighborOrientation::TOP, faceCenterPosition);
      } else {
        verticalProcess(top, bottom->children[0]);
        verticalProcess(top, bottom->children[1]);
      }
    }
  };
  verticalProcess(_root, nullptr);
  verticalProcess(nullptr, _root);
  std::function<void(Node *, Node *)> horizontalProcess = [&](Node *left,
                                                              Node *right) {
    if (!left && !right)
      return;
    if (left) {
      if (left->isLeaf()) {
        _leafs[left->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point2d faceCenterPosition =
            left->region().center() + Vector2d(left->region().size(0), 0.) * .5;
        _leafEdges.emplace_back(
            left->leafIndex(), Edge::NeighborOrientation::LEFT, -1,
            Edge::NeighborOrientation::RIGHT, faceCenterPosition);
      } else {
        horizontalProcess(left->children[1], right);
        horizontalProcess(left->children[3], right);
      }
    } else {
      if (right->isLeaf()) {
        _leafs[right->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point2d faceCenterPosition = right->region().center() -
                                     Vector2d(right->region().size(0), 0.) * .5;
        _leafEdges.emplace_back(
            right->leafIndex(), Edge::NeighborOrientation::RIGHT, -1,
            Edge::NeighborOrientation::LEFT, faceCenterPosition);
      } else {
        horizontalProcess(left, right->children[0]);
        horizontalProcess(left, right->children[2]);
      }
    }
  };
  horizontalProcess(_root, nullptr);
  horizontalProcess(nullptr, _root);
}

} // furoo namespace
