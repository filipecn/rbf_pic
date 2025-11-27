#include <algorithm>
#include <common/debug.h>
#include <geometry/point.h>
#include <structures/octree.h>

namespace furoo {

Octree::Octree() : _root(nullptr) {
  _height = 0;
  _count = 0;
}

Octree::Octree(const BBox3D &region, const std::function<bool(Node &node)> &f,
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

Octree::Octree(const BBox3D &region, unsigned maxLevel, bool buildNeighborhood)
    : Octree(region,
             [maxLevel](Node &node) -> bool {
               if (node._level < maxLevel)
                 return true;
               return false;
             },
             buildNeighborhood) {}

Octree::~Octree() {
  clear(_root);
  _height = _count = 0;
}

void Octree::traverse(const std::function<bool(Node &node)> &f) {
  traverse(_root, f);
}

void Octree::traverse(const std::function<bool(const Node &node)> &f) const {
  traverse(_root, f);
}

void Octree::traverse(Node *node, const std::function<bool(Node &node)> &f) {
  if (!node)
    return;
  if (!f(*node))
    return;
  for (auto c : node->children)
    traverse(c, f);
}

void Octree::traverse(Node *node,
                      const std::function<bool(const Node &node)> &f) const {
  if (!node)
    return;
  if (!f(*node))
    return;
  for (auto c : node->children)
    traverse(c, f);
}

void Octree::iterateLeafs(const std::function<void(Leaf &)> &f) {
  for (auto &leaf : _leafs)
    f(leaf);
}

int Octree::neighborIndex(unsigned i, unsigned e) const {
  return _leafEdges[_leafs[i].edges[e]].neighborIndex(i);
}

typename Octree::Edge::NeighborOrientation
Octree::neighborOrientation(unsigned i, unsigned e) const {
  return _leafEdges[_leafs[i].edges[e]].neighborOrientation(i);
}

Point3d Octree::neighborCenter(unsigned i, unsigned e) const {
  ASSERT_FATAL(neighborIndex(i, e) >= 0);
  return _leafs[_leafEdges[_leafs[i].edges[e]].neighborIndex(i)]
      .node->region()
      .center();
}

void Octree::refine(const std::function<bool(Node &)> &f,
                    const std::function<bool(Node &)> &sf) {
  refine(_root, f, sf);
}

int Octree::intersect(Point3d p) const {
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

void Octree::clear(Node *node) {
  if (!node)
    return;
  for (int i = 0; i < 8; i++)
    clear(node->children[i]);
  delete node;
}

void Octree::refine(Node *node, const std::function<bool(Node &node)> &f,
                    std::function<void(Node &)> sf) {
  if (!node)
    return;
  _height = std::max(_height, node->_level);
  if (f(*node)) {
    Point3d p[3] = {node->_bbox.lower(), node->_bbox.center(),
                    node->_bbox.upper()};
    int k = 0;
    for (int z = 0; z < 2; z++)
      for (int y = 0; y < 2; y++)
        for (int x = 0; x < 2; x++) {
          node->children[k] =
              new Node(BBox3D(Point3d(p[x][0], p[y][1], p[z][2]),
                              Point3d(p[x + 1][0], p[y + 1][1], p[z + 1][2])),
                       node);
          node->children[k++]->id = _count++;
        }
    if (sf)
      sf(*node);
    for (int i = 0; i < 8; i++) {
      node->children[i]->_childNumber = i;
      refine(node->children[i], f, sf);
    }
  } else {
    node->_leafIndex = _leafs.size();
    _leafs.emplace_back(node);
  }
}

void Octree::buildLeafNeighborhood() {
  std::function<void(Node *, Node *)> xProcess = [&](Node *left, Node *right) {
    if (!left || !right)
      return;
    if (left->isLeaf() && right->isLeaf()) {
      _leafs[left->leafIndex()].edges.emplace_back(_leafEdges.size());
      _leafs[right->leafIndex()].edges.emplace_back(_leafEdges.size());
      Point3d faceCenterPosition;
      if (left->level() >= right->level())
        faceCenterPosition = left->region().center() +
                             Vector3d(left->region().size(0), 0., 0.) * .5;
      else
        faceCenterPosition = right->region().center() -
                             Vector3d(right->region().size(0), 0., 0.) * .5;
      _leafEdges.emplace_back(
          left->leafIndex(), Edge::NeighborOrientation::LEFT,
          right->leafIndex(), Edge::NeighborOrientation::RIGHT,
          faceCenterPosition);
    } else if (left->isLeaf()) {
      xProcess(left, right->children[0]);
      xProcess(left, right->children[2]);
      xProcess(left, right->children[4]);
      xProcess(left, right->children[6]);
    } else if (right->isLeaf()) {
      xProcess(left->children[1], right);
      xProcess(left->children[3], right);
      xProcess(left->children[5], right);
      xProcess(left->children[7], right);
    } else {
      xProcess(left->children[1], right->children[0]);
      xProcess(left->children[3], right->children[2]);
      xProcess(left->children[5], right->children[4]);
      xProcess(left->children[7], right->children[6]);
    }
  };
  std::function<void(Node *, Node *)> yProcess = [&](Node *top, Node *bottom) {
    if (!top || !bottom)
      return;
    if (top->isLeaf() && bottom->isLeaf()) {
      _leafs[top->leafIndex()].edges.emplace_back(_leafEdges.size());
      _leafs[bottom->leafIndex()].edges.emplace_back(_leafEdges.size());
      Point3d faceCenterPosition;
      if (top->level() >= bottom->level())
        faceCenterPosition = top->region().center() -
                             Vector3d(0., top->region().size(1), 0.) * .5;
      else
        faceCenterPosition = bottom->region().center() +
                             Vector3d(0., bottom->region().size(1), 0.) * .5;

      _leafEdges.emplace_back(
          top->leafIndex(), Edge::NeighborOrientation::TOP, bottom->leafIndex(),
          Edge::NeighborOrientation::BOTTOM, faceCenterPosition);
    } else if (top->isLeaf()) {
      yProcess(top, bottom->children[0]);
      yProcess(top, bottom->children[1]);
      yProcess(top, bottom->children[4]);
      yProcess(top, bottom->children[5]);
    } else if (bottom->isLeaf()) {
      yProcess(top->children[2], bottom);
      yProcess(top->children[3], bottom);
      yProcess(top->children[6], bottom);
      yProcess(top->children[7], bottom);
    } else {
      yProcess(top->children[0], bottom->children[2]);
      yProcess(top->children[1], bottom->children[3]);
      yProcess(top->children[4], bottom->children[6]);
      yProcess(top->children[5], bottom->children[7]);
    }
  };
  std::function<void(Node *, Node *)> zProcess = [&](Node *front, Node *back) {
    if (!front || !back)
      return;
    if (front->isLeaf() && back->isLeaf()) {
      _leafs[front->leafIndex()].edges.emplace_back(_leafEdges.size());
      _leafs[back->leafIndex()].edges.emplace_back(_leafEdges.size());
      Point3d faceCenterPosition;
      if (front->level() >= back->level())
        faceCenterPosition = front->region().center() -
                             Vector3d(0., 0., front->region().size(2)) * .5;
      else
        faceCenterPosition = back->region().center() +
                             Vector3d(0., 0., back->region().size(2)) * .5;

      _leafEdges.emplace_back(
          front->leafIndex(), Edge::NeighborOrientation::FRONT,
          back->leafIndex(), Edge::NeighborOrientation::BACK,
          faceCenterPosition);
    } else if (front->isLeaf()) {
      zProcess(front, back->children[4]);
      zProcess(front, back->children[5]);
      zProcess(front, back->children[6]);
      zProcess(front, back->children[7]);
    } else if (back->isLeaf()) {
      zProcess(front->children[0], back);
      zProcess(front->children[1], back);
      zProcess(front->children[2], back);
      zProcess(front->children[3], back);
    } else {
      zProcess(front->children[0], back->children[4]);
      zProcess(front->children[1], back->children[5]);
      zProcess(front->children[2], back->children[6]);
      zProcess(front->children[3], back->children[7]);
    }
  };
  std::function<void(Node *)> faceProcess = [&](Node *node) {
    if (!node)
      return;
    for (auto child : node->children)
      faceProcess(child);
    xProcess(node->children[0], node->children[1]);
    xProcess(node->children[2], node->children[3]);
    xProcess(node->children[4], node->children[5]);
    xProcess(node->children[6], node->children[7]);
    yProcess(node->children[2], node->children[0]);
    yProcess(node->children[3], node->children[1]);
    yProcess(node->children[6], node->children[4]);
    yProcess(node->children[7], node->children[5]);
    zProcess(node->children[4], node->children[0]);
    zProcess(node->children[5], node->children[1]);
    zProcess(node->children[6], node->children[2]);
    zProcess(node->children[7], node->children[3]);
  };
  faceProcess(_root);
}

void Octree::buildLeafPhantomNeighborhood() {
  std::function<void(Node *, Node *)> xProcess = [&](Node *left, Node *right) {
    if (!left && !right)
      return;
    if (left) {
      if (left->isLeaf()) {
        _leafs[left->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            left->region().center() +
            Vector3d(left->region().size(0), 0., 0.) * .5;
        _leafEdges.emplace_back(
            left->leafIndex(), Edge::NeighborOrientation::LEFT, -1,
            Edge::NeighborOrientation::RIGHT, faceCenterPosition);
      } else {
        xProcess(left->children[1], right);
        xProcess(left->children[3], right);
        xProcess(left->children[5], right);
        xProcess(left->children[7], right);
      }
    } else {
      if (right->isLeaf()) {
        _leafs[right->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            right->region().center() -
            Vector3d(right->region().size(0), 0., 0.) * .5;
        _leafEdges.emplace_back(
            right->leafIndex(), Edge::NeighborOrientation::RIGHT, -1,
            Edge::NeighborOrientation::LEFT, faceCenterPosition);
      } else {
        xProcess(left, right->children[0]);
        xProcess(left, right->children[2]);
        xProcess(left, right->children[4]);
        xProcess(left, right->children[6]);
      }
    }
  };
  xProcess(_root, nullptr);
  xProcess(nullptr, _root);
  std::function<void(Node *, Node *)> yProcess = [&](Node *top, Node *bottom) {
    if (!top && !bottom)
      return;
    if (top) {
      if (top->isLeaf()) {
        _leafs[top->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            top->region().center() -
            Vector3d(0., top->region().size(1), 0.) * .5;
        _leafEdges.emplace_back(
            top->leafIndex(), Edge::NeighborOrientation::TOP, -1,
            Edge::NeighborOrientation::BOTTOM, faceCenterPosition);
      } else {
        yProcess(top->children[0], bottom);
        yProcess(top->children[1], bottom);
        yProcess(top->children[4], bottom);
        yProcess(top->children[5], bottom);
      }
    } else {
      if (bottom->isLeaf()) {
        _leafs[bottom->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            bottom->region().center() +
            Vector3d(0., bottom->region().size(1), 0.) * .5;
        _leafEdges.emplace_back(
            bottom->leafIndex(), Edge::NeighborOrientation::BOTTOM, -1,
            Edge::NeighborOrientation::TOP, faceCenterPosition);
      } else {
        yProcess(top, bottom->children[2]);
        yProcess(top, bottom->children[3]);
        yProcess(top, bottom->children[6]);
        yProcess(top, bottom->children[7]);
      }
    }
  };
  yProcess(_root, nullptr);
  yProcess(nullptr, _root);
  std::function<void(Node *, Node *)> zProcess = [&](Node *front, Node *back) {
    if (!front && !back)
      return;
    if (front) {
      if (front->isLeaf()) {
        _leafs[front->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            front->region().center() -
            Vector3d(0., 0., front->region().size(1)) * .5;
        _leafEdges.emplace_back(
            front->leafIndex(), Edge::NeighborOrientation::FRONT, -1,
            Edge::NeighborOrientation::BACK, faceCenterPosition);
      } else {
        zProcess(front->children[0], back);
        zProcess(front->children[1], back);
        zProcess(front->children[2], back);
        zProcess(front->children[3], back);
      }
    } else {
      if (back->isLeaf()) {
        _leafs[back->leafIndex()].edges.emplace_back(_leafEdges.size());
        Point3d faceCenterPosition =
            back->region().center() +
            Vector3d(0., 0., back->region().size(1)) * .5;
        _leafEdges.emplace_back(
            back->leafIndex(), Edge::NeighborOrientation::BACK, -1,
            Edge::NeighborOrientation::FRONT, faceCenterPosition);
      } else {
        zProcess(front, back->children[4]);
        zProcess(front, back->children[5]);
        zProcess(front, back->children[6]);
        zProcess(front, back->children[7]);
      }
    }
  };
  zProcess(_root, nullptr);
  zProcess(nullptr, _root);
}

} // furoo namespace
