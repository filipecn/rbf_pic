#ifndef FUROO_BVH_H
#define FUROO_BVH_H

#include <geometry/mesh.h>

namespace furoo {

template <int D> class BVH {
public:
  BVH(Mesh<D> &mesh);
  virtual ~BVH() {}

  int intersect(const Point<double, D> &origin,
                const Vector<double, D> &direction, float *t = nullptr);
  bool isInside(const Point<double, D> &p);

private:
  struct BVHElement {
    BVHElement(size_t i, const BBox<double, D> &b) : ind(i), bounds(b) {
      centroid = b.center();
    }
    size_t ind;
    BBox<double, D> bounds;
    Point<double, D> centroid;
  };
  struct BVHNode {
    BVHNode() { children[0] = children[1] = nullptr; }
    void initLeaf(uint32_t first, uint32_t n, const BBox<double, D> &b) {
      firstElementOffset = first;
      nElements = n;
      bounds = b;
    }
    void initInterior(uint32_t axis, BVHNode *c0, BVHNode *c1) {
      children[0] = c0;
      children[1] = c1;
      bounds = BBox<double, D>::combine(c0->bounds, c1->bounds);
      splitAxis = axis;
      nElements = 0;
    }
    BBox<double, D> bounds;
    BVHNode *children[2];
    uint32_t splitAxis = 0, firstElementOffset = 0, nElements = 0;
  };
  struct LinearBVHNode {
    BBox<double, D> bounds;
    union {
      uint32_t elementsOffset;
      uint32_t secondChildOffset;
    };
    uint8_t nElements;
    uint8_t axis;
    uint8_t pad[2];
  };
  struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHElement &a, const BVHElement &b) const {
      return a.centroid[dim] < b.centroid[dim];
    }
  };
  std::vector<uint32_t> orderedElements;
  std::vector<LinearBVHNode> nodes;
  BVHNode *root;
  BVHNode *recursiveBuild(std::vector<BVHElement> &buildData, uint32_t start,
                          uint32_t end, uint32_t *totalNodes,
                          std::vector<uint32_t> &orderedElements);
  uint32_t flattenBVHTree(BVHNode *node, uint32_t *offset);

  bool intersect(const BBox<double, D> &bounds, const Point<double, D> &origin,
                 const Vector<double, D> &direction) const;
  Mesh<D> &_mesh;
};

#include "bvh.inl"

} // namespace furoo

#endif // FUROO_BVH_H
