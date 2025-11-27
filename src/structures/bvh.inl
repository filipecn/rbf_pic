template <int D> BVH<D>::BVH(Mesh<D> &mesh) : _mesh(mesh) {
  std::vector<BVHElement> buildData;
  for (size_t i = 0; i < _mesh.elementCount(); i++)
    buildData.emplace_back(BVHElement(i, _mesh.elementBBox(i)));
  uint32_t totalNodes = 0;
  orderedElements.reserve(_mesh.elementCount());
  root = recursiveBuild(buildData, 0, _mesh.elementCount(), &totalNodes,
                        orderedElements);
  nodes.resize(totalNodes);
  for (uint32_t i = 0; i < totalNodes; i++)
    new (&nodes[i]) LinearBVHNode;
  uint32_t offset = 0;
  flattenBVHTree(root, &offset);
}

template <int D>
typename BVH<D>::BVHNode *
BVH<D>::recursiveBuild(std::vector<BVHElement> &buildData, uint32_t start,
                       uint32_t end, uint32_t *totalNodes,
                       std::vector<uint32_t> &orderedElements) {
  (*totalNodes)++;
  BVHNode *node = new BVHNode();
  BBox<double, D> bbox;
  for (uint32_t i = start; i < end; ++i)
    bbox = BBox<double, D>::combine(bbox, buildData[i].bounds);
  // compute all bounds
  uint32_t nElements = end - start;
  if (nElements == 1) {
    // create leaf node
    uint32_t firstElementOffset = orderedElements.size();
    for (uint32_t i = start; i < end; i++) {
      uint32_t elementNum = buildData[i].ind;
      orderedElements.emplace_back(elementNum);
    }
    node->initLeaf(firstElementOffset, nElements, bbox);
  } else {
    // compute bound of primitives
    BBox<double, D> centroidBounds;
    for (uint32_t i = start; i < end; i++)
      centroidBounds =
          BBox<double, D>::combine(centroidBounds, buildData[i].centroid);
    auto bsize = centroidBounds.size();
    int dim = 0;
    for (int d = 0; d < D; d++)
      if (bsize[d] > bsize[dim])
        dim = d;
    // partition primitives
    uint32_t mid = (start + end) / 2;
    if (centroidBounds.upper()[dim] == centroidBounds.lower()[dim]) {
      node->initInterior(
          dim,
          recursiveBuild(buildData, start, mid, totalNodes, orderedElements),
          recursiveBuild(buildData, mid, end, totalNodes, orderedElements));
      return node;
    }
    // partition into equally sized subsets
    std::nth_element(&buildData[start], &buildData[mid],
                     &buildData[end - 1] + 1, ComparePoints(dim));
    node->initInterior(
        dim, recursiveBuild(buildData, start, mid, totalNodes, orderedElements),
        recursiveBuild(buildData, mid, end, totalNodes, orderedElements));
  }
  return node;
}

template <int D>
uint32_t BVH<D>::flattenBVHTree(BVHNode *node, uint32_t *offset) {
  LinearBVHNode *linearNode = &nodes[*offset];
  linearNode->bounds = node->bounds;
  uint32_t myOffset = (*offset)++;
  if (node->nElements > 0) {
    linearNode->elementsOffset = node->firstElementOffset;
    linearNode->nElements = node->nElements;
  } else {
    linearNode->axis = node->splitAxis;
    linearNode->nElements = 0;
    flattenBVHTree(node->children[0], offset);
    linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
  }
  return myOffset;
}

template <int D>
int BVH<D>::intersect(const Point<double, D> &origin,
                      const Vector<double, D> &direction, float *t) {
  UNUSED_VARIABLE(t);
  if (!nodes.size())
    return false;
  //  ponos::Transform inv = ponos::inverse(sceneMesh->transform);
  //  ponos::Ray3 r = inv(ray);
  int hit = 0;
  Vector<double, D> invDir;
  uint32_t dirIsNeg[D];
  for (int d = 0; d < D; d++) {
    invDir[d] = 1. / direction[d];
    dirIsNeg[d] = invDir[d] < 0;
  }
  uint32_t todoOffset = 0, nodeNum = 0;
  uint32_t todo[64];
  while (true) {
    LinearBVHNode *node = &nodes[nodeNum];
    if (intersect(node->bounds, origin, direction)) {
      if (node->nElements > 0) {
        // intersect ray with primitives
        for (uint32_t i = 0; i < node->nElements; i++) {
          if (_mesh.intersectElement(orderedElements[node->elementsOffset + i],
                                     origin, direction))
            hit++;
        }
        if (todoOffset == 0)
          break;
        nodeNum = todo[--todoOffset];
      } else {
        if (dirIsNeg[node->axis]) {
          todo[todoOffset++] = nodeNum + 1;
          nodeNum = node->secondChildOffset;
        } else {
          todo[todoOffset++] = node->secondChildOffset;
          nodeNum++;
        }
      }
    } else {
      if (todoOffset == 0)
        break;
      nodeNum = todo[--todoOffset];
    }
  }
  return hit;
}

template <int D>
bool BVH<D>::intersect(const BBox<double, D> &bounds,
                       const Point<double, D> &origin,
                       const Vector<double, D> &direction) const {
  return bounds.intersect(origin, direction);
}

template <int D> bool BVH<D>::isInside(const Point<double, D> &p) {
  return intersect(p, Vector<double, D>(1.2), nullptr) % 2 &&
         intersect(p, Vector<double, D>(-0.8), nullptr) % 2;
}
