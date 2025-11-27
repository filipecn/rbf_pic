#ifndef FUROO_CELL_GRAPH_H
#define FUROO_CELL_GRAPH_H

#include <blas/interpolation.h>
#include <common/definitions.h>
#include <functional>
#include <geometry/bbox.h>
#include <geometry/point.h>
#include <queue>
#include <set>
#include <structures/structure_interface.h>

namespace furoo {

/***
 * CellGraph structure Class
 * This structure works like a graph-based quadtree structure. Has min
 * and max levels, as in quadtree structures.
 * Used for adaptive applications. It is able to adapt according to refining
 * parameter/function. If min level and max level are the same, it works as
 * a slower regular grid.
 ***/
class CellGraph2 : public StructureInterface<2> {
public:
  class Edge;
  // ========================================================
  class Node {
  public:
    /// \param vertices node vertices in ccw order from lower left corner
    /// \param r node region
    /// \param id node id (accessed as cellId)
    /// \param level node level
    explicit Node(const size_t vertices[], const BBox2d &r, size_t id = 0,
                  size_t level = 0);

    /// Sets node data
    /// \param vertices node vertices in ccw order from lower left corner
    /// \param r node region
    /// \param level node tree level
    void set(const size_t vertices[], const BBox2d &r, size_t level);

    /// \return node's region
    BBox2d region() const;

    /// \return node's level
    size_t level() const;

    /// \param i vertex index (ccw starting from lower left corner)
    /// \return vertex id
    size_t vertex(size_t i) const;

    /// \return edges connected to the node
    const std::set<size_t> &edges() const;

    ///
    /// \return vertices in the node corner
    std::vector<size_t> vertices() const;

    /// \return edges connected to the node
    std::set<size_t> edgesCopy() const;

    /// Connect edge to node
    /// \param e edge index
    void addEdge(size_t e);

    /// Disconnect edge from node
    /// \param e edge index
    void removeEdge(size_t e);

    /// \return true if node id is valid
    bool isValid() const;

    /// \return node id
    size_t id() const;

    /// make invalid and clear all edges
    void destroy();

    Point2u addressCode() const;

    void setAddressCode(Point2u addressCode);

  private:
    std::set<size_t> _edges; //!< list of edge indices connected to this node
    BBox2d _region;          //!< node's region
    size_t _level;           //!< region's level of hierarchy
    size_t _id;              //!< node id (accessed as cellId)
    Point2u _addressCode; //!< relative unique id to maxLevel, used to compute
                          //!< siblings
    bool _isValid;        //!< true if node id is valid
    size_t _cornerVertices[4]{}; //!< list of vertex indices connected to this
                                 //!< node,
    //!< vertices are always stored in ccw order from lower left corner
  };
  // ========================================================
  class Edge {
  public:
    /// bSide is computed automatically
    /// \param aSide position of vertexA relative to edge center
    /// \param vertexA node a index
    /// \param vertexB node b index
    /// \param center face center position
    /// \param id edge id (faceId)
    Edge(Definitions::Side aSide, size_t vertexA, size_t vertexB,
         Point2d center, size_t id);

    /// bSide is computed automatically
    /// \param aSide position of vertexA relative to edge center
    /// \param vertexA node a index
    /// \param vertexB node b index
    /// \param center face center position
    void set(Definitions::Side aSide, size_t vertexA, size_t vertexB,
             Point2d center);

    /// \return edge's orientation
    Definitions::Orientation orientation() const;

    /// Given one of the edge vertices id, returns the side and id of the other
    /// vertex \param nodeId node id \return node neighbor [side and id]
    Definitions::Neighbor nodeNeighbor(size_t nodeId) const;

    /// \return true if node id is valid
    bool isValid() const;

    /// \return edge's id
    size_t id() const;

    /// makes edge invalid
    void destroy();

    std::vector<size_t> nodes() const;

    Point2d center() const;

  private:
    Definitions::Orientation _orientation =
        Definitions::Orientation::CUSTOM; //!< edge orientation (horizontal or
                                          //!< vertical)
    Definitions::Side
        _nodeSides[2]{};  //!< vertex position related to edge's center
    size_t _nodeIds[2]{}; //!< nodes connected by the edge
    Point2d _center;      //!< used as face center position
    bool _isValid = true;
    size_t _id;
  };
  // ========================================================
  class Vertex {
  public:
    Vertex(Point2d position, size_t id);

    void set(Point2d position, size_t id);

    void destroy();

    void removeCell();

    void addCell();

    size_t cellCount() const;

    Point2d position() const;

    bool isValid() const;

    size_t id() const;

  private:
    size_t _cellCount = 0;
    Point2d _position; //!< vertex location
    bool _isValid = true;
    size_t _id;
  };
  // ========================================================
  // Cellgraph Class methods
  // ========================================================
  ~CellGraph2() override;

  /// \param region domain region
  /// \param [optional] maxLevel maximum node level
  explicit CellGraph2(BBox2d region, size_t maxLevel = 13);

  /// \param region domain region
  /// \param f node callback, recursion stops if f returns false
  /// \param [optional] maxLevel maximum node level
  CellGraph2(const BBox2d &region, const std::function<bool(const Node &)> &f,
             size_t maxLevel = 13);

  /// Splits node into 4 children of node.level + 1
  /// \param nodeId node id
  /// \return list of new ids
  std::vector<size_t> refine(size_t nodeId);

  /// Checks if nodes can be coarsed
  /// \param nodesId list of nodes id
  /// \return true if nodes are children of the same parent
  bool canCoarse(std::vector<size_t> nodesId) const;

  /// Searches for node siblings. The case where only diagonal sibling exists is
  /// not handled.
  /// \param nodeId node id
  /// \return list of nodes, including nodeId
  std::vector<size_t> nodeSiblings(size_t nodeId) const;

  /// \param nodeId node to be coarsed
  /// \return new node id, 0 otherwise
  size_t coarse(size_t nodeId);

  /// \param nodeId node id
  /// \return list of ids from neighbor nodes
  std::vector<Definitions::Neighbor> neighbors(size_t nodeId) const;

  /// \param nodeId node id
  /// \param side neighbors side
  /// \return list of ids from neighbor nodes of side <parameter side>
  std::vector<size_t> neighbors(size_t nodeId, Definitions::Side side) const;

  /// \param nodeId node id (vector index)
  /// \param f node callback, recursion stops if f returns false
  void refine(size_t nodeId, std::function<bool(const Node &)> f);

  /// Refine all cells present in the structure by one level. It does not take
  /// in consideration differences between face levels, refining them
  /// accordingly.
  void refineAll();

  /// \return number of active nodes
  size_t nodeCount() const;

  /// \return number of active edges
  size_t edgeCount() const;

  /// \param nodeId
  /// \return node struct reference
  const Node &node(size_t nodeId) const;

  /// \param edgeId
  /// \return edge struct reference
  const Edge &edge(size_t edgeId) const;

  /// \param vertexId  vertex index
  /// \return vertex struct reference
  const Vertex &vertex(size_t vertexId) const;

  /// Vertices, nodes and edges must respect their order!
  /// \param vertices
  /// \param nodes
  /// \param edges
  void setGraph(const std::vector<Vertex> &vertices,
                const std::vector<Node> &nodes, const std::vector<Edge> &edges);

  /// Dump into stdout nodes, edges and vertices information
  void print() const;

  // ==============================================================
  // STRUCTURE INTERFACE
  size_t cellCount() const override;

  size_t faceCount() const override;

  size_t vertexCount() const override;

  double maxScalarFieldValue(size_t fieldId) const override;

  double minScalarFieldValue(size_t fieldId) const override;

  Point2d cellCenterPosition(size_t cellId) const override;

  Point2d faceCenterPosition(size_t faceId) const override;

  Point2d vertexPosition(size_t vertexId) const override;

  BBox2d cellRegion(size_t cellId) const override;

  BBox2d domainRegion() const override;

  void iterateCells(const std::function<void(size_t)> &f) const override;

  void iterateFaces(const std::function<void(size_t)> &f) const override;

  void iterateVertices(const std::function<void(size_t)> &f) const override;

#ifdef _USE_OPENMP
  void iterateCells_par(const std::function<void(size_t)> &f) const override;

  void iterateFaces_par(const std::function<void(size_t)> &f) const override;

  void iterateVertices_par(const std::function<void(size_t)> &f) const override;
#endif // _USE_OPENMP

  int cellId(Point2d target) const override;

  Definitions::Orientation faceOrientation(size_t faceId) const override;

  std::vector<Definitions::Neighbor>
  cellNeighbors(size_t cellId) const override;

  std::vector<Definitions::Neighbor> faceCells(size_t faceId) const override;

  std::vector<size_t> cellFaces(size_t cellId) const override;

  virtual std::vector<size_t> cellVertices(size_t cellId) const override;

  void
  iterateNeighborFaces(Point2d center, size_t star,
                       Definitions::Orientation orientationFilter,
                       const std::function<void(size_t)> &f) const override;

  void
  iterateNeighborCells(size_t cellId, size_t star,
                       const std::function<void(size_t)> &f) const override;

  void iterateCellRing(size_t cellId,
                       const std::function<void(size_t)> &f) const override;

  double
  sampleFaceField(size_t fieldId, Point2d target,
                  Definitions::Orientation orientationFilter,
                  const std::function<bool(size_t)> &faceIsValid = nullptr,
                  const std::function<bool(size_t)> &cellIsValid = nullptr,
                  Definitions::BoundaryVelocity boundaryVelocity =
                      Definitions::BoundaryVelocity::NONE) const override;

  double sampleCellField(
      size_t fieldId, Point2d target,
      const std::function<bool(size_t)> &isValid = nullptr) const override;

  double sampleField(
      size_t fieldId, Point2d target,
      const std::function<bool(size_t)> &isValid = nullptr) const override;

  double sampleVertexField(size_t fieldId, Point2d target) const override;

  void transferFaceFieldToVertexField(
      size_t faceFieldId, size_t vertexFieldId,
      Definitions::Orientation orientationFilter,
      const std::function<bool(size_t)> &isValid) override;

  void transferCellFieldToFaceField(
      size_t cellFieldId, size_t faceFieldId,
      const std::function<bool(size_t)> &cellIsValid,
      const std::function<bool(size_t)> &faceIsValid) override;

protected:
  size_t allCellDataCount() const override;

  size_t allFaceDataCount() const override;

  size_t allVertexDataCount() const override;

  /// checks avaliable node ids queue first or allocates a new one. Its added to
  /// each vertex as well. \param r vertices (in ccw order starting from lower
  /// left corner) \param level node level \return the index of a new node
  size_t newNode(const size_t vertices[], size_t level = 0);

  /// checks avaliable edge ids queue first or allocates a new one. bSide is
  /// computed automatically. \param aSide a relative side to center \param
  /// vertexA node a index \param vertexB node b index \param center face center
  /// position \return the index of a new edge
  size_t newEdge(Definitions::Side aSide, size_t vertexA, size_t vertexB,
                 Point2d center);

  /// checks avaliable vertex ids queue first or allocates a new one
  /// \param position vertex position
  /// \return the index of the new vertex
  size_t newVertex(const Point2d &position);

  /// send node id to available node ids queue
  /// \param nodeId
  void destroyNode(size_t nodeId);

  /// send node id to available edge ids queue
  /// \param edgeId
  void destroyEdge(size_t edgeId);

  /// send node id to available vertex ids queue
  /// \param vertexId
  void destroyVertex(size_t vertexId);

  /// - Allocates a new node to each child of the given node.
  /// - New vertices are allocated as well
  /// - The address code of each child is computed
  /// \param nodeId node index to be refined
  /// \return children indices in the following order:
  /// 0 1
  /// 2 3
  std::vector<size_t> createChildren(size_t nodeId);

  /// Given 4 nodes ids (generated by refining a parent node), this method
  /// creates edges between each of them
  /// \param children children ids expected in the following order:
  /// 0 1
  /// 2 3
  void interconnectChildren(const std::vector<size_t> &children);

  /// Finds the middle vertex of each face of a given node. The middle face
  /// vertex is the vertex that is located in the center position of the face.
  /// (this method is used in the refine process, when we have to assign new/old
  /// vertices to node children. \param node node id \param vertices indices
  /// from sides L R B T, -1 if vertex does not exists
  void middleNodeFaceVertices(size_t node, int *vertices);

  /// Creates edges between children to a neighbor with region size <= to child
  /// region size \param children list of children in order 0 1 2 3 \param
  /// neighbor neighbor object
  void connectChildrenToSmallNeighbor(const std::vector<size_t> &children,
                                      Definitions::Neighbor neighbor);
  /// Creates edges between children to a neighbor with region size > to child
  /// region size \param children list of children in order 0 1 2 3 \param
  /// neighbor neighbor object
  void connectChildrenToBigNeighbor(const std::vector<size_t> &children,
                                    Definitions::Neighbor neighbor);

  std::vector<Node> _nodes;
  std::vector<Edge> _edges;
  std::vector<Vertex> _vertices;
  std::queue<size_t> _availableNodeIds;
  std::queue<size_t> _availableEdgeIds;
  std::queue<size_t> _availableVertexIds;
  size_t _maxLevel;
  BBox2d _domainRegion;
  std::shared_ptr<Interpolant<Point2d>> _interpolator;
};

class CellGraph3 : public StructureInterface<3> {
public:
  class Edge;
  // ========================================================
  class Node {
  public:
    /// \param vertices node vertices in ccw order from lower left corner
    /// \param r node region
    /// \param id node id (accessed as cellId)
    /// \param level node level
    explicit Node(const size_t vertices[], const BBox3d &r, size_t id = 0,
                  size_t level = 0);

    /// Sets node data and makes it valid
    /// \param vertices node vertices sorted respecting x < y < z, where z has
    /// highest precedence \param r node region \param level node tree level
    void set(const size_t vertices[], const BBox3d &r, size_t level);

    /// \return node's region
    BBox3d region() const;

    /// \return node's level
    size_t level() const;

    /// Vertex index of corner i
    /// \param i vertex index (z order x < y < z)
    /// \return vertex id
    size_t vertexAt(size_t i) const;

    /// Reference of edges set
    /// \return edges connected to the node
    const std::set<size_t> &edges() const;

    ///
    /// \return vertices in the node corner
    std::vector<size_t> vertices() const;

    /// Copy of edges set
    /// \return edges connected to the node
    std::set<size_t> edgesCopy() const;

    /// Connect edge to node
    /// \param e edge index
    void addEdge(size_t e);

    /// Disconnect edge from node, just erases its index from the set of edges
    /// \param e edge index
    void removeEdge(size_t e);

    /// \return true if node id is valid
    bool isValid() const;

    /// \return node id
    size_t id() const;

    /// make invalid and clear all edges
    void destroy();

    /// The addres code of the first corner vertex
    /// \return node's address code
    Point3u addressCode() const;

    /// Computes the address code of a vertex corner of the current node
    /// \param vertexNodeId corner index (0 - 7)
    /// \param maxLevel max cell graph level
    /// \return the address code of that vertexNodeId
    Point3u addressCode(size_t vertexNodeId, size_t maxLevel) const;

    /// Sets node address code
    /// \param addressCode new address code
    void setAddressCode(Point3u addressCode);

  private:
    std::set<size_t> _edges; //!< list of edge indices connected to this node
    BBox3d _region;          //!< node's region
    size_t _level;           //!< region's level of hierarchy
    size_t _id;              //!< node id (accessed as cellId)
    // TODO: improve address code description
    Point3u _addressCode; //!< relative unique id to maxLevel, used to compute
                          //!< siblings
    bool _isValid;        //!< true if node id is valid
    size_t _cornerVertices[8]{}; //!< list of vertex indices connected to this
                                 //!< node,
    //!< vertices are always stored in ccw order from lower left corner
  };
  // ========================================================
  class Edge {
  public:
    /// Constructs a valid edge based on 2 vertices and their side relative to
    /// each other. From aSide, bSide is computed automatically. \param aSide
    /// position of vertexA relative to edge center \param vertexA node a index
    /// \param vertexB node b index
    /// \param center face center position
    /// \param id edge id (faceId)
    Edge(Definitions::Side aSide, size_t vertexA, size_t vertexB,
         Point3d center, size_t id);

    /// Sets an edge based on 2 vertices and their side relative to each other.
    /// From aSide, bSide is computed automatically. Makes the edge valid.
    /// \param aSide position of vertexA relative to edge center
    /// \param vertexA node a index
    /// \param vertexB node b index
    /// \param center face center position
    void set(Definitions::Side aSide, size_t vertexA, size_t vertexB,
             Point3d center);

    /// Edge's orientation relative to the cartesian vectors
    /// \return edge's orientation
    Definitions::Orientation orientation() const;

    /// Given one vertex id of this edge, returns the side and id of the other
    /// vertex \param nodeId node id \return node neighbor [side and id]
    Definitions::Neighbor nodeNeighbor(size_t nodeId) const;

    /// A valid node is a node that is currently present in the graph
    /// \return true if node id is valid
    bool isValid() const;

    /// \return edge's id
    size_t id() const;

    /// makes edge invalid
    void destroy();

    std::vector<size_t> nodes() const;

    Point3d center() const;

  private:
    Definitions::Orientation _orientation =
        Definitions::Orientation::CUSTOM; //!< edge orientation (horizontal or
                                          //!< vertical)
    Definitions::Side
        _nodeSides[2]{};  //!< vertex position related to edge's center
    size_t _nodeIds[2]{}; //!< nodes connected by the edge
    Point3d _center;      //!< used as face center position
    bool _isValid = true;
    size_t _id;
  };
  // ========================================================
  class Vertex {
  public:
    /// Constructs a valid vertex from its position
    /// \param position vertex's position
    /// \param id vertex's index
    Vertex(Point3d position, size_t id);

    /// Sets a vertex from its position. Makes it valid.
    /// \param position vertex's position
    /// \param id vertex's index
    void set(Point3d position, size_t id);

    /// Clears its cell counter and makes it invalid.
    void destroy();

    /// Decrements by one its cell counter
    void removeCell();

    /// Increments by one its cell counter
    void addCell();

    /// Number of cells currently connected to it
    /// \return cells counter value
    size_t cellCount() const;

    /// \return vertex position
    Point3d position() const;

    /// A valid vertex is a vertex currently present in the graph
    /// \return true if valid
    bool isValid() const;

    /// \return vertex index
    size_t id() const;

  private:
    size_t _cellCount = 0;
    Point3d _position; //!< vertex location
    bool _isValid = true;
    size_t _id;
  };
  // ========================================================
  ~CellGraph3() override;

  /// Setup a single node graph covering the input region
  /// \param region domain region
  /// \param [optional] maxLevel maximum node level
  explicit CellGraph3(BBox3d region, size_t maxLevel = 13);

  /// Setup a graph based on a refinement criteria
  /// \param region domain region
  /// \param f node callback, recursion stops if f returns false
  /// \param [optional] maxLevel maximum node level
  CellGraph3(const BBox3d &region, const std::function<bool(const Node &)> &f,
             size_t maxLevel = 13);

  /// \param nodeId
  /// \return node struct reference
  const Node &node(size_t nodeId) const;

  /// \param edgeId
  /// \return edge struct reference
  const Edge &edge(size_t edgeId) const;

  /// \param vertexId  vertex index
  /// \return vertex struct reference
  const Vertex &vertex(size_t vertexId) const;

  /// The neighbors from all node edges (does not modify any id)
  /// \param nodeId node id
  /// \return list of ids from neighbor nodes
  // TODO: this method could be protected?..
  std::vector<Definitions::Neighbor> neighbors(size_t nodeId) const;

  /// Splits node into 8 children of node.level + 1
  /// \param nodeId node id
  /// \return list of new ids
  std::vector<size_t> refine(size_t nodeId);

  /// \param nodeId node id (vector index)
  /// \param f node callback, recursion stops if f returns false
  void refine(size_t nodeId, std::function<bool(const Node &)> f);

  /// Searches for node siblings. The case where only diagonal sibling exists is
  /// not handled.
  /// \param nodeId node id
  /// \return list of nodes, including nodeId
  std::vector<size_t> nodeSiblings(size_t nodeId) const;

  /// Merges, if possible, 8 sibling nodes into their respective parent node
  /// \param nodeId node to be coarsed
  /// \return new node id, 0 otherwise
  size_t coarse(size_t nodeId);

  /// Vertices, nodes and edges must respect their order!
  /// \param vertices
  /// \param nodes
  /// \param edges
  void setGraph(const std::vector<Vertex> &vertices,
                const std::vector<Node> &nodes, const std::vector<Edge> &edges);

  /// Dump into stdout nodes, edges and vertices information
  void print() const;

  // ==============================================================
  // STRUCTURE INTERFACE
  size_t cellCount() const override;

  size_t faceCount() const override;

  size_t vertexCount() const override;

  Point3d cellCenterPosition(size_t cellId) const override;

  Point3d faceCenterPosition(size_t faceId) const override;

  Point3d vertexPosition(size_t vertexId) const override;

  BBox3d cellRegion(size_t cellId) const override;

  BBox3d domainRegion() const override;

  void iterateCells(const std::function<void(size_t)> &f) const override;

  void iterateFaces(const std::function<void(size_t)> &f) const override;

  void iterateVertices(const std::function<void(size_t)> &f) const override;

#ifdef _USE_OPENMP
  void iterateCells_par(const std::function<void(size_t)> &f) const override;

  void iterateFaces_par(const std::function<void(size_t)> &f) const override;

  void iterateVertices_par(const std::function<void(size_t)> &f) const override;
#endif // _USE_OPENMP

  int cellId(Point3d target) const override;

  Definitions::Orientation faceOrientation(size_t faceId) const override;

  std::vector<Definitions::Neighbor>
  cellNeighbors(size_t cellId) const override;

  std::vector<Definitions::Neighbor> faceCells(size_t faceId) const override;

  std::vector<size_t> cellFaces(size_t cellId) const override;

  virtual std::vector<size_t> cellVertices(size_t cellId) const override;

  void
  iterateNeighborFaces(Point3d center, size_t star,
                       Definitions::Orientation orientationFilter,
                       const std::function<void(size_t)> &f) const override;

  void
  iterateNeighborCells(size_t cellId, size_t star,
                       const std::function<void(size_t)> &f) const override;

  void iterateCellRing(size_t cellId,
                       const std::function<void(size_t)> &f) const override;

  double
  sampleFaceField(size_t fieldId, Point3d target,
                  Definitions::Orientation orientationFilter,
                  const std::function<bool(size_t)> &faceIsValid = nullptr,
                  const std::function<bool(size_t)> &cellIsValid = nullptr,
                  Definitions::BoundaryVelocity boundaryVelocity =
                      Definitions::BoundaryVelocity::NONE) const override;

  double sampleCellField(
      size_t fieldId, Point3d target,
      const std::function<bool(size_t)> &isValid = nullptr) const override;

  double sampleField(
      size_t fieldId, Point3d target,
      const std::function<bool(size_t)> &isValid = nullptr) const override;

  double sampleVertexField(size_t fieldId, Point3d target) const override;

  void transferFaceFieldToVertexField(
      size_t faceFieldId, size_t vertexFieldId,
      Definitions::Orientation orientationFilter,
      const std::function<bool(size_t)> &isValid) override;

  void transferCellFieldToFaceField(
      size_t cellFieldId, size_t faceFieldId,
      const std::function<bool(size_t)> &cellIsValid,
      const std::function<bool(size_t)> &faceIsValid) override;

private:
  size_t allCellDataCount() const override;

  size_t allFaceDataCount() const override;

  size_t allVertexDataCount() const override;

  /// checks avaliable node ids queue first or allocates a new one. Its added to
  /// each vertex as well.
  /// \param r vertices (in ccw order starting from lower left corner)
  /// \param level node level
  /// \return the index of a new node
  size_t newNode(const size_t vertices[], size_t level = 0);

  /// checks avaliable edge ids queue first or allocates a new one. bSide is
  /// computed automatically.
  /// \param aSide a relative side to center
  /// \param vertexA node a index
  /// \param vertexB node b index
  /// \param center face center position
  /// \return the index of a new edge
  size_t newEdge(Definitions::Side aSide, size_t vertexA, size_t vertexB,
                 Point3d center);

  /// checks avaliable vertex ids queue first or allocates a new one
  /// \param position vertex position
  /// \return the index of the new vertex
  size_t newVertex(const Point3d &position);

  /// send node id to available node ids queue
  /// \param nodeId
  void destroyNode(size_t nodeId);

  /// send node id to available edge ids queue
  /// \param edgeId
  void destroyEdge(size_t edgeId);

  /// send node id to available vertex ids queue
  /// \param vertexId
  void destroyVertex(size_t vertexId);

  /// - Allocates a new node to each child of the given node.
  /// - New vertices are allocated as well
  /// - The address code of each child is computed
  /// \param nodeId node index to be refined
  /// \return children indices in the z order
  std::vector<size_t> createChildren(size_t nodeId);

  /// Given 8 nodes ids (generated by refining a parent node), this method
  /// creates edges between each of them
  /// \param children children ids expected in z order
  void interconnectSiblings(const std::vector<size_t> &children);

  /// Lists the indices of all vertices that belong to the children of a given
  /// node. If any vertex does not exists, this method creates a new one.
  /// \param node node id
  /// \return vertex indices in sorted order by x < y < z
  std::vector<size_t> splitVertices(size_t node);

  /// Creates edges between children to a neighbor with region size <= to child
  /// region size
  /// \param children list of children in order x < y < z
  /// \param neighbor neighbor object
  void connectChildrenToSmallNeighbor(const std::vector<size_t> &children,
                                      Definitions::Neighbor neighbor);
  /// Creates edges between children to a neighbor with region size > to child
  /// region size
  /// \param children list of children in order x < y < z
  /// \param neighbor neighbor object
  void connectChildrenToBigNeighbor(const std::vector<size_t> &children,
                                    Definitions::Neighbor neighbor);

  /// Connects a node to its neighbors
  /// \param nodeId parent node index
  /// \param _neighbors list of neighbor cells
  void connectParentToNeighbors(
      size_t nodeId, const std::vector<Definitions::Neighbor> &_neighbors);

  /// Checks if nodes can be coarsed
  /// \param nodesId list of nodes id
  /// \return true if nodes are children of the same parent
  bool canCoarse(std::vector<size_t> nodesId) const;

  std::vector<Node> _nodes;
  std::vector<Edge> _edges;
  std::vector<Vertex> _vertices;
  std::queue<size_t> _availableNodeIds;
  std::queue<size_t> _availableEdgeIds;
  std::queue<size_t> _availableVertexIds;
  size_t _maxLevel;
  BBox3d _domainRegion;
  std::shared_ptr<Interpolant<Point3d>> _interpolator;
};

} // namespace furoo

#endif // FUROO_CELL_GRAPH_H
