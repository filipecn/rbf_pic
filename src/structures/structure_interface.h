#ifndef FUROO_STRUCTURE_INTERFACE_H
#define FUROO_STRUCTURE_INTERFACE_H

#include <common/definitions.h>
#include <geometry/bbox.h>
#include <geometry/vector.h>
#include <memory>
#include <type_traits>

namespace furoo {

/// Encapsulates a spatial data structure that breaks the domain into discrete
/// elements and allows fast iteration between cells, faces and vertices.
/// Also stores and provides access to some types of data (scalar, byte, ...)
/// arranged by different fields localized in cell centers, face centers and
/// vertices.
/// \tparam D data dimension
template <int D> class StructureInterface {
public:
  class FieldInterface {
  public:
    virtual void clear() = 0;
    /// increases field size by size positions (same as append)
    /// \param size number of new positions
    virtual void increaseSize(size_t size) = 0;

    /// resets a field position to the default value
    /// \param i position index, if < 0 reset all positions values
    virtual void reset(int i) = 0;

    /// \return string containing the type of data the field holds
    virtual std::string fieldDataTypeName() const = 0;

    /// sets field string identifier
    /// \param name identifier
    void setName(const std::string &name);

    /// \return field string identifier
    std::string name() const;

  protected:
    std::string _name;
  };
  template <typename T> class Field : public FieldInterface {
  public:
    typedef T FieldType;
    Field(size_t size, T value);

    Field(size_t size, T value, std::string fieldName);

    void setAll(T value);

    /// \param i element index
    /// \return copy of value at index i
    T operator[](size_t i) const;

    /// \param i element index
    /// \return reference of value at index i
    T &operator[](size_t i);

    /// remove field data, sizes becomes 0
    void clear() override;

    /// increases field size by <size>
    /// \param size number of new elements
    void increaseSize(size_t size) override;

    /// Sets field's element value to default value (if element index < 0, sets
    /// all elements)
    /// \param i element index
    void reset(int i) override;

    std::string fieldDataTypeName() const override;

  private:
    T _defaultValue;
    std::vector<T> _data;
  };

  StructureInterface();

  virtual ~StructureInterface();

  /// Adds a cell field to the structure
  /// \tparam T field data type
  /// \param value field default value
  /// \param fieldName field string identifier
  /// \return field id
  template <typename T>
  size_t addCellField(T value, std::string fieldName = "");

  /// Adds a face field to the structure
  /// \tparam T field data type
  /// \param fieldName field string identifier
  /// \param value field default value
  /// \return field id
  template <typename T>
  size_t addFaceField(T value, std::string fieldName = "",
                      Definitions::MeshLocation location =
                          Definitions::MeshLocation::FACE_CENTER);

  /// Adds a vertex field to the structure
  /// \tparam T field data type
  /// \param value field default value
  /// \param fieldName field string identifier
  /// \return field id
  template <typename T>
  size_t addVertexField(T value, std::string fieldName = "");

  /// Adds a field to the structure
  /// \tparam T field data type
  /// \param location field location
  /// \param value field default value
  /// \param fieldName field string identifier
  /// \return field id
  template <typename T>
  size_t addField(Definitions::MeshLocation location, T value,
                  std::string fieldName = "");

  /// Remove all fields
  void removeFields();

  /// Changes the type of a field
  /// \tparam T new type
  /// \param fieldId field index
  /// \param defaultValue new field default value
  template <typename T> void changeFieldType(size_t fieldId, T defaultValue);

  /// \tparam T field type
  /// \param fieldId field index
  /// \return raw pointer to field object
  template <typename T> Field<T> *field(size_t fieldId);

  /// \tparam T field type
  /// \param fieldId field index
  /// \return raw pointer to field object
  template <typename T> const Field<T> *field(size_t fieldId) const;

  /// \param fieldId field index
  /// \return field data type name string
  std::string fieldDataTypeName(size_t fieldId) const;

  /// \param fieldId field index
  /// \return field string identifier
  std::string fieldName(size_t fieldId) const;

  /// \return number of stored fields
  size_t fieldCount() const;

  /// \param fieldId field index
  /// \return field location (cell centers, face centers, etc...)
  Definitions::MeshLocation fieldLocation(size_t fieldId) const;

  /// Maximum value of field stored in a valid point
  /// \param fieldId field id
  /// \return the maximum value
  virtual double maxScalarFieldValue(size_t fieldId) const;

  /// Minimum value of field stored in a valid point
  /// \param fieldId field id
  /// \return the minimum value
  virtual double minScalarFieldValue(size_t fieldId) const;

  /// \return number of active cells
  virtual size_t cellCount() const = 0;

  /// \return number of active faces
  virtual size_t faceCount() const = 0;

  /// \return number of active vertices
  virtual size_t vertexCount() const = 0;

  /// \param cellId cell id
  /// \return bounding box of region covered by cell id
  virtual BBox<double, D> cellRegion(size_t cellId) const = 0;

  /// \return domain region
  virtual BBox<double, D> domainRegion() const = 0;

  /// Iterate over valid cells
  /// \param f callback to each cell (cellId)
  virtual void iterateCells(const std::function<void(size_t)> &f) const = 0;

  /// Iterate over valid faces
  /// \param f callback to each face (faceId)
  virtual void iterateFaces(const std::function<void(size_t)> &f) const = 0;

  /// Iterate over valid vertices
  /// \param f callback to each vertex (vertexId)
  virtual void iterateVertices(const std::function<void(size_t)> &f) const = 0;

#ifdef _USE_OPENMP
  /// Iterate over valid cells in parallel using OpenMP
  /// \param f callback (no race conditions) to each cell (cellId)
  virtual void iterateCells_par(const std::function<void(size_t)> &f) const = 0;

  /// Iterate over valid faces in parallel using OpenMP
  /// \param f callback (no race conditions) to each face (faceId)
  virtual void iterateFaces_par(const std::function<void(size_t)> &f) const = 0;

  /// Iterate over valid vertices in parallel using OpenMP
  /// \param f callback (no race conditions) to each vertex (vertexId)
  virtual void
  iterateVertices_par(const std::function<void(size_t)> &f) const = 0;
#endif // _USE_OPENMP

  /// Faces orientations are orthogonal to structure's edge orientations.
  /// Example: On a MAC grid, velocities aligned to x axis will lie on VERTICAL
  /// orientation
  /// \param faceId
  /// \return face's orientation
  virtual Definitions::Orientation faceOrientation(size_t faceId) const = 0;

  /// Filter neighbors by byte value
  /// \param cellId
  /// \return list of neighbors with byte value
  virtual std::vector<Definitions::Neighbor>
  cellNeighbors(size_t cellId) const = 0;

  /// \param p target point
  /// \return id of cell containing target, -1 if p is outside domain
  virtual int cellId(Point<double, D> p) const = 0;

  /// cell center position
  /// \param cellId cell index
  /// \return center position of cell with index cellId
  virtual Point<double, D> cellCenterPosition(size_t cellId) const = 0;

  /// face center position
  /// \param faceId face index
  /// \return center position of face with index faceId
  virtual Point<double, D> faceCenterPosition(size_t faceId) const = 0;

  /// vertex position
  /// \param vertexId vertex index
  /// \return position of vertex with index vertexId
  virtual Point<double, D> vertexPosition(size_t vertexId) const = 0;

  /// Cells that share face faceId.
  /// Note: Don't assume the result will always have size 2, since boundary
  /// faces have only one incident cell
  /// \param faceId face index
  /// \return list of neighbors cells
  virtual std::vector<Definitions::Neighbor> faceCells(size_t faceId) const = 0;

  /// \param cellId
  /// \return list of faces conected to cell
  virtual std::vector<size_t> cellFaces(size_t cellId) const = 0;

  ///
  /// \param cellId
  /// \return std::vector<size_t>
  virtual std::vector<size_t> cellVertices(size_t cellId) const = 0;

  /// Iterates over faces that lie in the topological neighborhood of the cell
  /// containing center. Examples of stars:
  /// star = 0: faces of center cell
  /// star = 1: faces of star 0 + faces of first neighboring cells
  /// star = 2: faces o star 1 + faces of neighbors of first neighboring cells
  /// \param center query point
  /// \param star size of topological stencil
  /// \param orientationFilter FACE ORIENTATION (orthogonal to structure's edge
  /// orientation) \param f callback for each face
  virtual void
  iterateNeighborFaces(Point<double, D> center, size_t star,
                       Definitions::Orientation orientationFilter,
                       const std::function<void(size_t)> &f) const = 0;

  /// Iterate cell neighbors of a cell following the topological neighborhood of
  /// size star.
  /// \param cellId center cell
  /// star = 0: just cellId
  /// star = 1: neighbors of 0
  /// star = 2: neighbors of 1
  /// \param star star size
  /// \param f callback for each cell
  virtual void
  iterateNeighborCells(size_t cellId, size_t star,
                       const std::function<void(size_t)> &f) const = 0;

  /// Iterate the topological ring of cells centered on cellId
  /// \param cellId center cell
  /// \param f callback for each valid cell
  virtual void iterateCellRing(size_t cellId,
                               const std::function<void(size_t)> &f) const = 0;

  /// Interpolates/Approximates the value of field fieldId at target
  /// \param fieldId field index
  /// \param target query point
  /// \param orientationFilter FACE ORIENTATION (orthogonal to structure's edge
  /// orientation)
  /// \param faceIsValid callback function for each candidate face, returns true
  /// if face enters on the computation \param cellIsValid callback function for
  /// each cell, returns true if cell can be used as source (used on slip
  /// conditions) \return interpolated/approximated value at target
  // TODO: will be removed later
  virtual double
  sampleFaceField(size_t fieldId, Point<double, D> target,
                  Definitions::Orientation orientationFilter,
                  const std::function<bool(size_t)> &faceIsValid = nullptr,
                  const std::function<bool(size_t)> &cellIsValid = nullptr,
                  Definitions::BoundaryVelocity boundaryVelocity =
                      Definitions::BoundaryVelocity::NONE) const = 0;

  /// Interpolates/Approximates the value of field fieldId at target
  /// \param fieldId field index
  /// \param target query point
  /// \param isValid callback function for each candidate face, returns true if
  /// face enters on the computation
  /// \return interpolated/approximated value at target
  // TODO: will be removed later
  virtual double sampleCellField(
      size_t fieldId, Point<double, D> target,
      const std::function<bool(size_t)> &isValid = nullptr) const = 0;

  /// Interpolates/Approximates the value of field fieldId at target
  /// \param fieldId  field index
  /// \param target query point
  /// \return interpolated/approximated value at target
  // TODO: will be removed later
  virtual double sampleVertexField(size_t fieldId,
                                   Point<double, D> target) const = 0;

  /// Interpolates/Approximates the value of field fieldId at target
  /// \param fieldId field index
  /// \param target query point
  /// \param isValid callback function for each candidate element, returns true
  /// if element is used on interpolation.
  /// For example: If fieldId is an index of a cell centered field, then isValid
  /// will receive cell ids.
  virtual double
  sampleField(size_t fieldId, Point<double, D> target,
              const std::function<bool(size_t)> &isValid = nullptr) const = 0;

  /// Interpolates values on a face field to vertices positions and stores in a
  /// vertex field
  /// \param faceFieldId face field index
  /// \param vertexFieldId vertex field index
  /// \param orientationFilter FACE ORIENTATION (orthogonal to structure's edge
  /// orientation)
  /// \param isValid callback function for each candidate face, returns true if
  /// face enters on the computation
  // TODO: to be removed
  virtual void transferFaceFieldToVertexField(
      size_t faceFieldId, size_t vertexFieldId,
      Definitions::Orientation orientationFilter,
      const std::function<bool(size_t)> &isValid) = 0;

  /// Interpolates values on a cell field to face center positions and stores in
  /// a face field
  /// \param cellFieldId cell field index
  /// \param faceFieldId face field index
  /// \param cellIsValid callback function for each candidate cell, returns true
  /// if face enters on the computation
  /// \param faceIsValid callback function for each candidate face, returns true
  /// if face enters on the computation
  // TODO: to be removed
  virtual void transferCellFieldToFaceField(
      size_t cellFieldId, size_t faceFieldId,
      const std::function<bool(size_t)> &cellIsValid,
      const std::function<bool(size_t)> &faceIsValid) = 0;

protected:
  /// Some structures may vary the number of cells over time meaning that
  /// cells may have indices greater than cellCount(). It is important
  /// when adding a new cell field to the structure, since we must store
  /// enough data for these indices as well
  /// \return cellCount() + invalid cell count
  virtual size_t allCellDataCount() const = 0;

  /// Some structures may vary the number of faces over time meaning that
  /// faces may have indices greater than faceCount(). It is important
  /// when adding a new vertex field to the structure, since we must store
  /// enough data for these indices as well
  /// \return faceCount() + invalid face count
  virtual size_t allFaceDataCount() const = 0;

  /// Some structures may vary the number of vertices over time meaning that
  /// vertices may have indices greater than vertexCount(). It is important
  /// when adding a new vertex field to the structure, since we must store
  /// enough data for these indices as well
  /// \return vertexCount() + invalid vertex count
  virtual size_t allVertexDataCount() const = 0;

  std::vector<Definitions::MeshLocation>
      _fieldsLocations; //!< location of each field
  std::vector<std::shared_ptr<FieldInterface>> _fields; //!< all fields data
  std::vector<size_t> _cellFieldIds; //!< ids of fields that are stored in cells
  std::vector<size_t> _faceFieldIds; //!< ids of fields that are stored in faces
  std::vector<size_t>
      _vertexFieldIds; //!< ids of fields that are stored in vertices
  std::vector<Definitions::Boundary> _boundaryTypes;
  std::vector<Vector<double, D>> _boundaryValues;
};

#include "structure_interface.inl"
#include "structure_interface_field.inl"

} // namespace furoo

#endif // FUROO_STRUCTURE_INTERFACE_H
