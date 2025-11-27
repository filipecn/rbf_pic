#ifndef FUROO_COMMON_IO_H
#define FUROO_COMMON_IO_H

#include <fstream>
#include <geometry/mesh.h>
#include <iostream>
#include <istream>
#include <physics/particle_system.h>
#include <physics/simulation_domain.h>
#include <physics/simulation_quad_tree.h>
#include <sstream>
#include <structures/cell_graph.h>
// #include <tiny_obj_loader.h>

namespace furoo {

class IO {
public:
  /// Loads an obj file into a mesh structure. Assumes triangles for 3D and
  /// line segments for 2D objects
  /// NOTE: only the first shape of the file is loaded
  /// \tparam D dimensions
  /// \param filename path/to/file.obj
  /// \param mesh mesh structure pointer
  /// \return true if success
  template <int D> static bool loadOBJ(std::string filename, Mesh<D> *mesh);

  /// store fields into a file following:
  /// # of fields FIELD_1_TYPE FIELD_1_LOC FIELD_1_NAME ... FIELD_N_TYPE
  /// FIELD_N_LOC FIELD_N_NAME # of cells # of cell fields CELL_INDEX FIELD_ID
  /// FIELD_VALUE ... FIELD_ID FIELD_VALUE # of faces # of face fields
  /// FACE_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  /// # of vertices # of vertex fields
  /// VERTEX_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  /// \tparam D structure dimension
  /// \param filename path/to/file
  /// \param structure structure raw pointer
  /// \return true if success
  template <int D>
  static bool saveFields(std::string filename,
                         StructureInterface<D> *structure);
  /// loads fields from a file following:
  /// # of fields FIELD_1_TYPE FIELD_1_LOC FIELD_1_NAME ... FIELD_N_TYPE
  /// FIELD_N_LOC FIELD_N_NAME # of cells # of cell fields CELL_INDEX FIELD_ID
  /// FIELD_VALUE ... FIELD_ID FIELD_VALUE # of faces # of face fields
  /// FACE_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  /// # of vertices # of vertex fields
  /// VERTEX_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  /// \tparam D structure dimension
  /// \param filename path/to/file
  /// \param structure structure raw pointer
  /// \return true if success
  template <int D>
  static bool loadFields(std::string filename,
                         StructureInterface<D> *structure);

  /// stores cell geometry / topology into a file following:
  /// # of vertices
  /// vertexId positionX_i positionY_j
  /// # of cells
  /// cellId v0 v1 v2 v3 cellLevel addressCodeI addressCodeJ # of edges e0 e1
  /// ... # of edges edgeId aSide a b centerX centerY
  /// \param filename path/to/file
  /// \param cellGraph raw pointer to cell graph data
  /// \return true if success
  template <int D>
  static bool saveCellGraph(std::string filename,
                            StructureInterface<D> *cellGraph);

  /// loads cell geometry / topology following the file template:
  /// # of vertices
  /// vertexId positionX_i positionY_j
  /// # of cells
  /// cellId v0 v1 v2 v3 cellLevel addressCodeI addressCodeJ # of edges e0 e1
  /// ... # of edges edgeId aSide a b centerX centerY
  /// \param filename path/to/file
  /// \param cellGraph raw pointer to cell graph data \return true if success
  template <int D>
  static bool loadCellGraph(std::string filename,
                            StructureInterface<D> *cellGraph);

  static bool saveToKim(const char *, const ParticleSystem2 *);

  /// The pv file type stores a list of particle positions and velocities as
  /// # of particles
  /// id px py vx vy
  /// The properties used are 0 to vx and 1 to vy
  /// \param filename path/to/file
  /// \param ps raw pointer to particle system data
  /// \return true if success
  template <typename T, int D>
  static bool saveToPV(std::string filename, ParticleSystem<T, D> *ps);

  /// The pv file type stores a list of particle positions and velocities as
  /// # of particles
  /// id px py vx vy
  /// The properties used are 0 to vx and 1 to vy
  /// \param filename path/to/file
  /// \param ps raw pointer to particle system data
  /// \return true if success
  template <typename T, int D>
  static bool loadFromPV(std::string filename, ParticleSystem<T, D> *ps);

  static bool saveToCF(const char *, const SimulationDomain2 *, unsigned);

  static bool loadFromCF(const char *, SimulationDomain2 *, unsigned);

  static bool saveToQT(const char *, const SimQuadTree *);

  static bool loadFromQT(const char *, SimQuadTree *);
};

#include "io.inl"

} // namespace furoo

#endif // FUROO_COMMON_IO_H
