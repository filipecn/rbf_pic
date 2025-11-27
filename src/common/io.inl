template <int D>
bool IO::saveFields(std::string filename, StructureInterface<D> *structure) {
  FILE *fp = fopen(filename.c_str(), "w+");
  if (!fp) {
    std::cerr << "IO::saveFields failed to open file " << filename << std::endl;
    return false;
  }
  size_t fieldCount = structure->fieldCount();
  /// # of fields FIELD_1_TYPE ... FIELD_N_TYPE
  fprintf(fp, "%lu", fieldCount);
  for (size_t i = 0; i < fieldCount; i++)
    fprintf(fp, " %s %lu %s ", structure->fieldDataTypeName(i).c_str(),
            static_cast<size_t>(structure->fieldLocation(i)),
            structure->fieldName(i).c_str());
  fprintf(fp, "\n");
  // FIELDS IDS
  std::vector<size_t> cellFields, faceFields, vertexFields;
  for (size_t i = 0; i < fieldCount; i++)
    switch (structure->fieldLocation(i)) {
    case Definitions::MeshLocation::VERTEX_CENTER:
      vertexFields.emplace_back(i);
      break;
    case Definitions::MeshLocation::CELL_CENTER:
      cellFields.emplace_back(i);
      break;
    case Definitions::MeshLocation::FACE_CENTER:
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      faceFields.emplace_back(i);
      break;
    default:
      std::cerr << "IO::saveFields field " << i
                << " location will not be stored";
    }
  std::vector<size_t> &fieldIds = cellFields;
  std::function<void(size_t)> elementCallback = [&](size_t id) {
    fprintf(fp, "%lu", id);
    std::ostringstream s;
    for (size_t i : fieldIds) {
      s << " " << i << " ";
      if (structure->fieldDataTypeName(i) == "double")
        s << (*structure->template field<double>(i))[id];
      else if (structure->fieldDataTypeName(i) == "int")
        s << (*structure->template field<int>(i))[id];
      else if (structure->fieldDataTypeName(i) == "unsigned_char")
        s << static_cast<size_t>(
            (*structure->template field<unsigned char>(i))[id]);
      else if (structure->fieldDataTypeName(i) == "Definitions::Material")
        s << static_cast<size_t>(
            (*structure->template field<Definitions::Material>(i))[id]);
      else if (structure->fieldDataTypeName(i) == "Definitions::Boundary")
        s << static_cast<size_t>(
            (*structure->template field<Definitions::Boundary>(i))[id]);
      else if (structure->fieldDataTypeName(i) == "Point<double,D>") {
        auto p = (*structure->template field<Point<double, D>>(i))[id];
        for (int d = 0; d < D; d++)
          s << p[d] << " ";
      }
    }
    fprintf(fp, "%s\n", s.str().c_str());
  };
  // store fields
  /// CELL_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  fprintf(fp, "%lu %lu\n", structure->cellCount(), fieldIds.size());
  structure->iterateCells(elementCallback);
  /// FACE_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  fieldIds = faceFields;
  fprintf(fp, "%lu %lu\n", structure->faceCount(), fieldIds.size());
  structure->iterateFaces(elementCallback);
  /// VERTEX_INDEX FIELD_ID FIELD_VALUE ... FIELD_ID FIELD_VALUE
  fieldIds = vertexFields;
  fprintf(fp, "%lu %lu\n", structure->vertexCount(), fieldIds.size());
  structure->iterateVertices(elementCallback);
  fclose(fp);
  return true;
}

template <int D>
bool IO::loadFields(std::string filename, StructureInterface<D> *structure) {
  std::ifstream fp(filename, std::ifstream::in);
  if (!fp.good()) {
    std::cerr << "IO::loadFields failed to open file " << filename << std::endl;
    return false;
  }
  size_t fieldCount = 0;
  /// # of fields FIELD_1_TYPE ... FIELD_N_TYPE
  fp >> fieldCount;
  structure->removeFields();
  for (size_t i = 0; i < fieldCount; i++) {
    std::string fieldType, fieldName;
    size_t locationId = 0;
    fp >> fieldType >> locationId >> fieldName;
    if (fieldType == "double")
      structure->template addField<double>(
          Definitions::meshLocationFromId(locationId), 0., fieldName);
    else if (fieldType == "unsigned_char")
      structure->template addField<unsigned char>(
          Definitions::meshLocationFromId(locationId), 0, fieldName);
    else if (fieldType == "int")
      structure->template addField<int>(
          Definitions::meshLocationFromId(locationId), 0, fieldName);
    else if (fieldType == "Definitions::Material")
      structure->template addField<Definitions::Material>(
          Definitions::meshLocationFromId(locationId),
          Definitions::Material::AIR, fieldName);
    else if (fieldType == "Definitions::Boundary")
      structure->template addField<Definitions::Boundary>(
          Definitions::meshLocationFromId(locationId),
          Definitions::Boundary::NONE, fieldName);
    else if (fieldType == "Point<double,D>")
      structure->template addField<Point<double, D>>(
          Definitions::meshLocationFromId(locationId), Point<double, D>(),
          fieldName);
  }
  auto readElementFields = [&](size_t nFields) {
    size_t elementId = 0;
    fp >> elementId;
    for (size_t f = 0; f < nFields; f++) {
      size_t fieldId = 0;
      fp >> fieldId;
      std::string fieldType = structure->fieldDataTypeName(fieldId);
      if (fieldType == "double") {
        auto &field = *structure->template field<double>(fieldId);
        fp >> field[elementId];
      } else if (fieldType == "unsigned_char") {
        auto &field = *structure->template field<unsigned char>(fieldId);
        size_t value = 0;
        fp >> value;
        field[elementId] = value;
      } else if (fieldType == "int") {
        auto &field = *structure->template field<int>(fieldId);
        fp >> field[elementId];
      } else if (fieldType == "Definitions::Material") {
        auto &field =
            *structure->template field<Definitions::Material>(fieldId);
        size_t materialId = 0;
        fp >> materialId;
        field[elementId] = Definitions::materialFromId(materialId);
      } else if (fieldType == "Definitions::Boundary") {
        auto &field =
            *structure->template field<Definitions::Boundary>(fieldId);
        size_t boundaryId = 0;
        fp >> boundaryId;
        field[elementId] = Definitions::boundaryFromId(boundaryId);
      } else if (fieldType == "Point<double,D>") {
        auto &field = *structure->template field<Point<double, D>>(fieldId);
        Point<double, D> p;
        for (int d = 0; d < D; d++) {
          std::string a;
          fp >> a;
          if (a == "inf")
            p[d] = INFINITY;
        }
        field[elementId] = p;
      }
    }
  };
  size_t cellCount = 0, nFields = 0;
  fp >> cellCount >> nFields;
  if (cellCount != structure->cellCount()) {
    std::cerr << "IO::loadFields incompatible cell graph (cells " << cellCount
              << " x " << structure->cellCount() << ")" << filename
              << std::endl;
    return false;
  }
  for (size_t i = 0; i < structure->cellCount(); i++)
    readElementFields(nFields);
  size_t faceCount = 0;
  fp >> faceCount >> nFields;
  if (faceCount != structure->faceCount()) {
    std::cerr << "IO::loadFields incompatible cell graph (faces " << faceCount
              << " x " << structure->faceCount() << ")" << filename
              << std::endl;
    return false;
  }
  for (size_t i = 0; i < structure->faceCount(); i++)
    readElementFields(nFields);
  size_t vertexCount = 0;
  fp >> vertexCount >> nFields;
  if (vertexCount != structure->vertexCount()) {
    std::cerr << "IO::loadFields incompatible cell graph (vertices "
              << vertexCount << " x " << structure->vertexCount() << ")"
              << filename << std::endl;
    return false;
  }
  for (size_t i = 0; i < structure->vertexCount(); i++)
    readElementFields(nFields);
  fp.close();
  return true;
}

template <>
inline bool IO::saveCellGraph(std::string filename,
                              StructureInterface<2> *structure) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  FILE *fp = fopen(filename.c_str(), "w+");
  if (!fp) {
    std::cerr << "IO::saveCellGraph failed to open " << filename << std::endl;
    return false;
  }
  // store vertices
  fprintf(fp, "%lu\n", cellGraph->vertexCount());
  cellGraph->iterateVertices([&](size_t vertexId) {
    /// vertexId positionX_i positionY_j
    auto &vertex = cellGraph->vertex(vertexId);
    fprintf(fp, "%lu %lf %lf\n", vertex.id(), vertex.position().x(),
            vertex.position().y());
  });
  // store cells
  fprintf(fp, "%lu\n", cellGraph->cellCount());
  cellGraph->iterateCells([&](size_t i) {
    auto &node = cellGraph->node(i);
    const auto &edges = node.edges();
    /// cellId v0 v1 v2 v3 cellLevel addressCodeI addressCodeJ # of edges e0 e1
    /// ...
    fprintf(fp, "%lu %lu %lu %lu %lu %lu %d %d %lu", node.id(), node.vertex(0),
            node.vertex(1), node.vertex(2), node.vertex(3), node.level(),
            node.addressCode().x(), node.addressCode().y(), edges.size());
    for (size_t e : edges)
      fprintf(fp, " %lu", e);
    fprintf(fp, "\n");
  });
  // store edges
  fprintf(fp, "%lu\n", cellGraph->faceCount());
  cellGraph->iterateFaces([&](size_t e) {
    auto &edge = cellGraph->edge(e);
    /// edgeId aSide a b centerX centerY
    fprintf(fp, "%lu %lu %lu %lu %lf %lf\n", edge.id(),
            static_cast<size_t>(edge.nodeNeighbor(edge.nodes()[1]).side),
            edge.nodes()[0], edge.nodes()[1], edge.center().x(),
            edge.center().y());
  });
  fclose(fp);
  return true;
}

template <>
inline bool IO::saveCellGraph(std::string filename,
                              StructureInterface<3> *structure) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  std::ofstream fp(filename, std::ofstream::out);
  if (!fp.good()) {
    std::cerr << "IO::saveCellGraph failed to open " << filename << std::endl;
    return false;
  }
  // store vertices
  fp << cellGraph->vertexCount() << std::endl;
  cellGraph->iterateVertices([&](size_t vertexId) {
    /// vertexId positionX_i positionY_j
    auto &vertex = cellGraph->vertex(vertexId);
    fp << vertex.id() << " " << vertex.position().x() << " "
       << vertex.position().y() << " " << vertex.position().z() << std::endl;
  });
  // store cells
  fp << cellGraph->cellCount() << std::endl;
  cellGraph->iterateCells([&](size_t i) {
    auto &node = cellGraph->node(i);
    const auto &edges = node.edges();
    /// cellId v0 v1 v2 v3 cellLevel addressCodeI addressCodeJ # of edges e0 e1
    /// ...
    fp << node.id() << " ";
    for (size_t i = 0; i < 8; i++)
      fp << node.vertexAt(i) << " ";
    fp << node.level() << " ";
    fp << node.addressCode().x() << " ";
    fp << node.addressCode().y() << " ";
    fp << node.addressCode().z() << " ";
    fp << edges.size() << " ";
    for (size_t e : edges)
      fp << e << " ";
    fp << std::endl;
  });
  // store edges
  fp << cellGraph->faceCount() << std::endl;
  cellGraph->iterateFaces([&](size_t e) {
    auto &edge = cellGraph->edge(e);
    /// edgeId aSide a b centerX centerY
    fp << edge.id() << " ";
    fp << static_cast<size_t>(edge.nodeNeighbor(edge.nodes()[1]).side) << " ";
    fp << edge.nodes()[0] << " ";
    fp << edge.nodes()[1] << " ";
    fp << edge.center().x() << " ";
    fp << edge.center().y() << " ";
    fp << edge.center().z() << std::endl;
  });
  fp.close();
  return true;
}

template <>
inline bool IO::loadCellGraph(std::string filename,
                              StructureInterface<2> *structure) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  std::ifstream fp(filename, std::ifstream::in);
  if (!fp.good())
    return false;
  std::vector<CellGraph2::Vertex> vertices;
  std::vector<CellGraph2::Node> nodes;
  std::vector<CellGraph2::Edge> edges;
  size_t cellCount = 0, vertexCount = 0, edgeCount = 0;
  fp >> vertexCount;
  for (size_t vertexId = 0; vertexId < vertexCount; vertexId++) {
    double x, y;
    size_t id;
    fp >> id >> x >> y;
    if (vertices.size() <= id)
      for (size_t k = vertices.size(); k <= id; k++) {
        vertices.emplace_back(Point2d(x, y), k);
        vertices[k].destroy();
      }
    vertices[id].set(Point2d(x, y), id);
    // comentario muito importante//
  }
  fp >> cellCount;
  { // add outside domain node first
    size_t v[4] = {0, 1, 2, 3};
    nodes.emplace_back(v, BBox2d::squareBox(1.), 0u, 0u);
  }
  for (size_t i = 0; i < cellCount; i++) {
    size_t v[4] = {0, 0, 0, 0};
    size_t id, level, n;
    int ax, ay;
    fp >> id >> v[0] >> v[1] >> v[2] >> v[3] >> level >> ax >> ay >> n;
    BBox2d region(vertices[v[0]].position(), vertices[v[2]].position());
    if (id >= nodes.size())
      for (size_t k = nodes.size(); k <= id; k++) {
        nodes.emplace_back(v, region, k, level);
        nodes[k].destroy();
      }
    nodes[id].set(v, region, level);
    nodes[id].setAddressCode(Point2u(ax, ay));
    for (size_t e = 0; e < n; e++) {
      size_t edgeId;
      fp >> edgeId;
      nodes[id].addEdge(edgeId);
    }
  }
  fp >> edgeCount;
  for (size_t i = 0; i < edgeCount; i++) {
    double x(0), y(0);
    size_t sideId(0), nodeA(0), nodeB(0), id(0);
    fp >> id >> sideId >> nodeA >> nodeB >> x >> y;
    if (edges.size() <= id)
      for (size_t k = edges.size(); k <= id; k++) {
        edges.emplace_back(Definitions::sideFromId(sideId), nodeA, nodeB,
                           Point2d(x, y), k);
        edges[k].destroy();
      }
    edges[id].set(Definitions::sideFromId(sideId), nodeA, nodeB, Point2d(x, y));
  }
  cellGraph->setGraph(vertices, nodes, edges);
  fp.close();
  return true;
}

template <>
inline bool IO::loadCellGraph(std::string filename,
                              StructureInterface<3> *structure) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  std::ifstream fp(filename, std::ifstream::in);
  if (!fp.good())
    return false;
  std::vector<CellGraph3::Vertex> vertices;
  std::vector<CellGraph3::Node> nodes;
  std::vector<CellGraph3::Edge> edges;
  size_t cellCount = 0, vertexCount = 0, edgeCount = 0;
  fp >> vertexCount;
  for (size_t vertexId = 0; vertexId < vertexCount; vertexId++) {
    double x, y, z;
    size_t id;
    fp >> id >> x >> y >> z;
    if (vertices.size() <= id)
      for (size_t k = vertices.size(); k <= id; k++) {
        vertices.emplace_back(Point3d(x, y, z), k);
        vertices[k].destroy();
      }
    vertices[id].set(Point3d(x, y, z), id);
    // comentario muito importante//
  }
  fp >> cellCount;
  { // add outside domain node first
    size_t v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    nodes.emplace_back(v, BBox3d::squareBox(1.), 0u, 0u);
  }
  for (size_t i = 0; i < cellCount; i++) {
    size_t v[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    size_t id, level, n;
    int ax, ay, az;
    fp >> id >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >>
        level >> ax >> ay >> az >> n;
    BBox3d region(vertices[v[0]].position(), vertices[v[7]].position());
    if (id >= nodes.size())
      for (size_t k = nodes.size(); k <= id; k++) {
        nodes.emplace_back(v, region, k, level);
        nodes[k].destroy();
      }
    nodes[id].set(v, region, level);
    nodes[id].setAddressCode(Point3u(ax, ay, az));
    for (size_t e = 0; e < n; e++) {
      size_t edgeId;
      fp >> edgeId;
      nodes[id].addEdge(edgeId);
    }
  }
  fp >> edgeCount;
  for (size_t i = 0; i < edgeCount; i++) {
    double x(0), y(0), z(0);
    size_t sideId(0), nodeA(0), nodeB(0), id(0);
    fp >> id >> sideId >> nodeA >> nodeB >> x >> y >> z;
    if (edges.size() <= id)
      for (size_t k = edges.size(); k <= id; k++) {
        edges.emplace_back(Definitions::sideFromId(sideId), nodeA, nodeB,
                           Point3d(x, y, z), k);
        edges[k].destroy();
      }
    edges[id].set(Definitions::sideFromId(sideId), nodeA, nodeB,
                  Point3d(x, y, z));
  }
  cellGraph->setGraph(vertices, nodes, edges);
  fp.close();
  return true;
}

template <typename T, int D>
bool IO::saveToPV(std::string filename, ParticleSystem<T, D> *ps) {
  if (filename.find(".pv") == std::string::npos)
    filename += ".pv";
  std::ofstream fp(filename, std::ofstream::out);
  if (!fp.good()) {
    std::cerr << "IO::saveToPV failed to open file " << filename << std::endl;
    return false;
  }
  std::vector<std::pair<size_t, Point<double, D>>> points;
  fp << ps->size() << std::endl;
  ps->iterateParticles(
      [&](size_t i, Point<double, D> p) { points.emplace_back(i, p); });
  std::sort(points.begin(), points.end());
  // scalar properties count
  int scalar_properties_count = ps->scalarPropertiesCount();
  for (auto p : points) {
    fp << p.first;
    for (size_t d = 0; d < D; d++)
      fp << " " << p.second[d];
    for (size_t d = 0; d < scalar_properties_count; d++)
      fp << " " << ps->getScalarProperty(d, p.first);
    fp << std::endl;
  }
  fp.close();
  return true;
}

template <typename T, int D>
bool IO::loadFromPV(std::string filename, ParticleSystem<T, D> *ps) {
  if (filename.find(".pv") == std::string::npos)
    filename += ".pv";
  std::ifstream fp(filename, std::ifstream::in);
  if (!fp.good()) {
    std::cerr << "IO::loadFromPV failed to open file " << filename << std::endl;
    return false;
  }
  std::vector<double> p(D);
  int nParticles(0);
  fp >> nParticles;
  ps->clear();
  int scalar_properties_count = ps->scalarPropertiesCount();
  std::vector<T> v(scalar_properties_count);
  for (int size = 0; size < nParticles; size++) {
    size_t pid = 0;
    fp >> pid;
    size_t id = ps->addParticle(Point<double, D>());
    for (size_t d = 0; d < D; d++)
      fp >> p[d];
    for (size_t d = 0; d < scalar_properties_count; d++)
      fp >> v[d];
    ps->setPosition(id, Point<double, D>(p));
    for (size_t d = 0; d < scalar_properties_count; d++)
      ps->setScalarProperty(d, id, v[d]);
  }
  ps->size();
  fp.close();
  return true;
}

template <int D> bool IO::loadOBJ(std::string filename, Mesh<D> *mesh) {
  if (!mesh)
    return false;
  // tinyobj::attrib_t attrib;
  // std::vector<tinyobj::shape_t> shapes;
  // std::vector<tinyobj::material_t> materials;
  // std::string err;
  // bool r =
  //     tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str());
  // if (!r)
  //   return false;
  // if (attrib.vertices.size() % D) {
  //   std::cerr << "IO::loadOBJ OBJ with wrong dimensions\n";
  //   return false;
  // }
  // size_t nvertices = attrib.vertices.size() / D;
  // for (int i = 0; i < nvertices; i++) {
  //   Point<double, D> vertex;
  //   for (int d = 0; d < D; d++)
  //     vertex[d] = attrib.vertices[i * D + d];
  //   mesh->addVertex(vertex);
  // }
  // if (shapes[0].mesh.indices.size() % D) {
  //   std::cerr << "IO::loadOBJ OBJ with wrong dimensions\n";
  //   return false;
  // }
  // size_t nelements = shapes[0].mesh.indices.size() / D;
  // for (int i = 0; i < nvertices; i++) {
  //   Vector<int, D> element;
  //   for (int d = 0; d < D; d++)
  //     element[d] = shapes[0].mesh.indices[i * D + d].vertex_index;
  //   mesh->addElement(element);
  // }
  return true;
}
