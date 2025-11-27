template <int D>
void SimulationSteps::computeFaceMidPoint(StructureInterface<D> *structure,
                                          size_t faceMidPointFieldId,
                                          size_t faceMaterialFieldId) {
  auto &faceMidPoint =
      *structure->template field<Point<double, D>>(faceMidPointFieldId);
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);

  structure->iterateFaces([&](size_t faceId) {
    try {
      auto cells = structure->faceCells(faceId);
      Point<double, D> p;
      double count = 0;
      for (auto cell : cells) {
        p += (structure->cellCenterPosition(cell.id) - Point<double, D>());
        count += 1;
      }
      faceMidPoint[faceId] = p / count;
    } catch (std::string e) {
      std::cerr << e << std::endl;
    }
  });
}

template <int D>
void SimulationSteps::transferVelocitiesFromParticlesToMidPointFaces(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    const std::vector<size_t>& particleVelocityPropertyIds,
    const std::vector<size_t>& faceVelocityFieldIds, size_t faceMaterialFieldId,
    size_t cellMaterialFieldId, size_t faceMidPointFieldId) {
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &faceMidPoint =
      *structure->template field<Point<double, D>>(faceMidPointFieldId);

  std::vector<typename StructureInterface<D>::template Field<double> *>
      faceVelocityFields(D, nullptr);
  for (int d = 0; d < D; d++)
    faceVelocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);
  particles->size();
  structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] == Definitions::Material::AIR)
      return;
    // TODO: Change particles sampling from knn to radius search
    double radius = 0;
    auto faceCenter = faceMidPoint[faceId];
    auto cells = structure->faceCells(faceId);
    bool hasFluid = false;
    for (auto cell : cells)
      if (cellMaterialField[cell.id] == Definitions::Material::FLUID)
        hasFluid = true;
    if (!hasFluid)
      return;
    // evaluate cell's size
    // Cells from same level: radius from face center to farthest corner
    // Cells from differente level: radius from face center to farthest corner
    // of the smallest corner
    size_t smallestCell;
    if (cells.size() == 1)
      smallestCell = cells[0].id;
    else if (structure->cellRegion(cells[0].id).size() <
             structure->cellRegion(cells[1].id).size())
      smallestCell = cells[0].id;
    else
      smallestCell = cells[1].id;
    auto cellVertices = structure->cellVertices(smallestCell);
    for (auto vertex : cellVertices)
      radius = std::max(
          radius, (faceCenter - structure->vertexPosition(vertex)).length());

    auto velocity = particles->sampleVectorFromProperties(
        particleVelocityPropertyIds, radius, faceMidPoint[faceId]);

    if (faceMaterialField[faceId] == Definitions::Material::SOLID) {
      auto faceCells = structure->faceCells(faceId);
      if (faceCells.size() == 2 &&
          (cellMaterialField[faceCells[0].id] == Definitions::Material::FLUID ||
           cellMaterialField[faceCells[1].id] ==
               Definitions::Material::FLUID)) {
        auto orientation = Definitions::sideOrientation(faceCells[0].side);
        switch (orientation) {
        case Definitions::Orientation::HORIZONTAL:
          velocity[0] = 0.;
          break;
        case Definitions::Orientation::VERTICAL:
          velocity[1] = 0.;
          break;
        case Definitions::Orientation::DEPTH:
          velocity[2] = 0.;
          break;
        }
      } else
        return;
    }
    for (int d = 0; d < D; d++)
      (*faceVelocityFields[d])[faceId] = velocity[d];
  });
}

template <int D>
void SimulationSteps::applyExternalForceToMidPointFaces(
    StructureInterface<D> *structure,
    const std::vector<size_t> &faceVelocityFieldIds,
    Vector<double, D> acceleration, double dt, size_t faceMaterialFieldId,
    size_t cellMaterialFieldId) {
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      faceVelocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    faceVelocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);

  structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] == Definitions::Material::AIR)
      return;
    Vector<double, D> velocity = acceleration * dt;
    if (faceMaterialField[faceId] == Definitions::Material::SOLID) {
      auto faceCells = structure->faceCells(faceId);
      if (faceCells.size() == 2 &&
          (cellMaterialField[faceCells[0].id] == Definitions::Material::FLUID ||
           cellMaterialField[faceCells[1].id] ==
               Definitions::Material::FLUID)) {
        auto orientation = Definitions::sideOrientation(faceCells[0].side);
        switch (orientation) {
        case Definitions::Orientation::HORIZONTAL:
          velocity[0] = 0.;
          break;
        case Definitions::Orientation::VERTICAL:
          velocity[1] = 0.;
          break;
        case Definitions::Orientation::DEPTH:
          velocity[2] = 0.;
          break;
        }
      } else
        return;
    }
    for (int d = 0; d < D; d++)
      (*faceVelocityFields[d])[faceId] += velocity[d];
  });
}

template <>
inline void SimulationSteps::computeDivergenceFieldFromMidPointFaces(
    StructureInterface<2> *structure, size_t divergenceCellFieldId,
    size_t faceMidPointFieldId, const std::vector<size_t> &velocityFaceFieldIds,
    const std::function<bool(size_t)> &useCell, DifferentialRBF<2> *rbf) {
  int D = 2;
  auto cellGraph = dynamic_cast<CellGraph2 *>(structure);
  std::vector<typename StructureInterface<2>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &faceMidPoint =
      *structure->template field<Point<double, 2>>(faceMidPointFieldId);
  divergenceField.setAll(0.);
  // TODO: change to flux method
  structure->iterateCells([&](size_t cellId) {
    if (!useCell(cellId))
      return;
    auto cellCenter = structure->cellCenterPosition(cellId);
    auto neighborFaces = structure->cellFaces(cellId);
    std::vector<Point<double, 2>> points;
    std::vector<double> values[2];
    bool sameLevel = true;
    std::vector<double> flows;
    for (auto face : neighborFaces) {
      auto cells = structure->faceCells(face);
      if (cells.size() == 2 && cells[0].id >= 0 && cells[1].id >= 0 &&
          cellGraph->node(cells[0].id).level() !=
              cellGraph->node(cells[1].id).level())
        sameLevel = false;
      points.emplace_back(faceMidPoint[face]);
      Vector<double, 2> faceVelocity;
      for (int d = 0; d < D; d++) {
        values[d].emplace_back((*velocityFields[d])[face]);
        faceVelocity[d] = (*velocityFields[d])[face];
      }
      Vector<double, 2> direction = faceMidPoint[face] - cellCenter;
      double flow = faceVelocity.dot(2. * direction);
      flows.emplace_back(flow);
    }
    divergenceField[cellId] = 0.;
    if (!sameLevel)
      divergenceField[cellId] =
          0.25 * rbf->laplacianAt(cellCenter, points, flows);
    else
      for (int d = 0; d < D; d++)
        divergenceField[cellId] +=
            rbf->gradientAt(cellCenter, d, points, values[d]);
    // }
  });
}

template <>
inline void SimulationSteps::computeDivergenceFieldFromMidPointFaces(
    StructureInterface<3> *structure, size_t divergenceCellFieldId,
    size_t faceMidPointFieldId, const std::vector<size_t> &velocityFaceFieldIds,
    const std::function<bool(size_t)> &useCell, DifferentialRBF<3> *rbf) {
  int D = 3;
  auto cellGraph = dynamic_cast<CellGraph3 *>(structure);
  std::vector<typename StructureInterface<3>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &faceMidPoint =
      *structure->template field<Point<double, 3>>(faceMidPointFieldId);
  divergenceField.setAll(0.);
  // TODO: change to flux method
  structure->iterateCells([&](size_t cellId) {
    if (!useCell(cellId))
      return;
    auto cellCenter = structure->cellCenterPosition(cellId);
    auto neighborFaces = structure->cellFaces(cellId);
    std::vector<Point<double, 3>> points;
    std::vector<double> values[3];
    bool sameLevel = true;
    std::vector<double> flows;
    for (auto face : neighborFaces) {
      auto cells = structure->faceCells(face);
      if (cells.size() == 2 && cells[0].id >= 0 && cells[1].id >= 0 &&
          cellGraph->node(cells[0].id).level() !=
              cellGraph->node(cells[1].id).level())
        sameLevel = false;
      points.emplace_back(faceMidPoint[face]);
      Vector<double, 3> faceVelocity;
      for (int d = 0; d < D; d++) {
        values[d].emplace_back((*velocityFields[d])[face]);
        faceVelocity[d] = (*velocityFields[d])[face];
      }
      Vector<double, 3> direction = faceMidPoint[face] - cellCenter;
      double flow = faceVelocity.dot(2. * direction);
      flows.emplace_back(flow);
    }
    divergenceField[cellId] = 0.;
    if (!sameLevel)
      divergenceField[cellId] =
          0.25 * rbf->laplacianAt(cellCenter, points, flows);
    else
      for (int d = 0; d < D; d++)
        divergenceField[cellId] +=
            rbf->gradientAt(cellCenter, d, points, values[d]);
    // }
  });
}

template <int D>
void SimulationSteps::transferMidPointFacesFlowToCellVelocities(
    StructureInterface<D> *structure, size_t cellPressureFieldId,
    size_t cellMaterialFieldId, size_t faceMaterialFieldId,
    size_t faceMidPointFieldId, const std::vector<size_t> &cellVelocityIds,
    const std::vector<size_t> &faceVelocityFieldIds, DifferentialRBF<D> *rbf,
    double dt) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellPressureField =
      *structure->template field<double>(cellPressureFieldId);
  auto &faceMidPoint =
      *structure->template field<Point<double, D>>(faceMidPointFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      cellVelocityFields(D, nullptr), faceVelocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++) {
    cellVelocityFields[d] =
        structure->template field<double>(cellVelocityIds[d]);
    faceVelocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);
  }

  structure->iterateCells([&](size_t cellId) {
    for (int d = 0; d < D; d++)
      (*cellVelocityFields[d])[cellId] = 0;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      // setup stencil
      std::vector<size_t> faceIds;
      std::vector<Point<double, D>> points;
      auto cellFaces = structure->cellFaces(cellId);
      for (auto f : cellFaces)
        faceIds.emplace_back(f);
      std::vector<double> flows;
      auto cellCenter = structure->cellCenterPosition(cellId);
      for (size_t i = 0; i < faceIds.size(); i++) {
        size_t faceId = faceIds[i];
        // Compute pressure gradient to flow correction
        auto faceCells = structure->faceCells(faceId);
        if (faceCells.size() < 2)
          continue;
        if (faceCells[0].id < 0 || faceCells[1].id < 0)
          continue;
        if (cellMaterialField[faceCells[0].id] ==
                Definitions::Material::SOLID ||
            cellMaterialField[faceCells[1].id] ==
                Definitions::Material::SOLID) {
          flows.emplace_back(0.);
          points.emplace_back(faceMidPoint[faceId]);
          continue;
        }
        double pressureGradient = cellPressureField[faceCells[1].id] -
                                  cellPressureField[faceCells[0].id];
        if (faceCells[0].id != cellId)
          pressureGradient = -pressureGradient;

        // Project the current face velocity to obtain flow value
        auto facePosition = faceMidPoint[faceId];
        Vector<double, D> direction = facePosition - cellCenter;
        Vector<double, D> velocity;
        for (int d = 0; d < D; d++)
          velocity[d] = (*faceVelocityFields[d])[faceId];
        double currentFlow =
            velocity.dot(direction * 2) - dt * pressureGradient;
        flows.emplace_back(currentFlow);
        points.emplace_back(faceMidPoint[faceId]);
      }
      for (int d = 0; d < D; d++)
        (*cellVelocityFields[d])[cellId] =
            0.5 * rbf->gradientAt(structure->cellCenterPosition(cellId), d,
                                  points, flows);
    }
  });
}
