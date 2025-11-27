template <int D>
void SimulationSteps::applyFunctionToParticles(
    ParticleSystem<double, D> *particles, size_t testFunctionPropertyId,
    const std::function<double(Point<double, D>)> &f) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    particles->setScalarProperty(testFunctionPropertyId, id, f(p));
  });
}

template <int D>
void SimulationSteps::applyFunctionToCells(
    StructureInterface<D> *structure, size_t testFunctionFieldId,
    size_t cellParticleCenterFieldId,
    const std::function<double(Point<double, D>)> &f) {
  auto &field = *structure->template field<double>(testFunctionFieldId);
  structure->iterateCells([&](size_t cellId) {
    field[cellId] = f(structure->cellCenterPosition(cellId));
    std::cerr << structure->cellCenterPosition(cellId) << "\t" << field[cellId]
              << std::endl;
  });
}

template <int D>
void SimulationSteps::computeError(StructureInterface<D> *structure,
                                   size_t resultId, size_t solutionId) {
  auto &resultField = *structure->template field<double>(resultId);
  auto &solutionField = *structure->template field<double>(solutionId);
  structure->iterateCells([&](size_t cellId) {
    resultField[cellId] = solutionField[cellId] - resultField[cellId];
  });
}

template <int D>
void SimulationSteps::computeGradient(
    int dimension, StructureInterface<D> *cellGraph, size_t functionFieldId,
    size_t gradientFieldId, size_t cellMaterialFieldId,
    size_t cellParticleCenterFieldId, DifferentialRBF<D> *rbf) {
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *cellGraph->template field<Point<double, D>>(cellParticleCenterFieldId);
  auto &pressureField = *cellGraph->template field<double>(functionFieldId);
  auto &gradientField = *cellGraph->template field<double>(gradientFieldId);
  std::set<size_t> airSurfaceCells;
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto cellCenter = cellParticleCenterField[cellId];
    auto neighbors = cellGraph->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0 ||
          cellMaterialField[cell.id] == Definitions::Material::SOLID)
        continue;
      points.emplace_back(cellParticleCenterField[cell.id]);
      pressureValues.emplace_back(pressureField[cell.id]);
      if (cellMaterialField[cell.id] == Definitions::Material::AIR)
        airSurfaceCells.insert(cell.id);
    }
    points.emplace_back(cellCenter);
    pressureValues.emplace_back(pressureField[cellId]);
    gradientField[cellId] =
        rbf->gradientAt(cellCenter, dimension, points, pressureValues);
  });
  // now compute pressure gradient on close air as well
  for (size_t cellId : airSurfaceCells) {
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto cellCenter = cellParticleCenterField[cellId];
    auto neighbors = cellGraph->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0 ||
          cellMaterialField[cell.id] == Definitions::Material::SOLID)
        continue;
      points.emplace_back(cellParticleCenterField[cell.id]);
      pressureValues.emplace_back(pressureField[cell.id]);
    }
    points.emplace_back(cellCenter);
    pressureValues.emplace_back(pressureField[cellId]);
    gradientField[cellId] =
        rbf->gradientAt(cellCenter, dimension, points, pressureValues);
  }
}

template <int D>
void SimulationSteps::computeLaplacian(StructureInterface<D> *structure,
                                       size_t functionFieldId,
                                       size_t laplacianFieldId,
                                       size_t cellMaterialFieldId,
                                       size_t cellParticleCenterFieldId,
                                       DifferentialRBF<D> *rbf) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *structure->template field<Point<double, D>>(cellParticleCenterFieldId);
  auto &functionField = *structure->template field<double>(functionFieldId);
  auto &laplacianField = *structure->template field<double>(laplacianFieldId);
  size_t fluidCount = 0;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      fluidCount++;
    }
  });
  if (!fluidCount)
    return;
  // Since the matrix will have rows only to fluid cells, we need a map
  // between cell id and matrix row id
  // structure index -> matrix index
  std::map<size_t, size_t> indexMap;
  std::function<size_t(size_t)> cellIndex = [&](size_t cellId) {
    auto it = indexMap.find(cellId);
    if (it == indexMap.end()) {
      size_t nextIndex = indexMap.size();
      indexMap[cellId] = nextIndex;
      return indexMap.size() - 1;
    }
    return it->second;
  };

  // Lets compute the RHS side of the system, which is = divergence / dt
  Eigen::VectorXd rhs;
  rhs = Eigen::VectorXd::Zero(fluidCount); // +1 is for air cells
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      rhs[cellIndex(cellId)] = functionField[cellId];
  });
  // And now the triplets that form the matrix
  std::vector<Eigen::Triplet<double>> triplets;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    auto neighbors = structure->cellNeighbors(cellId);
    auto weights = computePressureMatrixStencilWeights(
        structure, cellId, cellMaterialFieldId, cellParticleCenterFieldId, rbf,
        false);
    // auto weights = computePressureMatrixStencilWeights(
    //     structure, cellId, cellMaterialFieldId, rbf, false);
    double centerWeight = weights[0];
    for (size_t i = 0; i < neighbors.size(); i++) {
      if (neighbors[i].id >= 0 &&
          cellMaterialField[neighbors[i].id] == Definitions::Material::AIR) {
        continue;
        triplets.emplace_back(Eigen::Triplet<double>(
            cellIndex(cellId), cellIndex(neighbors[i].id), weights[i + 1]));
        triplets.emplace_back(cellIndex(neighbors[i].id),
                              cellIndex(neighbors[i].id), 1);
      }
      // continue;
      else if (neighbors[i].id >= 0 && cellMaterialField[neighbors[i].id] ==
                                           Definitions::Material::FLUID)
        triplets.emplace_back(Eigen::Triplet<double>(
            cellIndex(cellId), cellIndex(neighbors[i].id), weights[i + 1]));
      else
        centerWeight += weights[i + 1];
    }
    triplets.emplace_back(Eigen::Triplet<double>(
        cellIndex(cellId), cellIndex(cellId), centerWeight));
  });

  // and finally solve the system...
  Eigen::SparseMatrix<double> A(indexMap.size(), indexMap.size());
  A.setFromTriplets(triplets.begin(), triplets.end());
  std::cerr << "A.size() " << A.rows() << "x" << A.cols() << '\n';
  std::cerr << "rhs.size() " << rhs.size() << '\n';
  std::cerr << A << std::endl << rhs << std::endl;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver;
  solver.compute(A);
  Eigen::VectorXd result;
  result = solver.solve(rhs);
  structure->iterateCells([&](size_t cellId) {
    laplacianField[cellId] = 0.;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      laplacianField[cellId] = result[cellIndex(cellId)];
  });
}

template <typename T, int D>
using PointFunction = std::function<void(size_t, Point<T, D>)>;

template <typename T, int D>
inline static void iterateAllParticles(ParticleSystem<T, D> *p,
                                       const PointFunction<T, D> &f) {
#ifdef _USE_OPENMP
  p->iterateParticles_par(f);
#else
  p->iterateParticles(f);
#endif // _USE_OPENMP
}

template <int D>
void SimulationSteps::advectParticles(
    ParticleSystem<double, D> *particleSystem,
    const std::vector<size_t> &velocityProperties,
    const StructureInterface<D> *domain,
    const std::vector<std::shared_ptr<SolidInterface<D>>> &solids,
    const TemporalIntegratorInterface<D> *integrator, double dt) {
  std::vector<size_t> removeCanditate;
  iterateAllParticles<double, D>(
      particleSystem,
      [&](size_t particleId, Point<double, D> particlePosition) {
        Vector<double, D> v =
            particleSystem->getVector(velocityProperties, particleId);
        auto np = integrator->integrate(particlePosition, v, dt);
        // Check if point is still inside domain
        auto region = domain->domainRegion();
        if (!region.contains(np)) {
          // particleSystem->remove(i);
          for (int i = 0; i < D; i++) {
            if (np[i] < region.lower()[i] || np[i] > region.upper()[i]) {
              v[i] = 0.;
              np[i] = max(region.lower()[i] + 0.0000001,
                          min(region.upper()[i] - 0.0000001, np[i]));
              particlePosition[i] =
                  max(region.lower()[i] + 0.0000001,
                      min(region.upper()[i] - 0.0000001, particlePosition[i]));
            }
          }
        }
        bool overlap = false;
        for (auto &solid : solids) {
          if (solid->contains(np)) {
            // particleSystem->remove(particleId);
            // return;
            auto intersection = solid->intersectionPoint(particlePosition, np);
            if (!intersection.isValid) {
              particleSystem->remove(particleId);
              break;
            }
            // np = intersection.point;
            for (int d = 0; d < D; d++)
              if (intersection.normal[d] != 0.) {
                v[d] = 0.;
                np[d] =
                    intersection.point[d] + 0.0000001 * intersection.normal[d];
              }
            if (solid->contains(np)) {
              std::cerr << "Intersection normal: " << intersection.normal
                        << std::endl;
              std::cerr << "Particleposition: " << particlePosition
                        << std::endl;
              std::cerr << "New Position: " << np << std::endl;
              std::cerr << "Velocity: " << v << std::endl;
              std::cerr << "Solid box: "
                        << dynamic_cast<SolidBox<D> *>(solid.get())->box()
                        << std::endl;
              particleSystem->remove(particleId);
              std::cerr << "A particle was removed!\n";
#ifdef FUROO_PROFILE
              CodeProfiler::count("#PariclesRemovedOnAdvection");
#endif
              break;
              // THROW(false, "particle still in solid?!");
            }
            // particleSystem->remove(i);
            // return;
            // removeCanditate.emplace_back(particleId);
          }
        }
        // TODO: why do we clamp if particle is inside domain?
        // np = Point2d(clamp(np.x(), 0.0, 0.99999), clamp(np.y(), 0.0,
        // 0.99999));
        particleSystem->setPosition(particleId, np);
        for (int i = 0; i < D; i++)
          particleSystem->setScalarProperty(velocityProperties[i], particleId,
                                            v[i]);
        if (!domain->domainRegion().contains(np))
          particleSystem->remove(particleId);
        // THROW(domain->domainRegion().contains(np),
        //      "SimulationSteps::advectParticles Particles outside domain!");
      });
  std::vector<size_t> toRemove;
  for (auto id : removeCanditate) {
    Point<double, D> np = (*particleSystem)[id];
    particleSystem->iterateClosestParticles(
        np, 1, [&](size_t pId, Point<double, D> pPosition) {
          // TODO: mark particle to be removed
          if ((np - pPosition).length() < 1e-6 && pId < id)
            toRemove.emplace_back(id);
        });
  }
  for (auto particle : toRemove) {
    std::cout << particle << " removed.\n";
    particleSystem->remove(particle);
  }
#ifdef FUROO_PROFILE
  CodeProfiler::count("#PariclesRemovedOnAdvection", toRemove.size());
#endif
}

template <int D>
void SimulationSteps::computeDivergenceOnParticles(
    ParticleSystem<double, D> *particles,
    const std::vector<size_t> &velocityPropertyIds, size_t divergencePropertyId,
    DifferentialRBF<D> *rbf) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    std::vector<Point<double, D>> points;
    std::vector<double> values[D];
    particles->iterateClosestParticles(
        p, 7, [&](size_t nid, Point<double, D> np) {
          if (nid == id)
            return;
          points.emplace_back(np);
          for (int d = 0; d < D; d++)
            values[d].emplace_back(
                particles->getScalarProperty(velocityPropertyIds[d], nid));
        });
    double div = 0;
    for (int d = 0; d < D; d++)
      div += rbf->gradientAt(p, d, points, values[d]);
    particles->setScalarProperty(divergencePropertyId, id, div);
  });
  std::cerr << "DIVERGENCE ON PARTICLES" << std::endl;
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    std::cerr << particles->getScalarProperty(divergencePropertyId, id)
              << std::endl;
  });
}

template <int D>
void SimulationSteps::computeDivergenceField(
    StructureInterface<D> *structure, size_t divergenceCellFieldId,
    const std::vector<size_t> &velocityFaceFieldIds,
    const std::function<bool(size_t)> &useCell, DifferentialRBF<D> *rbf,
    bool adaptiveRBFRadius) {
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  divergenceField.setAll(0.);
  structure->iterateCells([&](size_t cellId) {
    if (!useCell(cellId))
      return;
    auto cellCenter = structure->cellCenterPosition(cellId);
    // std::cerr << "\ncellCenter: " << cellCenter << std::endl;
    auto neighborFaces = structure->cellFaces(cellId);
    std::vector<Point<double, D>> points[D];
    std::vector<double> values[D];
    double div = 0.;
    // std::cerr << "Face positions: " << std::endl;
    for (auto face : neighborFaces) {
      switch (structure->faceOrientation(face)) {
      case Definitions::Orientation::HORIZONTAL:
        // std::cerr << "\tY value: " << structure->faceCenterPosition(face)
        //           << std::endl;
        points[1].emplace_back(structure->faceCenterPosition(face));
        values[1].emplace_back((*velocityFields[1])[face]);
        break;
      case Definitions::Orientation::VERTICAL:
        // std::cerr << "\tX value: " << structure->faceCenterPosition(face)
        //           << std::endl;
        points[0].emplace_back(structure->faceCenterPosition(face));
        values[0].emplace_back((*velocityFields[0])[face]);
        break;
      case Definitions::Orientation::DEPTH:
        // std::cerr << "\tZ value: " << structure->faceCenterPosition(face)
        //           << std::endl;
        points[2].emplace_back(structure->faceCenterPosition(face));
        values[2].emplace_back((*velocityFields[2])[face]);
        break;
      default:
        THROW(false, "SimulationSteps::computeDivergenceField invalid face "
                     "orientation.");
      }
    }
    // X-direction divergent
    for (size_t d = 0; d < D; d++) {
      THROW(points[d].size() == values[d].size(),
            "SimulationSteps::computeDivergence points.size() != "
            "values.size()");
      // std::cerr << "Values: \n";
      // for (auto value : values[d])
      //   std::cerr << value << ' ';
      // std::cerr << std::endl;
      div += rbf->gradientAt(cellCenter, d, points[d], values[d]);
      // std::cerr << "Grad: " << div << std::endl;
    }
    // std::cerr << "Divergent: " << div << std::endl;
    divergenceField[cellId] = div;
  });
}
template <int D>
void SimulationSteps::computeDivergenceField(
    StructureInterface<D> *structure, size_t divergenceCellFieldId,
    const std::vector<size_t> &velocityFaceFieldIds,
    const std::vector<size_t> &velocityCellFieldIds,
    const std::function<bool(size_t)> &useCell, DifferentialRBF<D> *rbf,
    bool adaptiveRBFRadius) {
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFaceFields(D, nullptr), velocityCellFields(D, nullptr);
  for (size_t d = 0; d < D; d++) {
    velocityFaceFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
    velocityCellFields[d] =
        structure->template field<double>(velocityCellFieldIds[d]);
  }
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  divergenceField.setAll(0.);
  structure->iterateCells([&](size_t cellId) {
    if (!useCell(cellId))
      return;
    auto cellCenter = structure->cellCenterPosition(cellId);
    auto neighborFaces = structure->cellFaces(cellId);
    auto neighborCells = structure->cellNeighbors(cellId);
    std::vector<Point<double, D>> points[D];
    std::vector<double> values[D];
    for (auto face : neighborFaces) {
      switch (structure->faceOrientation(face)) {
      case Definitions::Orientation::HORIZONTAL:
        points[1].emplace_back(structure->faceCenterPosition(face));
        values[1].emplace_back((*velocityFaceFields[1])[face]);
        break;
      case Definitions::Orientation::VERTICAL:
        points[0].emplace_back(structure->faceCenterPosition(face));
        values[0].emplace_back((*velocityFaceFields[0])[face]);
        break;
      case Definitions::Orientation::DEPTH:
        points[2].emplace_back(structure->faceCenterPosition(face));
        values[2].emplace_back((*velocityFaceFields[2])[face]);
        break;
      default:
        THROW(false, "SimulationSteps::computeDivergenceField invalid face "
                     "orientation.");
      }
    }
    double div = 0.;
    // X-direction divergent
    for (size_t d = 0; d < D; d++) {
      points[d].emplace_back(cellCenter);
      values[d].emplace_back((*velocityCellFields[d])[cellId]);
      THROW(points[d].size() == values[d].size(),
            "SimulationSteps::computeDivergence points.size() != "
            "values.size()");
      div += rbf->gradientAt(cellCenter, d, points[d], values[d]);
    }
    divergenceField[cellId] = div;
  });
}

template <int D>
void SimulationSteps::computeDivergenceFieldFromParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t divergenceCellFieldId, size_t cellMaterialFieldId,
    size_t cellParticleCenterFieldId,
    const std::vector<size_t> &velocityParticleFieldIds,
    const std::function<bool(size_t)> &useCell, DifferentialRBF<D> *rbf) {
  size_t numberOfWallSamples = (D == 2) ? 3 : 9;
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &materialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &particleCenterField =
      *structure->template field<Point<double, D>>(cellParticleCenterFieldId);
  divergenceField.setAll(0.);
  structure->iterateCells([&](size_t cellId) {
    if (!useCell(cellId))
      return;
    // check if cell touches solid
    auto neighbors = structure->cellNeighbors(cellId);
    auto cellRegion = structure->cellRegion(cellId);
    std::vector<Point<double, D>> points;
    std::vector<size_t> particleIds;
    bool interfaceWithAir = false;
    for (auto neighbor : neighbors) {
      if (neighbor.id >= 0 &&
          materialField[neighbor.id] != Definitions::Material::FLUID) {
        interfaceWithAir = true;
        break;
      }
    }
    Point<double, D> centroid(0.);
    int count = 0;
    if (interfaceWithAir) {
      particles->iterateParticles(cellRegion.expanded(Vector<double, D>(1.95)),
                                  [&](size_t id, Point<double, D> p) {
                                    particleIds.emplace_back(id);
                                    points.emplace_back(p);
                                  });
    } else
      particles->iterateParticles(cellRegion.expanded(Vector<double, D>(1.5)),
                                  [&](size_t id, Point<double, D> p) {
                                    particleIds.emplace_back(id);
                                    points.emplace_back(p);
                                  });

    THROW(particleIds.size(), "SimulationSteps::computeDivergence empty cell");

    auto cellCenter = particleCenterField[cellId];
    std::vector<double> gradient(D, 0.0);
    std::vector<bool> alreadyComputed(D, false);
    auto falloff = [](double d) {
      return d; /*-d * d + 2 * d;*/
    };
    for (auto neighbor : neighbors) {
      if (neighbor.id <= 0 ||
          materialField[neighbor.id] == Definitions::Material::SOLID) {
        std::vector<double> values(points.size(), 0.);
        size_t dimension = 0;
        switch (neighbor.side) {
        case Definitions::Side::LEFT:
          dimension = 0;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);
            values[i] = clamp(falloff((points[i].x() - cellRegion.lower().x()) /
                                      cellRegion.size().x()),
                              0., 1.) *
                        velocity[0];
            // particles->setScalarProperty(velocityParticleFieldIds[0],
            //  particleIds[i], values[i]);
          }
          break;
        case Definitions::Side::RIGHT:
          dimension = 0;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);

            values[i] = clamp(falloff((cellRegion.upper().x() - points[i].x()) /
                                      cellRegion.size().x()),
                              0., 1.) *
                        velocity[0];
            // particles->setScalarProperty(velocityParticleFieldIds[0],
            //  particleIds[i], values[i]);
          }
          break;
        case Definitions::Side::BOTTOM:
          dimension = 1;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);

            values[i] = clamp(falloff((points[i].y() - cellRegion.lower().y()) /
                                      cellRegion.size().y()),
                              0., 1.) *
                        velocity[1];
            // particles->setScalarProperty(velocityParticleFieldIds[1],
            //  particleIds[i], values[i]);
          }
          break;
        case Definitions::Side::TOP:
          dimension = 1;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);

            values[i] = clamp(falloff((cellRegion.upper().y() - points[i].y()) /
                                      cellRegion.size().y()),
                              0., 1.) *
                        velocity[1];
            // particles->setScalarProperty(velocityParticleFieldIds[1],
            //  particleIds[i], values[i]);
          }
          break;
        case Definitions::Side::BACK:
          dimension = 2;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);

            values[i] = clamp(falloff((points[i].z() - cellRegion.lower().z()) /
                                      cellRegion.size().z()),
                              0., 1.) *
                        velocity[2];
            // particles->setScalarProperty(velocityParticleFieldIds[2],
            //  particleIds[i], values[i]);
          }
          break;
        case Definitions::Side::FRONT:
          dimension = 2;
          for (size_t i = 0; i < particleIds.size(); i++) {
            auto velocity =
                particles->getVector(velocityParticleFieldIds, particleIds[i]);
            Vector<double, D> direction = cellCenter - points[i];
            velocity.projectInto(direction);

            values[i] = clamp(falloff((cellRegion.upper().z() - points[i].z()) /
                                      cellRegion.size().z()),
                              0., 1.) *
                        velocity[2];
            // particles->setScalarProperty(velocityParticleFieldIds[2],
            //  particleIds[i], values[i]);
          }
          break;
        }
        gradient[dimension] =
            rbf->gradientAt(cellCenter, dimension, points, values);
        alreadyComputed[dimension] = true;
      }
    }
    // search for all particles in the cell and compute divergence using
    // them
    std::cerr << "************************************************\n";
    std::cerr << "DIVERGENCE ON CELL " << cellId << std::endl;
    divergenceField[cellId] = 0.;
    for (size_t d = 0; d < D; d++) {
      if (!alreadyComputed[d]) {
        std::vector<double> values;
        for (int i = 0; i < particleIds.size(); i++) {
          auto velocity =
              particles->getVector(velocityParticleFieldIds, particleIds[i]);
          Vector<double, D> direction = cellCenter - points[i];
          velocity.projectInto(direction);
          values.emplace_back(velocity[d]);
        }
        gradient[d] = rbf->gradientAt(cellCenter, d, points, values);
      }
      divergenceField[cellId] += gradient[d];
    }
    std::cerr << "Gradients: " << gradient[0] << " " << gradient[1]
              << std::endl;
    std::cerr << "DIVERGENCE VALUE " << divergenceField[cellId] << std::endl;
    std::cerr << "************************************************\n";
    // if (std::fabs(divergenceField[cellId]) < 1e-8)
    //   divergenceField[cellId] = 0.0;
  });
}

template <int D>
void SimulationSteps::applyFalloffToParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialFieldId,
    const std::vector<size_t> &velocityParticleFieldIds) {
  auto &materialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);

  structure->iterateCells([&](size_t cellId) {
    if (materialField[cellId] != Definitions::Material::FLUID)
      return;
    // check if cell touches solid
    auto neighbors = structure->cellNeighbors(cellId);
    auto cellRegion = structure->cellRegion(cellId);
    std::vector<double> gradient(D, 0.0);
    std::vector<size_t> particleIds;
    std::vector<Point<double, D>> points;
    particles->iterateParticles(
        structure->cellRegion(cellId), //.expanded(Vector<double, D>(1.5)),
        [&](size_t id, Point<double, D> p) {
          particleIds.emplace_back(id);
          points.emplace_back(p);
        });
    auto falloff = [](double d) {
      return d; /*-d * d + 2 * d;*/
    };
    for (auto neighbor : neighbors) {
      if (neighbor.id <= 0 ||
          materialField[neighbor.id] == Definitions::Material::SOLID) {
        std::vector<double> values(points.size(), 0.);
        switch (neighbor.side) {
        case Definitions::Side::LEFT:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((points[i].x() - cellRegion.lower().x()) /
                                   cellRegion.size().x()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[0], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[0],
            //  particleIds[i], value);
          }
          break;
        case Definitions::Side::RIGHT:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((cellRegion.upper().x() - points[i].x()) /
                                   cellRegion.size().x()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[0], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[0],
            //  particleIds[i], value);
          }
          break;
        case Definitions::Side::BOTTOM:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((points[i].y() - cellRegion.lower().y()) /
                                   cellRegion.size().y()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[1], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[1],
            //  particleIds[i], value);
          }
          break;
        case Definitions::Side::TOP:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((cellRegion.upper().y() - points[i].y()) /
                                   cellRegion.size().y()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[1], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[1],
            //  particleIds[i], value);
          }
          break;
        case Definitions::Side::BACK:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((points[i].z() - cellRegion.lower().z()) /
                                   cellRegion.size().z()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[2], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[2],
            //  particleIds[i], value);
          }
          break;
        case Definitions::Side::FRONT:
          for (size_t i = 0; i < particleIds.size(); i++) {
            double value = falloff((cellRegion.upper().z() - points[i].z()) /
                                   cellRegion.size().z()) *
                           particles->getScalarProperty(
                               velocityParticleFieldIds[2], particleIds[i]);
            // particles->setScalarProperty(velocityParticleFieldIds[2],
            //  particleIds[i], value);
          }
          break;
        }
      }
    }
  });
}

template <int D>
void SimulationSteps::solvePressure(StructureInterface<D> *structure,
                                    size_t pressureCellFieldId,
                                    size_t divergenceCellFieldId,
                                    size_t cellMaterialFieldId, double dt,
                                    DifferentialRBF<D> *rbf,
                                    bool adaptiveRBFRadius) {
  Timer assembly_timer;
  // setup field references
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &pressureField = *structure->template field<double>(pressureCellFieldId);
  // In order to setup the system, we need to track only the fluid cells
  size_t fluidCount = 0;
  std::map<size_t, size_t> indexMap;
  std::vector<int> fluidCells;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      fluidCount++;
      int lastId = indexMap.size();
      indexMap[cellId] = lastId;
      fluidCells.push_back(cellId);
    }
  });

  if (!fluidCount)
    return;
#ifdef FUROO_PROFILE
  CodeProfiler::count(concat("#fluidCells"), fluidCount);
#endif

  // Lets compute the RHS side of the system, which is = divergence / dt
  Eigen::VectorXd rhs;
  rhs = Eigen::VectorXd::Zero(fluidCount);

  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd pressure;
  std::cerr << "BUILD SYSTEM" << std::endl;

#pragma omp parallel
  {
    std::vector<Eigen::Triplet<double>> threadTriplets;
#pragma omp for
    for (int fi = 0; fi < fluidCells.size(); fi++) {
      size_t cellId = fluidCells[fi];
      rhs[indexMap[cellId]] = divergenceField[cellId] / dt;
      std::vector<size_t> neighbors;
      auto directNeighbors = structure->cellNeighbors(cellId);
      for (auto n : directNeighbors)
        if (n.id >= 0)
          //          &&             cellMaterialField[n.id] ==
          //          Definitions::Material::SOLID)
          neighbors.emplace_back(n.id);
      //      structure->iterateCellRing(cellId, [&](size_t nid) {
      //        if (cellMaterialField[nid] != Definitions::Material::SOLID)
      //          neighbors.emplace_back(nid);
      //      });
      //      auto weights = computePressureRingMatrixStencilWeights(structure,
      //      cellId,
      //                                                             rbf,
      //                                                             neighbors);
      std::vector<double> weights;
      {
        auto cellCenter = structure->cellCenterPosition(cellId);
        // First of all, we need to construct the set of stencil points
        std::vector<Point<double, D>> points;
        points.emplace_back(cellCenter);
        for (auto nid : neighbors) {
          points.emplace_back(structure->cellCenterPosition(nid));
        }
        weights = rbf->laplacianWeights(distance(points[0], points[1]));
      }
      double centerWeight = weights[0];
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] >= 0 &&
            cellMaterialField[neighbors[i]] == Definitions::Material::AIR)
          continue;
        if (neighbors[i] >= 0 &&
            cellMaterialField[neighbors[i]] == Definitions::Material::FLUID)

          threadTriplets.emplace_back(Eigen::Triplet<double>(
              indexMap[cellId], indexMap[neighbors[i]], weights[i + 1]));

        else
          centerWeight += weights[i + 1];
      }
      threadTriplets.emplace_back(Eigen::Triplet<double>(
          indexMap[cellId], indexMap[cellId], centerWeight));
    }
#pragma omp critical
    {
      // Join all triplets
      triplets.insert(triplets.end(), threadTriplets.begin(),
                      threadTriplets.end());
    }
  }
  std::cerr << "SOLVE SYSTEM" << std::endl;
  CodeProfiler::profile("PPE_Assembly", assembly_timer.ellapsedTimeInSeconds());
  // and finally solve the system...
  Timer timer;
  timer.reset();
  auto timenow =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::cout << ctime(&timenow) << std::endl;
  Eigen::setNbThreads(omp_get_max_threads());
  Eigen::SparseMatrix<double, Eigen::RowMajor> A(indexMap.size(),
                                                 indexMap.size());
  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
  // Eigen::Lower|Eigen::Upper> solver;
  solver.setTolerance(1e-6);
  solver.compute(A);
  pressure = solver.solve(rhs);

  std::cout << "Error: " << solver.error() << std::endl;
  std::cout << "Iterations: " << solver.iterations() << " out of "
            << solver.maxIterations() << std::endl;
  std::cout << "Epsilon: " << Eigen::NumTraits<double>::epsilon() << std::endl;

  if (pressure.isZero(1e-15)) {
    std::cerr << A << std::endl;
    std::cerr << rhs << std::endl;
  }
  CodeProfiler::profile("PPE", timer.ellapsedTimeInSeconds());
  std::cerr << "Psolver: " << timer.ellapsedTimeInSeconds() << std::endl;
  timenow =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::cout << ctime(&timenow) << std::endl;

#pragma omp parallel
  {
    iterateGridCells(structure, [&](size_t cellId) {
      pressureField[cellId] = 0.;
      if (cellMaterialField[cellId] == Definitions::Material::FLUID)
        pressureField[cellId] = pressure[indexMap[cellId]];
    });
  }
}

// BOTTLENECK
template <int D>
void SimulationSteps::solvePressureRing(StructureInterface<D> *structure,
                                        size_t pressureCellFieldId,
                                        size_t divergenceCellFieldId,
                                        size_t cellMaterialFieldId, double dt,
                                        DifferentialRBF<D> *rbf) {
  Timer assembly_timer;
  // setup field references
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &pressureField = *structure->template field<double>(pressureCellFieldId);
  // In order to setup the system, we need to track only the fluid cells
  size_t fluidCount = 0;
  std::map<size_t, size_t> indexMap;
  std::vector<int> fluidCells;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      fluidCount++;
      int lastId = indexMap.size();
      indexMap[cellId] = lastId;
      fluidCells.push_back(cellId);
    }
  });

  if (!fluidCount)
    return;
#ifdef FUROO_PROFILE
  CodeProfiler::count(concat("#fluidCells"), fluidCount);
#endif

  // Lets compute the RHS side of the system, which is = divergence / dt
  Eigen::VectorXd rhs;
  rhs = Eigen::VectorXd::Zero(fluidCount);

  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd pressure;
  std::cerr << "BUILD SYSTEM" << std::endl;

#pragma omp parallel
  {
    std::vector<Eigen::Triplet<double>> threadTriplets;
#pragma omp for
    for (int fi = 0; fi < fluidCells.size(); fi++) {
      size_t cellId = fluidCells[fi];
      rhs[indexMap[cellId]] = divergenceField[cellId] / dt;
      auto directNeighbors = structure->cellNeighbors(cellId);
      bool isRegular = false;

      if (directNeighbors.size() == 6) {
        // check if all neighbors have level equal to the cellId
        size_t cellLevel =
            dynamic_cast<CellGraph3 *>(structure)->node(cellId).level();
        bool allSameLevel = true;
        for (auto &neighbor : directNeighbors) {
          size_t neighLevel =
              dynamic_cast<CellGraph3 *>(structure)->node(neighbor.id).level();
          if (neighLevel != cellLevel)
            allSameLevel = false;
        }
        if (allSameLevel)
          isRegular = true;
      }
      // isRegular = false;

      std::vector<double> weights;
      std::vector<size_t> neighbors;
      if (!isRegular) {
        for (auto n : directNeighbors) {
          if (n.id >= 0 &&
              cellMaterialField[n.id] == Definitions::Material::SOLID)
            neighbors.emplace_back(n.id);
        }
        structure->iterateCellRing(cellId, [&](size_t nid) {
          if (cellMaterialField[nid] != Definitions::Material::SOLID)
            neighbors.emplace_back(nid);
        });
        weights = computePressureRingMatrixStencilWeights(structure, cellId,
                                                          rbf, neighbors);
      } else {
        auto cellCenter = structure->cellCenterPosition(cellId);
        // First of all, we need to construct the set of stencil points
        std::vector<Point<double, D>> points;
        points.emplace_back(cellCenter);
        for (auto neighbor : directNeighbors) {
          points.emplace_back(structure->cellCenterPosition(neighbor.id));
          neighbors.push_back(neighbor.id);
        }
        weights = rbf->laplacianWeights(distance(points[0], points[1]));
      }

      double centerWeight = weights[0];
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] >= 0 &&
            cellMaterialField[neighbors[i]] == Definitions::Material::AIR)
          continue;
        if (neighbors[i] >= 0 &&
            cellMaterialField[neighbors[i]] == Definitions::Material::FLUID)

          threadTriplets.emplace_back(Eigen::Triplet<double>(
              indexMap[cellId], indexMap[neighbors[i]], weights[i + 1]));

        else
          centerWeight += weights[i + 1];
      }
      threadTriplets.emplace_back(Eigen::Triplet<double>(
          indexMap[cellId], indexMap[cellId], centerWeight));
    }
#pragma omp critical
    {
      // Join all triplets
      triplets.insert(triplets.end(), threadTriplets.begin(),
                      threadTriplets.end());
    }
  }
  std::cerr << "SOLVE SYSTEM" << std::endl;
  CodeProfiler::profile("PPE_Assembly", assembly_timer.ellapsedTimeInSeconds());
  // and finally solve the system...
  Timer timer;
  timer.reset();
  auto timenow =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::cout << ctime(&timenow) << std::endl;
  Eigen::setNbThreads(omp_get_max_threads());
  Eigen::SparseMatrix<double, Eigen::RowMajor> A(indexMap.size(),
                                                 indexMap.size());
  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
  solver.setTolerance(1e-6);
  solver.compute(A);
  pressure = solver.solve(rhs);

  std::cout << "Tolerance: " << solver.error() << std::endl;
  std::cout << "Iterations: " << solver.iterations() << " out of "
            << solver.maxIterations() << std::endl;
  std::cout << "Epsilon: " << Eigen::NumTraits<double>::epsilon() << std::endl;

  if (pressure.isZero(1e-15)) {
    std::cerr << A << std::endl;
    std::cerr << rhs << std::endl;
  }
  CodeProfiler::profile("PPE", timer.ellapsedTimeInSeconds());
  std::cerr << "Psolver: " << timer.ellapsedTimeInSeconds() << std::endl;
  timenow =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::cout << ctime(&timenow) << std::endl;

#pragma omp parallel
  {
    iterateGridCells(structure, [&](size_t cellId) {
      pressureField[cellId] = 0.;
      if (cellMaterialField[cellId] == Definitions::Material::FLUID)
        pressureField[cellId] = pressure[indexMap[cellId]];
    });
  }
}

template <int D>
void SimulationSteps::solvePressure(
    StructureInterface<D> *structure,
    const ParticleSystem<double, D> *particles, size_t pressureCellFieldId,
    size_t divergenceCellFieldId, size_t cellMaterialFieldId,
    size_t surfaceMaskFieldId, size_t cellParticleCenterFieldId, double dt,
    DifferentialRBF<D> *rbf, bool adaptiveRBFRadius) {
  // setup field references
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &divergenceField =
      *structure->template field<double>(divergenceCellFieldId);
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskFieldId);
  auto &pressureField = *structure->template field<double>(pressureCellFieldId);
  // In order to setup the system, we need to track only the fluid cells
  size_t fluidCount = 0;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      fluidCount++;
      auto neighborCells = structure->cellNeighbors(cellId);
    }
  });
  if (!fluidCount)
    return;
  // Since the matrix will have rows only to fluid cells, we need a map
  // between cell id and matrix row id
  // structure index -> matrix index
  std::map<size_t, size_t> indexMap;
  std::function<size_t(size_t)> cellIndex = [&](size_t cellId) {
    auto it = indexMap.find(cellId);
    if (it == indexMap.end()) {
      size_t nextIndex = indexMap.size();
      indexMap[cellId] = nextIndex;
      return indexMap.size() - 1;
    }
    return it->second;
  };

  // Lets compute the RHS side of the system, which is = divergence / dt
  Eigen::VectorXd rhs;
  rhs = Eigen::VectorXd::Zero(fluidCount); // +1 is for air cells
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      rhs[cellIndex(cellId)] = divergenceField[cellId] / dt;
  });
  // And now the triplets that form the matrix
  std::vector<Eigen::Triplet<double>> triplets;
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    auto neighbors = structure->cellNeighbors(cellId);
    auto weights = computePressureMatrixStencilWeights(
        structure, cellId, cellMaterialFieldId, cellParticleCenterFieldId, rbf,
        adaptiveRBFRadius);
    double centerWeight = weights[0];
    for (size_t i = 0; i < neighbors.size(); i++) {
      if (neighbors[i].id >= 0 &&
          cellMaterialField[neighbors[i].id] == Definitions::Material::AIR) {
        continue;
        triplets.emplace_back(Eigen::Triplet<double>(
            cellIndex(cellId), cellIndex(neighbors[i].id), weights[i + 1]));
        triplets.emplace_back(cellIndex(neighbors[i].id),
                              cellIndex(neighbors[i].id), 1);
      }
      // continue;
      else if (neighbors[i].id >= 0 && cellMaterialField[neighbors[i].id] ==
                                           Definitions::Material::FLUID)
        triplets.emplace_back(Eigen::Triplet<double>(
            cellIndex(cellId), cellIndex(neighbors[i].id), weights[i + 1]));
      else {
        centerWeight += weights[i + 1];
#ifdef FUROO_PROFILE
        CodeProfiler::count("#Neumann");
#endif
      }
    }
    triplets.emplace_back(Eigen::Triplet<double>(
        cellIndex(cellId), cellIndex(cellId), centerWeight));
  });

  // and finally solve the system...
  Eigen::SparseMatrix<double> A(indexMap.size(), indexMap.size());
  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver;
  solver.compute(A);
  Eigen::VectorXd pressure;
  pressure = solver.solve(rhs);
  structure->iterateCells([&](size_t cellId) {
    pressureField[cellId] = 0.;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      pressureField[cellId] = pressure[cellIndex(cellId)];
  });
}

template <int D>
double SimulationSteps::cfl(const ParticleSystem<double, D> *particleSystem,
                            const std::vector<size_t> &fieldIds, double dx) {
  double maxVelocity = 0.;
  // for (size_t d = 0; d < D; d++) {
  //   maxVelocity = std::max(
  //       maxVelocity,
  //       std::fabs(particleSystem->maxPropertyValue(fieldIds[d])));
  //   maxVelocity = std::max(
  //       maxVelocity,
  //       std::fabs(particleSystem->minPropertyValue(fieldIds[d])));
  // }
  particleSystem->iterateParticles([&](size_t id, Point<double, D> p) {
    auto v = particleSystem->getVector(fieldIds, id);
    maxVelocity = std::max(maxVelocity, v.length());
  });
  if (maxVelocity == 0.)
    return INFINITY;
  return dx / maxVelocity;
}

template <int D>
void SimulationSteps::markCells(StructureInterface<D> *structure,
                                size_t cellMaterialFieldId,
                                ParticleSystem<double, D> *particles) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  structure->iterateCells([&](size_t cellId) {
    auto cellLevel = [&](size_t cellId) -> size_t {
      if (D == 2)
        return dynamic_cast<CellGraph2 *>(structure)->node(cellId).level();
      return dynamic_cast<CellGraph3 *>(structure)->node(cellId).level();
    };
    size_t level = cellLevel(cellId);

    cellMaterialField[cellId] =
        // (particles->containsParticle(structure->cellRegion(cellId)))
        (particles->containsParticleZ(structure->cellRegion(cellId), level))
            ? Definitions::Material::FLUID
            : Definitions::Material::AIR;
  });
}

template <>
inline void SimulationSteps::refineFluid(ParticleSystem<double, 2> *particles,
                                         StructureInterface<2> *structure,
                                         size_t cellMaterialFieldId,
                                         size_t fluidLevel) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  std::function<bool(size_t, size_t)> comp = [&](size_t cellA,
                                                 size_t cellB) -> bool {
    if (cellMaterialField[cellA] == cellMaterialField[cellB])
      return cellGraph->node(cellA).level() > cellGraph->node(cellB).level();
    return cellMaterialField[cellA] != Definitions::Material::FLUID;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  q.push(1);
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    auto curLevel = cellGraph->node(curCellId).level();
    if (cellMaterialField[curCellId] == Definitions::Material::FLUID) {
      if (curLevel >= fluidLevel)
        continue;
      auto children = cellGraph->refine(curCellId);
      for (auto child : children) {
        if (!child)
          continue;
        q.push(child);
        cellMaterialField[child] =
            // (particles->containsParticle(cellGraph->cellRegion(child)))
            (particles->containsParticleZ(cellGraph->cellRegion(child),
                                          cellGraph->node(child).level()))
                ? Definitions::Material::FLUID
                : Definitions::Material::AIR;
      }
    }
  }
}
// TODO: make cellGraph template!
template <>
inline void SimulationSteps::refineFluid(ParticleSystem<double, 3> *particles,
                                         StructureInterface<3> *structure,
                                         size_t cellMaterialFieldId,
                                         size_t fluidLevel) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  std::function<bool(size_t, size_t)> comp = [&](size_t cellA,
                                                 size_t cellB) -> bool {
    if (cellMaterialField[cellA] == cellMaterialField[cellB])
      return cellGraph->node(cellA).level() > cellGraph->node(cellB).level();
    return cellMaterialField[cellA] != Definitions::Material::FLUID;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  q.push(1);
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    auto curLevel = cellGraph->node(curCellId).level();
    if (cellMaterialField[curCellId] == Definitions::Material::FLUID) {
      if (curLevel >= fluidLevel)
        continue;
      auto children = cellGraph->refine(curCellId);
      for (auto child : children) {
        if (!child)
          continue;
        q.push(child);
        cellMaterialField[child] =
            // (particles->containsParticle(cellGraph->cellRegion(child)))
            (particles->containsParticleZ(cellGraph->cellRegion(child),
                                          cellGraph->node(child).level()))
                ? Definitions::Material::FLUID
                : Definitions::Material::AIR;
      }
    }
  }
}

template <>
inline std::set<size_t>
SimulationSteps::markFluidSurface(StructureInterface<2> *structure,
                                  size_t surfaceMaskCellFieldId,
                                  size_t cellMaterialFieldId) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  std::set<size_t> surfaceCells;
  cellGraph->iterateFaces([&](size_t faceId) {
    auto edge = cellGraph->edge(faceId);
    auto nodes = edge.nodes();
    if (!nodes[0] &&
        cellMaterialField[nodes[1]] == Definitions::Material::FLUID) {
      surfaceField[nodes[1]] = 1;
      surfaceCells.insert(nodes[1]);
    } else if (!nodes[1] &&
               cellMaterialField[nodes[0]] == Definitions::Material::FLUID) {
      surfaceField[nodes[0]] = 1;
      surfaceCells.insert(nodes[0]);
    } else if (cellMaterialField[nodes[0]] == Definitions::Material::FLUID &&
               cellMaterialField[nodes[1]] != Definitions::Material::FLUID) {
      surfaceField[nodes[0]] = 1;
      surfaceCells.insert(nodes[0]);
    } else if (cellMaterialField[nodes[0]] != Definitions::Material::FLUID &&
               cellMaterialField[nodes[1]] == Definitions::Material::FLUID) {
      surfaceField[nodes[1]] = 1;
      surfaceCells.insert(nodes[1]);
    }
  });
  return surfaceCells;
}

template <>
inline std::set<size_t>
SimulationSteps::markFluidSurface(StructureInterface<3> *structure,
                                  size_t surfaceMaskCellFieldId,
                                  size_t cellMaterialFieldId) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  std::set<size_t> surfaceCells;
  cellGraph->iterateFaces([&](size_t faceId) {
    auto edge = cellGraph->edge(faceId);
    auto nodes = edge.nodes();
    if (!nodes[0] &&
        cellMaterialField[nodes[1]] == Definitions::Material::FLUID) {
      surfaceField[nodes[1]] = 1;
      surfaceCells.insert(nodes[1]);
    } else if (!nodes[1] &&
               cellMaterialField[nodes[0]] == Definitions::Material::FLUID) {
      surfaceField[nodes[0]] = 1;
      surfaceCells.insert(nodes[0]);
    } else if (cellMaterialField[nodes[0]] == Definitions::Material::FLUID &&
               cellMaterialField[nodes[1]] != Definitions::Material::FLUID) {
      surfaceField[nodes[0]] = 1;
      surfaceCells.insert(nodes[0]);
    } else if (cellMaterialField[nodes[0]] != Definitions::Material::FLUID &&
               cellMaterialField[nodes[1]] == Definitions::Material::FLUID) {
      surfaceField[nodes[1]] = 1;
      surfaceCells.insert(nodes[1]);
    }
  });
  return surfaceCells;
}

template <>
inline std::queue<size_t> SimulationSteps::updateSurfaceCells(
    ParticleSystem<double, 2> *particles, StructureInterface<2> *structure,
    size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
    size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // First, gather old surface cells and find out new surface cells
  std::queue<size_t> surfaceCells;
  std::vector<size_t> oldSurface;
  // surface candidates are fluid cells that may become surface cells
  std::set<size_t> candidates;
  cellGraph->iterateCells([&](size_t id) {
    if (surfaceField[id]) {
      cellGraph->iterateCellRing(id, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::AIR) {
          // if (particles->containsParticle(cellGraph->cellRegion(nid))) {
          if (particles->containsParticleZ(cellGraph->cellRegion(nid),
                                           cellGraph->node(nid).level())) {
            cellMaterialField[nid] = Definitions::Material::FLUID;
            surfaceCells.push(nid);
            surfaceField[nid] = 1;
          }
        } else if (!surfaceField[nid] &&
                   cellMaterialField[nid] == Definitions::Material::FLUID) {
          candidates.insert(nid);
        }
      });
      // check if this surface cell remains
      // if (!particles->containsParticle(cellGraph->cellRegion(id))) {
      if (!particles->containsParticleZ(cellGraph->cellRegion(id),
                                        cellGraph->node(id).level())) {
        cellMaterialField[id] = Definitions::Material::AIR;
        surfaceField[id] = 0;
      } else
        oldSurface.emplace_back(id);
    }
  });
  // Now check if old surface cells remain
  for (auto c : oldSurface) {
    bool isSurface = false;
    size_t ringSize = 0;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      ringSize++;
      if (cellMaterialField[id] != Definitions::Material::FLUID)
        isSurface = true;
    });
    if (!isSurface && ringSize == 8)
      surfaceField[c] = 0;
    else if (isSurface || ringSize < 8) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  // Now check candidate cells
  for (auto c : candidates) {
    bool isSurface = false;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      if (cellMaterialField[id] != Definitions::Material::FLUID)
        isSurface = true;
    });
    if (isSurface) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  /*
    // iterate air cells fixing surface cells
    cellGraph->iterateCells([&](size_t id) {
      if (cellMaterialField[id] == Definitions::Material::AIR) {
        // TODO when we allow the existence of injectors, their region must be
    refined too!
        // We don't need to look for particles on coarse air cells, as the air
    near fluid
        // is refined and the cfl don't permit any particles to jump into a
    coarse air cell if (cellGraph->node(id).level() < surfaceLevel) return;
    bool isFluidCell = false;
    particles->iterateParticles(cellGraph->cellRegion(id),
                                    [&](size_t pid, Point<double, 2> p) {
                                      isFluidCell = true;
                                      UNUSED_VARIABLE(pid);
                                      UNUSED_VARIABLE(p);
                                    });
        if (!isFluidCell)
          return;
        cellMaterialField[id] = Definitions::Material::FLUID;
        bool amIASurface = false;
        for (auto n : cellGraph->neighbors(id)) {
          if (n.id > 0 && cellMaterialField[n.id] ==
    Definitions::Material::AIR) amIASurface = true; if (n.id > 0 &&
    surfaceField[n.id]) { bool isSurface = false; for (auto nn :
    cellGraph->neighbors(n.id)) if (nn.id > 0
                  && cellMaterialField[nn.id] == Definitions::Material::AIR) {
                isSurface = true;
                break;
              }
            if (!isSurface)
              surfaceField[n.id] = 0;
          }
        }
        if (amIASurface) {
          surfaceField[id] = 1;
          surfaceCells.push(id);
        }
      } else if (surfaceField[id]) {
        bool isAirCell = true;
        particles->iterateParticles(cellGraph->cellRegion(id),
                                    [&](size_t pid, Point<double, 2> p) {
                                      isAirCell = false;
                                      UNUSED_VARIABLE(pid);
                                      UNUSED_VARIABLE(p);
                                    });
        if (isAirCell) {
          surfaceField[id] = 0;
          cellMaterialField[id] = Definitions::Material::AIR;
          for (auto n : cellGraph->neighbors(id))
            if (n.id > 0
                && cellMaterialField[n.id] == Definitions::Material::FLUID) {
              surfaceField[n.id] = 1;
              surfaceCells.push(n.id);
            }
        } else {
          surfaceCells.push(id);
        }
      }
    });*/
  return surfaceCells;
}

template <>
inline std::queue<size_t> SimulationSteps::updateSurfaceCells(
    ParticleSystem<double, 3> *particles, StructureInterface<3> *structure,
    size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
    size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // First, gather old surface cells and find out new surface cells
  std::queue<size_t> surfaceCells;
  std::vector<size_t> oldSurface;
  // surface candidates are fluid cells that may become surface cells
  std::set<size_t> candidates;
  cellGraph->iterateCells([&](size_t id) {
    if (surfaceField[id]) {
      cellGraph->iterateCellRing(id, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::AIR) {
          // if (particles->containsParticle(cellGraph->cellRegion(nid))) {
          if (particles->containsParticleZ(cellGraph->cellRegion(nid),
                                           cellGraph->node(nid).level())) {
            cellMaterialField[nid] = Definitions::Material::FLUID;
            // surfaceCells.push(nid);
            // surfaceField[nid] = 1;
            candidates.insert(nid);
          }
        } else if (!surfaceField[nid] &&
                   cellMaterialField[nid] == Definitions::Material::FLUID) {
          candidates.insert(nid);
        }
      });
      // check if this surface cell remains
      // if (!particles->containsParticle(cellGraph->cellRegion(id))) {
      if (!particles->containsParticleZ(cellGraph->cellRegion(id),
                                        cellGraph->node(id).level())) {
        cellMaterialField[id] = Definitions::Material::AIR;
        surfaceField[id] = 0;
      } else
        oldSurface.emplace_back(id);
    }
  });
  // Now check if old surface cells remain
  for (auto c : oldSurface) {
    bool isSurface = false;
    size_t ringSize = 0;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      ringSize++;
      if (cellMaterialField[id] != Definitions::Material::FLUID)
        isSurface = true;
    });
    if (!isSurface /*&& ringSize == 26*/)
      surfaceField[c] = 0;
    else if (isSurface /* || ringSize < 26*/)
      surfaceCells.push(c);
  }
  // Now check candidate cells
  for (auto c : candidates) {
    bool isSurface = false;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      if (cellMaterialField[id] != Definitions::Material::FLUID)
        isSurface = true;
    });
    if (isSurface) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  return surfaceCells;
}

template <>
inline std::queue<size_t> SimulationSteps::updateSurfaceCellsCoarseBottom(
    ParticleSystem<double, 3> *particles, StructureInterface<3> *structure,
    size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
    size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // First, gather old surface cells and find out new surface cells
  std::queue<size_t> surfaceCells;
  std::vector<size_t> oldSurface;
  // surface candidates are fluid cells that may become surface cells
  std::set<size_t> candidates;
  cellGraph->iterateCells([&](size_t id) {
    if (surfaceField[id]) {
      cellGraph->iterateCellRing(id, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::AIR) {
          if (particles->containsParticle(cellGraph->cellRegion(nid))) {
            cellMaterialField[nid] = Definitions::Material::FLUID;
            surfaceCells.push(nid);
            surfaceField[nid] = 1;
          }
        } else if (!surfaceField[nid] &&
                   cellMaterialField[nid] == Definitions::Material::FLUID) {
          candidates.insert(nid);
        }
      });
      // check if this surface cell remains
      if (!particles->containsParticle(cellGraph->cellRegion(id))) {
        cellMaterialField[id] = Definitions::Material::AIR;
        surfaceField[id] = 0;
      } else
        oldSurface.emplace_back(id);
    }
  });
  // Now check if old surface cells remain
  for (auto c : oldSurface) {
    bool isSurface = false;
    size_t ringSize = 0;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      ringSize++;
      if (cellMaterialField[id] == Definitions::Material::AIR)
        isSurface = true;
    });
    if (!isSurface && ringSize == 26)
      surfaceField[c] = 0;
    else if (isSurface || ringSize < 26) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  // Now check candidate cells
  for (auto c : candidates) {
    bool isSurface = false;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      if (cellMaterialField[id] == Definitions::Material::AIR)
        isSurface = true;
    });
    if (isSurface) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  return surfaceCells;
}

template <>
inline std::queue<size_t> SimulationSteps::updateSurfaceCellsCoarseBottom(
    ParticleSystem<double, 2> *particles, StructureInterface<2> *structure,
    size_t surfaceMaskFieldId, size_t cellMaterialFieldId,
    size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // First, gather old surface cells and find out new surface cells
  std::queue<size_t> surfaceCells;
  std::vector<size_t> oldSurface;
  // surface candidates are fluid cells that may become surface cells
  std::set<size_t> candidates;
  cellGraph->iterateCells([&](size_t id) {
    if (surfaceField[id]) {
      cellGraph->iterateCellRing(id, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::AIR) {
          if (particles->containsParticle(cellGraph->cellRegion(nid))) {
            cellMaterialField[nid] = Definitions::Material::FLUID;
            oldSurface.emplace_back(nid);
            // surfaceCells.push(nid);
            // surfaceField[nid] = 1;
          }
        } else if (!surfaceField[nid] &&
                   cellMaterialField[nid] == Definitions::Material::FLUID) {
          candidates.insert(nid);
        }
      });
      // check if this surface cell remains
      if (!particles->containsParticle(cellGraph->cellRegion(id))) {
        cellMaterialField[id] = Definitions::Material::AIR;
        surfaceField[id] = 0;
      } else
        oldSurface.emplace_back(id);
    }
  });
  // Now check if old surface cells remain
  for (auto c : oldSurface) {
    bool airSurface = false, solidSurface = false;
    size_t ringSize = 0;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      ringSize++;
      if (cellMaterialField[id] == Definitions::Material::AIR)
        airSurface = true;
      else if (cellMaterialField[id] == Definitions::Material::SOLID)
        solidSurface = true;
    });
    if (!airSurface && ringSize == 8)
      surfaceField[c] = 0;
    else if (!airSurface && solidSurface)
      surfaceField[c] = 0;
    else if (airSurface || ringSize < 8) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  // Now check candidate cells
  for (auto c : candidates) {
    bool airSurface = false;
    cellGraph->iterateCellRing(c, [&](size_t id) {
      if (cellMaterialField[id] == Definitions::Material::AIR)
        airSurface = true;
    });
    if (airSurface) {
      surfaceCells.push(c);
      surfaceField[c] = 1;
    }
  }
  return surfaceCells;
}

template <>
inline void SimulationSteps::coarseGraph(
    StructureInterface<2> *structure, size_t surfaceMaskCellFieldId,
    size_t cellMaterialFieldId,
    const std::function<bool(Definitions::Material, size_t)> &doCoarse) {
  auto *graph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField =
      *graph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *graph->template field<Definitions::Material>(cellMaterialFieldId);
  std::function<bool(size_t, size_t)> comp = [&](size_t cellA,
                                                 size_t cellB) -> bool {
    return graph->node(cellA).level() < graph->node(cellB).level();
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  graph->iterateCells([&](size_t c) { q.push(c); });
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    if (graph->node(curCellId).isValid()) {
      // TODO: do we need really need to check siblings before?
      auto siblings = graph->nodeSiblings(curCellId);
      if (siblings.size() != 4)
        continue;
      bool coarse = true;
      auto nodeMaterialType = cellMaterialField[siblings[0]];
      for (auto s : siblings)
        if (surfaceField[s] || cellMaterialField[s] != nodeMaterialType)
          coarse = false;
      if (coarse) {
        if (doCoarse(nodeMaterialType, graph->node(curCellId).level())) {
          try {
            size_t newNode = graph->coarse(curCellId);
            cellMaterialField[newNode] = nodeMaterialType;
            if (newNode)
              q.push(newNode);
          } catch (std::string e) {
            std::ostringstream errorMessage;
            errorMessage
                << "SimulationSteps::coarseGraph failed to coarse node "
                << curCellId << " with region " << graph->cellRegion(curCellId);
            THROW(false, errorMessage.str());
          }
        }
      }
    }
  }
}

template <>
inline void SimulationSteps::coarseGraph(
    StructureInterface<3> *structure, size_t surfaceMaskCellFieldId,
    size_t cellMaterialFieldId,
    const std::function<bool(Definitions::Material, size_t)> &doCoarse) {
  auto *graph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField =
      *graph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *graph->template field<Definitions::Material>(cellMaterialFieldId);
  std::function<bool(size_t, size_t)> comp = [&](size_t cellA,
                                                 size_t cellB) -> bool {
    return graph->node(cellA).level() < graph->node(cellB).level();
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  graph->iterateCells([&](size_t c) { q.push(c); });
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    if (graph->node(curCellId).isValid()) {
      // TODO: do we need really need to check siblings before?
      auto siblings = graph->nodeSiblings(curCellId);
      if (siblings.size() != 8)
        continue;
      bool coarse = true;
      auto nodeMaterialType = cellMaterialField[siblings[0]];
      for (auto s : siblings)
        if (surfaceField[s] || cellMaterialField[s] != nodeMaterialType)
          coarse = false;
      if (coarse) {
        if (doCoarse(nodeMaterialType, graph->node(curCellId).level())) {
          try {
            size_t newNode = graph->coarse(curCellId);
            cellMaterialField[newNode] = nodeMaterialType;
            if (newNode)
              q.push(newNode);
          } catch (std::string e) {
            std::ostringstream errorMessage;
            errorMessage
                << "SimulationSteps::coarseGraph failed to coarse node "
                << curCellId << " with region " << graph->cellRegion(curCellId);
            THROW(false, errorMessage.str());
          }
        }
      }
    }
  }
}

template <>
inline std::queue<size_t> SimulationSteps::refineNearSurface(
    StructureInterface<2> *structure, size_t surfaceMaskCellFieldId,
    size_t cellSurfaceTDistanceFieldId, size_t cellMaterialFieldId,
    std::queue<size_t> &surfaceCells) {
  if (surfaceCells.empty())
    return std::queue<size_t>();
  size_t narrowBandWidth = 2;
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &distanceField =
      *cellGraph->template field<int>(cellSurfaceTDistanceFieldId);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // first we mark all surface cells with topological distance to 0 and the
  // rest as infinity. The idea is to refine all cells with distance <
  // narrowBandWidth
  distanceField.setAll(INT_INFINITY);
  size_t surfaceLevel = 0;
  std::queue<size_t> narrowBand;
  std::queue<size_t> newSurface;
  while (!surfaceCells.empty()) {
    auto surfaceCell = surfaceCells.front();
    newSurface.push(surfaceCell);
    narrowBand.push(surfaceCell);
    surfaceCells.pop();
    distanceField[surfaceCell] = 0;
    // TODO: does this throw conditionis really necessary? Or the following if
    // statement makes sense?
    //    THROW(surfaceField[surfaceCell] == 1,
    //          "SimulationSteps::refineNearSurface cell of surface queue not
    //          marked as surface");
    if (!surfaceField[surfaceCell])
      continue;
    surfaceLevel = cellGraph->node(surfaceCell).level();
  }
  THROW(surfaceLevel > 0,
        "SimulationSteps::refineNearSurface could not define surface level");
  // The bfs works like this:
  // For a given cell, tries to update the ring neighborhood with its distance
  // In the case a big neighbor exists, it must be refined
  // refine near surface
  while (!narrowBand.empty()) {
    auto curCellId = narrowBand.front();
    narrowBand.pop();
    if (!cellGraph->node(curCellId).isValid() ||
        distanceField[curCellId] >= narrowBandWidth)
      continue;
    bool returnToQueue = false;
    cellGraph->iterateCellRing(curCellId, [&](size_t neighborId) {
      if (cellGraph->node(neighborId).level() < surfaceLevel) {
        auto children = cellGraph->refine(neighborId);
        auto parentMaterial = cellMaterialField[neighborId];
        for (auto c : children)
          cellMaterialField[c] = parentMaterial;
        returnToQueue = true;
        return;
      }
      if (distanceField[neighborId] == INT_INFINITY)
        narrowBand.push(neighborId);
      distanceField[neighborId] =
          std::min(distanceField[neighborId], distanceField[curCellId] + 1);
    });
    if (returnToQueue)
      narrowBand.push(curCellId);
  }
  return newSurface;
}

template <>
inline std::queue<size_t> SimulationSteps::refineNearSurface(
    StructureInterface<3> *structure, size_t surfaceMaskCellFieldId,
    size_t cellSurfaceTDistanceFieldId, size_t cellMaterialFieldId,
    std::queue<size_t> &surfaceCells) {
  if (surfaceCells.empty())
    return std::queue<size_t>();
  size_t narrowBandWidth = 2;
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &distanceField =
      *cellGraph->template field<int>(cellSurfaceTDistanceFieldId);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  // first we mark all surface cells with topological distance to 0 and the
  // rest as infinity. The idea is to refine all cells with distance <
  // narrowBandWidth
  distanceField.setAll(INT_INFINITY);
  size_t surfaceLevel = 0;
  std::queue<size_t> narrowBand;
  std::queue<size_t> newSurface;
  while (!surfaceCells.empty()) {
    auto surfaceCell = surfaceCells.front();
    newSurface.push(surfaceCell);
    narrowBand.push(surfaceCell);
    surfaceCells.pop();
    distanceField[surfaceCell] = 0;
    // TODO: does this throw conditionis really necessary? Or the following if
    // statement makes sense?
    //    THROW(surfaceField[surfaceCell] == 1,
    //          "SimulationSteps::refineNearSurface cell of surface queue not
    //          marked as surface");
    if (!surfaceField[surfaceCell])
      continue;
    surfaceLevel = cellGraph->node(surfaceCell).level();
  }
  THROW(surfaceLevel > 0,
        "SimulationSteps::refineNearSurface could not define surface level");
  // The bfs works like this:
  // For a given cell, tries to update the ring neighborhood with its distance
  // In the case a big neighbor exists, it must be refined
  // refine near surface
  while (!narrowBand.empty()) {
    auto curCellId = narrowBand.front();
    narrowBand.pop();
    if (!cellGraph->node(curCellId).isValid() ||
        distanceField[curCellId] >= narrowBandWidth)
      continue;
    bool returnToQueue = false;
    cellGraph->iterateCellRing(curCellId, [&](size_t neighborId) {
      if (cellGraph->node(neighborId).level() < surfaceLevel) {
        auto children = cellGraph->refine(neighborId);
        auto parentMaterial = cellMaterialField[neighborId];
        for (auto c : children)
          cellMaterialField[c] = parentMaterial;
        returnToQueue = true;
        return;
      }
      if (distanceField[neighborId] == INT_INFINITY)
        narrowBand.push(neighborId);
      distanceField[neighborId] =
          std::min(distanceField[neighborId], distanceField[curCellId] + 1);
    });
    if (returnToQueue)
      narrowBand.push(curCellId);
  }
  return newSurface;
}

template <>
inline void SimulationSteps::refineNearSolid(StructureInterface<2> *structure,
                                             size_t cellMaterialFieldId,
                                             size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);

  std::queue<size_t> solidCells;
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::SOLID)
      solidCells.push(cellId);
  });
  while (!solidCells.empty()) {
    auto cellId = solidCells.front();
    solidCells.pop();
    auto neighbors = cellGraph->cellNeighbors(cellId);
    bool shouldRefine = false;
    bool comeBack = false;
    for (auto neighbor : neighbors) {
      if (neighbor.id > 0 &&
          cellMaterialField[neighbor.id] == Definitions::Material::FLUID) {
        // auto neighborLevel = cellGraph->node(neighbor.id).level();.
        if (cellGraph->node(cellId).level() < surfaceLevel)
          shouldRefine = true;
        if (cellGraph->node(neighbor.id).level() < surfaceLevel) {
          comeBack = true;
          try {
            auto children = cellGraph->refine(neighbor.id);
            for (auto child : children)
              cellMaterialField[child] = Definitions::Material::FLUID;

          } catch (std::string &e) {
            std::cerr << e << std::endl;
            std::cerr << "Solid id: " << cellId << std::endl;
            THROW(false,
                  "SimulationsSteps::refineNearSolid not able to refine");
          }
        }
      }
    }
    if (shouldRefine) {
      auto children = cellGraph->refine(cellId);
      for (auto child : children) {
        cellMaterialField[child] = Definitions::Material::SOLID;
        solidCells.push(child);
      }
    } else if (comeBack)
      solidCells.push(cellId);
  }
}

template <>
inline void SimulationSteps::refineNearSolid(StructureInterface<3> *structure,
                                             size_t cellMaterialFieldId,
                                             size_t surfaceLevel) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);

  std::queue<size_t> solidCells;
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::SOLID)
      solidCells.push(cellId);
  });
  while (!solidCells.empty()) {
    auto cellId = solidCells.front();
    solidCells.pop();
    auto neighbors = cellGraph->cellNeighbors(cellId);
    bool shouldRefine = false;
    bool comeBack = false;
    for (auto neighbor : neighbors) {
      if (neighbor.id > 0 &&
          cellMaterialField[neighbor.id] == Definitions::Material::FLUID) {
        // auto neighborLevel = cellGraph->node(neighbor.id).level();.
        if (cellGraph->node(cellId).level() < surfaceLevel)
          shouldRefine = true;
        if (cellGraph->node(neighbor.id).level() < surfaceLevel) {
          comeBack = true;
          try {
            auto children = cellGraph->refine(neighbor.id);
            for (auto child : children)
              cellMaterialField[child] = Definitions::Material::FLUID;

          } catch (std::string &e) {
            std::cerr << e << std::endl;
            std::cerr << "Solid id: " << cellId << std::endl;
            THROW(false,
                  "SimulationsSteps::refineNearSolid not able to refine");
          }
        }
      }
    }
    if (shouldRefine) {
      auto children = cellGraph->refine(cellId);
      for (auto child : children) {
        cellMaterialField[child] = Definitions::Material::SOLID;
        solidCells.push(child);
      }
    } else if (comeBack)
      solidCells.push(cellId);
  }
}

// DEPRECATED
template <>
inline void SimulationSteps::fixSurface(StructureInterface<2> *structure,
                                        size_t surfaceMaskCellFieldId,
                                        size_t cellMaterialFieldId,
                                        std::queue<size_t> &surfaceCells) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField = *cellGraph->field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  while (!surfaceCells.empty()) {
    size_t curCell = surfaceCells.front();
    surfaceCells.pop();
    // for each fluid neighbor of a surface cell, check if all its ring is
    // fluid
    auto neighbors = cellGraph->cellNeighbors(curCell);
    for (auto n : neighbors) {
      // we don't have to process surfaces
      if (n.id < 0 || surfaceField[n.id] ||
          cellMaterialField[n.id] != Definitions::Material::FLUID)
        continue;
      bool allFluid = true;
      cellGraph->iterateCellRing(n.id, [&](size_t i) {
        if (cellMaterialField[i] != Definitions::Material::FLUID)
          allFluid = false;
      });
      if (!allFluid)
        surfaceField[n.id] = 1;
      else
        surfaceField[n.id] = 0;
    }
  }
}

// DEPRECATED
template <>
inline void SimulationSteps::fixSurface(StructureInterface<3> *structure,
                                        size_t surfaceMaskCellFieldId,
                                        size_t cellMaterialFieldId,
                                        std::queue<size_t> &surfaceCells) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField = *cellGraph->field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  while (!surfaceCells.empty()) {
    size_t curCell = surfaceCells.front();
    surfaceCells.pop();
    /* LEFT BOTTOM RIGHT TOP*/
    int neighborsId[4] = {0, 0, 0, 0};
    auto neighbors = cellGraph->cellNeighbors(curCell);
    for (auto n : neighbors)
      neighborsId[static_cast<size_t>(n.side)] = n.id;
    Definitions::Side sides[4] = {
        Definitions::Side::LEFT, Definitions::Side::BOTTOM,
        Definitions::Side::RIGHT, Definitions::Side::TOP};
    Definitions::Material materials[4];
    bool allFluid = true;
    for (size_t i = 0; i < 4; i++) {
      materials[i] =
          (neighborsId[static_cast<size_t>(sides[i])] > 0)
              ? cellMaterialField[neighborsId[static_cast<size_t>(sides[i])]]
              : Definitions::Material::AIR;
      if (materials[i] != Definitions::Material::FLUID)
        allFluid = false;
    }
    if (allFluid)
      surfaceField[curCell] = 0;
    for (size_t i = 0; i < 4; i++) {
      if (sides[i] != Definitions::oppositeSide(sides[(i + 1) % 4])) {
        if (materials[i] == Definitions::Material::FLUID &&
            materials[(i + 1) % 4] != Definitions::Material::FLUID) {
          surfaceField[neighborsId[static_cast<size_t>(sides[i])]] = 1;
        } else if (materials[i] != Definitions::Material::FLUID &&
                   materials[(i + 1) % 4] == Definitions::Material::FLUID)
          surfaceField[neighborsId[static_cast<size_t>(sides[(i + 1) % 4])]] =
              1;
      }
    }
  }
}

template <>
inline void SimulationSteps::makeGraded(StructureInterface<2> *structure,
                                        size_t cellMaterialFieldId) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  std::queue<size_t> fluidCells;
  cellGraph->iterateCells([&](size_t i) {
    if (cellMaterialField[i] == Definitions::Material::FLUID)
      fluidCells.push(i);
  });
  while (!fluidCells.empty()) {
    size_t curCellId = fluidCells.front();
    fluidCells.pop();
    auto neighbors = cellGraph->neighbors(curCellId);
    for (auto n : neighbors)
      if (n.id > 0 && cellGraph->node(n.id).isValid()) {
        if (Definitions::Material::FLUID == cellMaterialField[n.id]) {
          if (cellGraph->node(n.id).level() <
              cellGraph->node(curCellId).level() - 2) {
            auto children = cellGraph->refine(n.id);
            for (auto c : children) {
              fluidCells.push(c);
              cellMaterialField[c] = Definitions::Material::FLUID;
            }
          }
        }
      }
  }
}

template <>
inline void SimulationSteps::makeGraded(StructureInterface<3> *structure,
                                        size_t cellMaterialFieldId) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  std::queue<size_t> fluidCells;
  cellGraph->iterateCells([&](size_t i) {
    if (cellMaterialField[i] == Definitions::Material::FLUID)
      fluidCells.push(i);
  });
  while (!fluidCells.empty()) {
    size_t curCellId = fluidCells.front();
    fluidCells.pop();
    auto neighbors = cellGraph->neighbors(curCellId);
    for (auto n : neighbors)
      if (n.id > 0 && cellGraph->node(n.id).isValid()) {
        if (Definitions::Material::FLUID == cellMaterialField[n.id]) {
          if (cellGraph->node(n.id).level() <
              cellGraph->node(curCellId).level() - 1) {
            auto children = cellGraph->refine(n.id);
            for (auto c : children) {
              fluidCells.push(c);
              cellMaterialField[c] = Definitions::Material::FLUID;
            }
          }
        }
      }
  }
}

template <int D>
void SimulationSteps::updateTree(ParticleSystem<double, D> *particles,
                                 StructureInterface<D> *structure,
                                 size_t surfaceMaskFieldId,
                                 size_t surfaceTDistanceFieldId,
                                 size_t cellMaterialFieldId, size_t minLevel,
                                 size_t surfaceLevel, bool graded) {
  std::cerr << "updating surface cells\n";
  auto surfaceCells =
      updateSurfaceCells(particles, structure, surfaceMaskFieldId,
                         cellMaterialFieldId, surfaceLevel);
  std::cerr << "coarsing\n";
  coarseGraph(structure, surfaceMaskFieldId, cellMaterialFieldId,
              [&](Definitions::Material nodeMaterialType, size_t level) {
                return !(nodeMaterialType == Definitions::Material::FLUID &&
                         level <= minLevel);
              });
  std::cerr << "refining near surface\n";
  try {
    auto newSurface = refineNearSurface(structure, surfaceMaskFieldId,
                                        surfaceTDistanceFieldId,
                                        cellMaterialFieldId, surfaceCells);
  } catch (const char *e) {
    std::cerr << e << std::endl;
    THROW(false, "SimulationSteps::updateTree");
  }
  if (graded)
    makeGraded(structure, cellMaterialFieldId);
  //  fixSurface(structure, surfaceMaskFieldId, cellMaterialFieldId,
  //  newSurface);
}

template <int D>
void SimulationSteps::updateTreeCoarseBottom(
    ParticleSystem<double, D> *particles, StructureInterface<D> *structure,
    size_t surfaceMaskFieldId, size_t surfaceTDistanceFieldId,
    size_t cellMaterialFieldId, size_t minLevel, size_t surfaceLevel,
    bool graded) {
  std::cerr << "updating surface cells\n";
  auto surfaceCells =
      updateSurfaceCellsCoarseBottom(particles, structure, surfaceMaskFieldId,
                                     cellMaterialFieldId, surfaceLevel);
  std::cerr << "coarsing\n";
  coarseGraph(structure, surfaceMaskFieldId, cellMaterialFieldId,
              [&](Definitions::Material nodeMaterialType, size_t level) {
                return !(nodeMaterialType == Definitions::Material::FLUID &&
                         level <= minLevel);
              });
  try {
    std::cerr << "refining near surface\n";
    auto newSurface = refineNearSurface(structure, surfaceMaskFieldId,
                                        surfaceTDistanceFieldId,
                                        cellMaterialFieldId, surfaceCells);
    std::cerr << "refining near solid\n";
    refineNearSolid(structure, cellMaterialFieldId, surfaceLevel);
  } catch (const char *e) {
    std::cerr << e << std::endl;
    THROW(false, "SimulationSteps::updateTree not refining");
  }
  if (graded) {
    std::cerr << "Making graded\n";
    makeGraded(structure, cellMaterialFieldId);
  }
  //  fixSurface(structure, surfaceMaskFieldId, cellMaterialFieldId,
  //  newSurface);
}

template <int D>
void SimulationSteps::markSolids(
    StructureInterface<D> *structure, size_t cellMaterialFieldId,
    const std::vector<std::shared_ptr<SolidInterface<D>>> &solids) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      return;
    auto center = structure->cellCenterPosition(cellId);
    for (auto &solid : solids)
      if (solid->contains(center)) {
        cellMaterialField[cellId] = Definitions::Material::SOLID;
        return;
      }
    cellMaterialField[cellId] = Definitions::Material::AIR;
  });
}

template <int D>
void SimulationSteps::markFaces(StructureInterface<D> *structure,
                                size_t cellMaterialFieldId,
                                size_t faceMaterialFieldId,
                                size_t faceBoundaryFieldId,
                                Definitions::Boundary *domainBoundaries) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);
  auto &faceBoundaryField =
      *structure->template field<Definitions::Boundary>(faceBoundaryFieldId);
  structure->iterateFaces([&](size_t faceId) {
    // First, lets retrieve the cells it face connects
    auto cells = structure->faceCells(faceId);
    if (cells.size() == 1)
      cells.push_back({-1, Definitions::oppositeSide(cells[0].side)});
    THROW(cells.size() == 2u,
          "SimulationSteps::markFaces wrong number of face connecting cells.");
    // The face boundary type is assumed none
    Definitions::Boundary faceBoundaryType = Definitions::Boundary::NONE;
    // First we assign domain boundary information to the face in case it is
    // at one of the extreme sides of the simulation domain (faces that have
    // only one valid cell)
    for (size_t n = 0; n < 2; n++)
      if (cells[n].id < 0)
        switch (cells[n].side) {
        case Definitions::Side::LEFT:
          faceBoundaryType = domainBoundaries[0];
          break;
        case Definitions::Side::RIGHT:
          faceBoundaryType = domainBoundaries[1];
          break;
        case Definitions::Side::BOTTOM:
          faceBoundaryType = domainBoundaries[2];
          break;
        case Definitions::Side::TOP:
          faceBoundaryType = domainBoundaries[3];
          break;
        case Definitions::Side::BACK:
          faceBoundaryType = domainBoundaries[4];
          break;
        case Definitions::Side::FRONT:
          faceBoundaryType = domainBoundaries[5];
          break;
        default:
          NOT_IMPLEMENTED();
        }
    // We assume the domain boundaries are solid walls, that is why we assign
    // SOLID to cells that are outside the domain, existent cells receive
    // their stored material
    std::vector<Definitions::Material> cellType(2,
                                                Definitions::Material::SOLID);
    for (size_t i = 0; i < 2; i++)
      if (cells[i].id >= 0)
        cellType[i] = cellMaterialField[cells[i].id];
    // Now, the face receives its type based on the types of the adjacent
    // cells DIRICHLET type is assigned to fluid-air interfaces
    if (cells[0].id >= 0 && cells[1].id >= 0) {
      if ((cellType[0] == Definitions::Material::FLUID &&
           cellType[1] == Definitions::Material::AIR) ||
          (cellType[1] == Definitions::Material::FLUID &&
           cellType[0] == Definitions::Material::AIR))
        faceBoundaryType = Definitions::Boundary::DIRICHLET;
      if (cellType[0] == Definitions::Material::SOLID ||
          cellType[1] == Definitions::Material::SOLID)
        faceBoundaryType = Definitions::Boundary::NEUMANN;
    }
    faceBoundaryField[faceId] = faceBoundaryType;
    // Now that the boundary type has been determined, we can compute the face
    // material
    auto faceMaterialType = Definitions::Material::AIR;
    if (cells[0].id < 0)
      faceMaterialType = (faceBoundaryType == Definitions::Boundary::DIRICHLET)
                             ? cellType[1]
                             : Definitions::Material::SOLID;
    else if (cells[1].id < 0)
      faceMaterialType = (faceBoundaryType == Definitions::Boundary::DIRICHLET)
                             ? cellType[0]
                             : Definitions::Material::SOLID;
    else {
      // the precedence order is SOLID -> FLUID -> AIR
      if (cellType[0] == Definitions::Material::SOLID ||
          cellType[1] == Definitions::Material::SOLID)
        faceMaterialType = Definitions::Material::SOLID;
      else if (cellType[0] == Definitions::Material::FLUID ||
               cellType[1] == Definitions::Material::FLUID)
        faceMaterialType = Definitions::Material::FLUID;
      else
        faceMaterialType = Definitions::Material::AIR;
    }
    faceMaterialField[faceId] = faceMaterialType;
  });
}

using GridFunction = std::function<void(size_t)>;

template <typename Grid>
inline static void iterateGridCells(Grid *grid, const GridFunction &f) {
#ifdef _USE_OPENMP
  grid->iterateCells_par(f);
#else
  grid->iterateCells(f);
#endif // _USE_OPENMP
}

template <typename Grid>
inline static void iterateGridFaces(Grid *grid, const GridFunction &f) {
#ifdef _USE_OPENMP
  grid->iterateFaces_par(f);
#else
  grid->iterateFaces(f);
#endif // _USE_OPENMP
}

template <typename Grid>
inline static void iterateGridVertices(Grid *grid, const GridFunction &f) {
#ifdef _USE_OPENMP
  grid->iterateVertices_par(f);
#else
  grid->iterateVertices(f);
#endif // _USE_OPENMP
}

template <int D>
void SimulationSteps::transferFromParticlesToGrid(
    const ParticleSystem<double, D> *particles, size_t propertyId,
    StructureInterface<D> *grid, size_t fieldId,
    const std::function<bool(size_t)> &useElement,
    Interpolant<Point<double, D>> *_interpolator) {
  auto &field = *grid->template field<double>(fieldId);
  QuinticKernel<Point<double, D>> kernel;
  std::vector<size_t> useableFaces;
  switch (grid->fieldLocation(fieldId)) {
  case Definitions::MeshLocation::CELL_CENTER:
    iterateGridCells(grid, [&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    iterateGridVertices(grid, [&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    const_cast<ParticleSystem<double, D> *>(particles)
        ->update(); // invoked here to avoid race conditions with OpenMP
#pragma omp parallel
    {
      QuinticKernel<Point<double, D>> kernel;
      RBFInterpolant<Point<double, D>> interpolator =
          RBFInterpolant<Point<double, D>>(kernel);
      iterateGridFaces(grid, [&](size_t id) {
        field[id] = 0.;
        if ((grid->fieldLocation(fieldId) ==
                 Definitions::MeshLocation::HORIZONTAL_FACE_CENTER &&
             grid->faceOrientation(id) !=
                 Definitions::Orientation::HORIZONTAL) ||
            (grid->fieldLocation(fieldId) ==
                 Definitions::MeshLocation::VERTICAL_FACE_CENTER &&
             grid->faceOrientation(id) != Definitions::Orientation::VERTICAL) ||
            (grid->fieldLocation(fieldId) ==
                 Definitions::MeshLocation::DEPTH_FACE_CENTER &&
             grid->faceOrientation(id) != Definitions::Orientation::DEPTH) ||
            !useElement(id))
          return;
        std::vector<Point<double, D>> points;
        std::vector<double> values;
        auto neighbors = grid->faceCells(id);
        // Change search region if neighbors are from different size
        std::vector<BBox<double, D>> searchRegions;
        auto cellLevel = [&](size_t cellId) -> size_t {
          if (D == 2)
            return dynamic_cast<CellGraph2 *>(grid)->node(cellId).level();
          return dynamic_cast<CellGraph3 *>(grid)->node(cellId).level();
        };
        std::vector<size_t> levels;
        if (neighbors.size() == 1) {
          searchRegions.emplace_back(grid->cellRegion(neighbors[0].id));
          levels.emplace_back(cellLevel(neighbors[0].id));
        } else {
          BBox<double, D> region[2];
          region[0] = grid->cellRegion(neighbors[0].id);
          region[1] = grid->cellRegion(neighbors[1].id);
          levels.emplace_back(cellLevel(neighbors[0].id));
          searchRegions.emplace_back(region[0]);
          levels.emplace_back(cellLevel(neighbors[1].id));
          searchRegions.emplace_back(region[1]);
        }
        std::set<size_t> pointIds;
        for (size_t r = 0; r < searchRegions.size(); r++) {
          particles->iterateParticlesZ(
              searchRegions[r], levels[r],
              [&](size_t pid, Point<double, D> p) { pointIds.insert(pid); });
        }
        for (auto pid : pointIds) {
          values.emplace_back(particles->getScalarProperty(propertyId, pid));
          points.emplace_back((*particles)[pid]);
        }
        try {
          field[id] = interpolator.interpolateAt(grid->faceCenterPosition(id),
                                                 points, values);
        } catch (std::string &e) {
          std::cerr << e << std::endl;
          std::cerr
              << "SimulationSteps::transferFromParticlesToGrid: attempting "
                 "with bigger neighborhood at face "
              << grid->faceCenterPosition(id) << std::endl;
          for (auto v : values)
            std::cerr << v << " ";
          std::cerr << std::endl;
          std::set<int> neighborSet;
          for (auto cell : neighbors)
            for (auto n : grid->cellNeighbors(cell.id))
              neighborSet.insert(n.id);
          for (auto n : neighborSet)
            if (n >= 0)
              particles->iterateParticles(
                  grid->cellRegion(n), [&](size_t pid, Point<double, D> p) {
                    values.emplace_back(
                        particles->getScalarProperty(propertyId, pid));
                    points.emplace_back(p);
                  });

          try {
            field[id] = interpolator.interpolateAt(grid->faceCenterPosition(id),
                                                   points, values);
          } catch (std::string &e) {
            std::cerr << grid->faceCenterPosition(id) << "   <<< " << id
                      << std::endl;
            THROW(false, "SimulationSteps::transferFromParticlesToGrid: second "
                         "attempt failed");
          }
        }
      });
    }
    break;
  default:
    THROW(false, "SimulationSteps::transferFromParticlesToGrid invalid field "
                 "location");
  }
}

template <>
inline void SimulationSteps::transferFromParticlesToGridHandlingCellLevels(
    const ParticleSystem<double, 2> *particles, size_t propertyId,
    StructureInterface<2> *grid, size_t cellMaterialFieldId, size_t fieldId,
    const std::function<bool(size_t)> &useElement,
    Interpolant<Point<double, 2>> *interpolator) {
  auto &cellMaterialField =
      *grid->template field<Definitions::Material>(cellMaterialFieldId);
  int D = 2;
  auto &field = *grid->template field<double>(fieldId);
  switch (grid->fieldLocation(fieldId)) {
  case Definitions::MeshLocation::CELL_CENTER:
    grid->iterateCells([&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    grid->iterateVertices([&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    grid->iterateFaces([&](size_t id) {
      field[id] = 0.;
      if ((grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::HORIZONTAL_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::HORIZONTAL) ||
          (grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::VERTICAL_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::VERTICAL) ||
          (grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::DEPTH_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::DEPTH) ||
          !useElement(id))
        return;
      auto faceNeighbors = grid->faceCells(id);
      std::vector<Definitions::Neighbor> neighbors;
      for (auto fn : faceNeighbors)
        if (fn.id >= 0 &&
            cellMaterialField[fn.id] == Definitions::Material::FLUID)
          neighbors.emplace_back(fn);
      // Change search region if neighbors are from different size
      std::vector<BBox<double, 2>> searchRegions;
      std::vector<double> weights;
      std::vector<double> value(1, 0.);
      if (neighbors.size() == 1) {
        searchRegions.emplace_back(grid->cellRegion(neighbors[0].id));
        weights.emplace_back(1);
      } else {
        BBox<double, 2> region[2];
        region[0] = grid->cellRegion(neighbors[0].id);
        region[1] = grid->cellRegion(neighbors[1].id);
        double ratio = region[0].size()[0] / region[1].size()[0];
        if (ratio == 1.) {
          searchRegions.emplace_back(
              BBox<double, 2>::combine(region[0], region[1]));
          weights.emplace_back(1);
        } else {
          value.emplace_back(0.);
          searchRegions.emplace_back(region[0]);
          searchRegions.emplace_back(region[1]);
          double volume[2] = {1., 1.};
          for (int d = 0; d < D; d++) {
            volume[0] *= region[0].size()[d];
            volume[1] *= region[1].size()[d];
          }
          volume[0] = 1 / volume[0];
          volume[1] = 1 / volume[1];
          weights.emplace_back(volume[0] / (volume[0] + volume[1]));
          weights.emplace_back(volume[1] / (volume[0] + volume[1]));
        }
      }
      for (size_t i = 0; i < searchRegions.size(); i++) {
        std::vector<Point<double, 2>> points;
        std::vector<double> values;
        particles->iterateParticles(searchRegions[i], [&](size_t pid,
                                                          Point<double, 2> p) {
          values.emplace_back(particles->getScalarProperty(propertyId, pid));
          points.emplace_back(p);
        });
        try {
          value[i] = interpolator->interpolateAt(grid->faceCenterPosition(id),
                                                 points, values);
        } catch (const std::string &) {
          points.clear();
          values.clear();
          particles->iterateClosestParticles(
              grid->faceCenterPosition(id), (D == 2) ? 4 : 8,
              [&](size_t pid, Point<double, 2> p) {
                values.emplace_back(
                    particles->getScalarProperty(propertyId, pid));
                points.emplace_back(p);
              });

          value[i] = interpolator->interpolateAt(grid->faceCenterPosition(id),
                                                 points, values);
        }
      }
      field[id] = 0.;
      for (size_t i = 0; i < searchRegions.size(); i++)
        field[id] += weights[i] * value[i];
    });
    break;
  default:
    THROW(false, "SimulationSteps::transferFromParticlesToGrid invalid field "
                 "location");
  }
}

template <>
inline void SimulationSteps::transferFromParticlesToGridHandlingCellLevels(
    const ParticleSystem<double, 3> *particles, size_t propertyId,
    StructureInterface<3> *grid, size_t cellMaterialFieldId, size_t fieldId,
    const std::function<bool(size_t)> &useElement,
    Interpolant<Point<double, 3>> *interpolator) {
  auto &cellMaterialField =
      *grid->template field<Definitions::Material>(cellMaterialFieldId);
  int D = 3;
  auto &field = *grid->template field<double>(fieldId);
  switch (grid->fieldLocation(fieldId)) {
  case Definitions::MeshLocation::CELL_CENTER:
    grid->iterateCells([&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    grid->iterateVertices([&](size_t id) {
      field[id] = 0.;
      NOT_IMPLEMENTED();
    });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    grid->iterateFaces([&](size_t id) {
      field[id] = 0.;
      if ((grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::HORIZONTAL_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::HORIZONTAL) ||
          (grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::VERTICAL_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::VERTICAL) ||
          (grid->fieldLocation(fieldId) ==
               Definitions::MeshLocation::DEPTH_FACE_CENTER &&
           grid->faceOrientation(id) != Definitions::Orientation::DEPTH) ||
          !useElement(id))
        return;
      auto faceNeighbors = grid->faceCells(id);
      std::vector<Definitions::Neighbor> neighbors;
      for (auto fn : faceNeighbors)
        if (fn.id >= 0 &&
            cellMaterialField[fn.id] == Definitions::Material::FLUID)
          neighbors.emplace_back(fn);
      // Change search region if neighbors are from different size
      std::vector<BBox<double, 3>> searchRegions;
      std::vector<double> weights;
      std::vector<double> value(1, 0.);
      if (neighbors.size() == 1) {
        searchRegions.emplace_back(grid->cellRegion(neighbors[0].id));
        weights.emplace_back(1);
      } else {
        BBox<double, 3> region[2];
        region[0] = grid->cellRegion(neighbors[0].id);
        region[1] = grid->cellRegion(neighbors[1].id);
        double ratio = region[0].size()[0] / region[1].size()[0];
        if (ratio == 1.) {
          searchRegions.emplace_back(
              BBox<double, 3>::combine(region[0], region[1]));
          weights.emplace_back(1);
        } else {
          value.emplace_back(0.);
          searchRegions.emplace_back(region[0]);
          searchRegions.emplace_back(region[1]);
          double volume[2] = {1., 1.};
          for (int d = 0; d < D; d++) {
            volume[0] *= region[0].size()[d];
            volume[1] *= region[1].size()[d];
          }
          volume[0] = 1 / volume[0];
          volume[1] = 1 / volume[1];
          weights.emplace_back(volume[0] / (volume[0] + volume[1]));
          weights.emplace_back(volume[1] / (volume[0] + volume[1]));
        }
      }
      for (size_t i = 0; i < searchRegions.size(); i++) {
        std::vector<Point<double, 3>> points;
        std::vector<double> values;
        particles->iterateParticles(searchRegions[i], [&](size_t pid,
                                                          Point<double, 3> p) {
          values.emplace_back(particles->getScalarProperty(propertyId, pid));
          points.emplace_back(p);
        });
        try {
          value[i] = interpolator->interpolateAt(grid->faceCenterPosition(id),
                                                 points, values);
        } catch (const std::string &) {
          points.clear();
          values.clear();
          particles->iterateClosestParticles(
              grid->faceCenterPosition(id), (D == 2) ? 4 : 8,
              [&](size_t pid, Point<double, 3> p) {
                values.emplace_back(
                    particles->getScalarProperty(propertyId, pid));
                points.emplace_back(p);
              });

          value[i] = interpolator->interpolateAt(grid->faceCenterPosition(id),
                                                 points, values);
        }
      }
      field[id] = 0.;
      for (size_t i = 0; i < searchRegions.size(); i++)
        field[id] += weights[i] * value[i];
    });
    break;
  default:
    THROW(false, "SimulationSteps::transferFromParticlesToGrid invalid field "
                 "location");
  }
}

template <int D>
void SimulationSteps::transferFromParticlesToGrid(
    const ParticleSystem<double, D> *particles, size_t propertyId, size_t n,
    StructureInterface<D> *grid, size_t fieldId,
    const std::function<bool(size_t)> &useElement) {
  auto &field = *grid->template field<double>(fieldId);
  switch (grid->fieldLocation(fieldId)) {
  case Definitions::MeshLocation::CELL_CENTER:
    grid->iterateCells([&](size_t id) {
      field[id] = 0.;
      if (useElement(id))
        field[id] = particles->sampleScalarProperty(
            propertyId, grid->cellCenterPosition(id), n);
    });
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    grid->iterateVertices([&](size_t id) {
      field[id] = 0.;
      if (useElement(id))
        field[id] = particles->sampleScalarProperty(
            propertyId, grid->vertexPosition(id), n);
    });
    break;
  case Definitions::MeshLocation::FACE_CENTER:
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    grid->iterateFaces([&](size_t id) {
      field[id] = 0.;
      if (useElement(id))
        field[id] = particles->sampleScalarProperty(
            propertyId, grid->faceCenterPosition(id), n);
    });
    break;
  default:
    THROW(false, "SimulationSteps::transferFromParticlesToGrid invalid field "
                 "location");
  }
}

template <int D>
void SimulationSteps::applyExternalForce(
    StructureInterface<D> *structure, size_t fieldId, double acceleration,
    double dt, const std::function<bool(size_t)> &useElement) {
  auto &field = *structure->template field<double>(fieldId);
  std::function<void(size_t)> callback = [&](size_t i) {
    if (useElement && !useElement(i))
      return;
    field[i] += acceleration * dt;
  };
  switch (structure->fieldLocation(fieldId)) {
  case Definitions::MeshLocation::CELL_CENTER:
    structure->iterateCells(callback);
    break;
  case Definitions::MeshLocation::VERTEX_CENTER:
    structure->iterateVertices(callback);
    break;
  case Definitions::MeshLocation::FACE_CENTER:
  case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
  case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
  case Definitions::MeshLocation::DEPTH_FACE_CENTER:
    structure->iterateFaces([&](size_t id) {
      if ((structure->fieldLocation(fieldId) ==
               Definitions::MeshLocation::HORIZONTAL_FACE_CENTER &&
           structure->faceOrientation(id) !=
               Definitions::Orientation::HORIZONTAL) ||
          (structure->fieldLocation(fieldId) ==
               Definitions::MeshLocation::VERTICAL_FACE_CENTER &&
           structure->faceOrientation(id) !=
               Definitions::Orientation::VERTICAL) ||
          (structure->fieldLocation(fieldId) ==
               Definitions::MeshLocation::DEPTH_FACE_CENTER &&
           structure->faceOrientation(id) != Definitions::Orientation::DEPTH))
        return;
      if (useElement && !useElement(id))
        return;
      field[id] += acceleration * dt;
    });
    break;
  default:
    THROW(false, "SimulationSteps::transferFromParticlesToGrid invalid field "
                 "location");
  }
}

template <int D>
void SimulationSteps::applyExternalForce(ParticleSystem<double, D> *particles,
                                         size_t propertyId, double acceleration,
                                         double dt) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    particles->setScalarProperty(propertyId, id,
                                 particles->getScalarProperty(propertyId, id) +
                                     acceleration * dt);
  });
}

template <int D>
inline std::vector<double>
SimulationSteps::computePressureRingMatrixStencilWeights(
    const StructureInterface<D> *structure, size_t cellId,
    DifferentialRBF<D> *rbf, const std::vector<size_t> &neighbors) {
  auto cellCenter = structure->cellCenterPosition(cellId);
  // First of all, we need to construct the set of stencil points
  std::vector<Point<double, D>> points;
  points.emplace_back(cellCenter);
  for (auto nid : neighbors) {
    points.emplace_back(structure->cellCenterPosition(nid));
  }
  // Now the weights
  // if we have a regular stencil we can use direct computed weights
  // if (points.size() == 5 && sameLevelNeighbors == 4)
  //   return rbf->laplacianWeights(distance(points[0], points[1]));
  return rbf->laplacianWeights(points);
}

template <>
inline std::vector<double> SimulationSteps::computePressureMatrixStencilWeights(
    const StructureInterface<2> *structure, size_t cellId,
    size_t cellMaterialFieldId, size_t cellParticleCenterFieldId,
    DifferentialRBF<2> *rbf, bool adaptiveRBFRadius) {
  auto *cellGraph = dynamic_cast<const CellGraph2 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *cellGraph->field<Point<double, 2>>(cellParticleCenterFieldId);
  auto cellRegion = cellGraph->cellRegion(cellId);
  auto inf = Point<double, 2>(INFINITY);
  THROW(cellParticleCenterField[cellId] != inf,
        "SimulationSteps::computePressureMatrixStencilWeights invalid center "
        "particle");
  auto cellCenter = cellParticleCenterField[cellId];
  if (cellMaterialField[cellId] != Definitions::Material::FLUID)
    cellCenter = cellGraph->cellCenterPosition(cellId);
  auto neighborCells = cellGraph->cellNeighbors(cellId);
  // First of all, we need to construct the set of stencil points
  std::vector<Point2d> points;
  points.emplace_back(cellCenter);
  size_t sameLevelNeighbors = 0;
  size_t centerLevel = cellGraph->node(cellId).level();
  for (auto n : neighborCells)
    if (n.id >= 0) {
      if (cellGraph->node(n.id).level() == centerLevel)
        sameLevelNeighbors++;
      THROW(cellParticleCenterField[n.id] != inf,
            "SimulationSteps::computePressureMatrixStencilWeights invalid "
            "neighbor center "
            "particle");
      if (cellMaterialField[n.id] == Definitions::Material::FLUID)
        points.emplace_back(cellParticleCenterField[n.id]);
      else
        points.emplace_back(cellGraph->cellCenterPosition(n.id));
    } else {
      sameLevelNeighbors++;
      switch (n.side) {
      case Definitions::Side::LEFT:
        points.emplace_back(Point2d(
            2. * cellRegion.lower().x() - cellCenter.x(), cellCenter.y()));
        break;
      case Definitions::Side::RIGHT:
        points.emplace_back(Point2d(
            2. * cellRegion.upper().x() - cellCenter.x(), cellCenter.y()));
        break;
      case Definitions::Side::TOP:
        points.emplace_back(Point2d(
            cellCenter.x(), 2. * cellRegion.upper().y() - cellCenter.y()));
        break;
      case Definitions::Side::BOTTOM:
        points.emplace_back(Point2d(
            cellCenter.x(), 2. * cellRegion.lower().y() - cellCenter.y()));
        break;
      default:
        THROW(false, "SimulationSteps::computePressureMatrixWeights invalid "
                     "neighbor side");
      }
    }
    // Now the weights
#ifdef FUROO_PROFILE
  if (points.size() == 5 && sameLevelNeighbors == 4)
    CodeProfiler::count("#RegularStencils");
  else
    CodeProfiler::count(concat("#Stencil", points.size()));
#endif
  return rbf->laplacianWeights(points);
}

template <>
inline std::vector<double> SimulationSteps::computePressureMatrixStencilWeights(
    const StructureInterface<3> *structure, size_t cellId,
    size_t cellMaterialFieldId, size_t cellParticleCenterFieldId,
    DifferentialRBF<3> *rbf, bool adaptiveRBFRadius) {
  auto *cellGraph = dynamic_cast<const CellGraph3 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *cellGraph->field<Point<double, 3>>(cellParticleCenterFieldId);
  auto cellRegion = cellGraph->cellRegion(cellId);
  THROW(cellParticleCenterField[cellId] != Point3d(INFINITY),
        "SimulationSteps::computePressureMatrixStencilWeights invalid center "
        "particle");
  auto cellCenter = cellParticleCenterField[cellId];
  auto neighborCells = cellGraph->cellNeighbors(cellId);
  // First of all, we need to construct the set of stencil points
  std::vector<Point3d> points;
  points.emplace_back(cellCenter);
  size_t sameLevelNeighbors = 0;
  size_t centerLevel = cellGraph->node(cellId).level();
  for (auto n : neighborCells)
    if (n.id >= 0) {
      if (cellGraph->node(n.id).level() == centerLevel)
        sameLevelNeighbors++;
      THROW(cellParticleCenterField[n.id] != Point3d(INFINITY),
            "SimulationSteps::computePressureMatrixStencilWeights invalid "
            "neighbor center "
            "particle");
      points.emplace_back(cellParticleCenterField[n.id]);
    } else {
      sameLevelNeighbors++;
      switch (n.side) {
      case Definitions::Side::LEFT:
        points.emplace_back(
            Point3d(2. * cellRegion.lower().x() - cellCenter.x(),
                    cellCenter.y(), cellCenter.z()));
        break;
      case Definitions::Side::RIGHT:
        points.emplace_back(
            Point3d(2. * cellRegion.upper().x() - cellCenter.x(),
                    cellCenter.y(), cellCenter.z()));
        break;
      case Definitions::Side::BOTTOM:
        points.emplace_back(Point3d(
            cellCenter.x(), 2. * cellRegion.lower().y() - cellCenter.y(),
            cellCenter.z()));
        break;
      case Definitions::Side::TOP:
        points.emplace_back(Point3d(
            cellCenter.x(), 2. * cellRegion.upper().y() - cellCenter.y(),
            cellCenter.z()));
        break;
      case Definitions::Side::BACK:
        points.emplace_back(
            Point3d(cellCenter.x(), cellCenter.y(),
                    2. * cellRegion.lower().z() - cellCenter.z()));
        break;
      case Definitions::Side::FRONT:
        points.emplace_back(
            Point3d(cellCenter.x(), cellCenter.y(),
                    2. * cellRegion.upper().z() - cellCenter.z()));
        break;
      default:
        THROW(false,
              "AdaptiveCenteredPic2::solvePressure invalid neighbor side");
      }
    }
    // Now the weights
    // if we have a regular stencil we can use direct computed weights
    //  if (points.size() == 7 && sameLevelNeighbors == 6)
    //    return rbf.laplacianWeights(distance(points[0], points[1]));

#ifdef FUROO_PROFILE
  if (points.size() == 7 && sameLevelNeighbors == 6)
    CodeProfiler::count("#RegularStencils");
  else
    CodeProfiler::count(concat("#Stencil", points.size()));
#endif
  return rbf->laplacianWeights(points);
}

template <>
inline std::vector<double> SimulationSteps::computePressureMatrixStencilWeights(
    const StructureInterface<2> *structure, size_t cellId,
    size_t cellMaterialFieldId, DifferentialRBF<2> *rbf,
    bool adaptiveRBFRadius) {
  auto *cellGraph = dynamic_cast<const CellGraph2 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  auto cellRegionSize = cellGraph->cellRegion(cellId).size();
  auto cellCenter = cellGraph->cellCenterPosition(cellId);
  auto neighborCells = cellGraph->cellNeighbors(cellId);
  // First of all, we need to construct the set of stencil points
  std::vector<Point2d> points;
  points.emplace_back(cellCenter);
  size_t sameLevelNeighbors = 0;
  size_t centerLevel = cellGraph->node(cellId).level();
  for (auto n : neighborCells)
    if (n.id >= 0 || cellMaterialField[n.id] != Definitions::Material::SOLID) {
      if (cellGraph->node(n.id).level() == centerLevel)
        sameLevelNeighbors++;
      points.emplace_back(cellGraph->cellCenterPosition(n.id));
    } else {
      sameLevelNeighbors++;
      switch (n.side) {
      case Definitions::Side::LEFT:
        points.emplace_back(
            Point2d(cellCenter.x() - cellRegionSize.x(), cellCenter.y()));
        break;
      case Definitions::Side::RIGHT:
        points.emplace_back(
            Point2d(cellCenter.x() + cellRegionSize.x(), cellCenter.y()));
        break;
      case Definitions::Side::TOP:
        points.emplace_back(
            Point2d(cellCenter.x(), cellCenter.y() + cellRegionSize.y()));
        break;
      case Definitions::Side::BOTTOM:
        points.emplace_back(
            Point2d(cellCenter.x(), cellCenter.y() - cellRegionSize.y()));
        break;
      default:
        THROW(false,
              "AdaptiveCenteredPic2::solvePressure invalid neighbor side");
      }
    }
  // Now the weights
  // if we have a regular stencil we can use direct computed weights
  // if (points.size() == 5 && sameLevelNeighbors == 4)
  //   return rbf->laplacianWeights(distance(points[0], points[1]));
  return rbf->laplacianWeights(points);
}

template <>
inline std::vector<double> SimulationSteps::computePressureMatrixStencilWeights(
    const StructureInterface<3> *structure, size_t cellId,
    size_t cellMaterialFieldId, DifferentialRBF<3> *rbf,
    bool adaptiveRBFRadius) {
  auto *cellGraph = dynamic_cast<const CellGraph3 *>(structure);
  auto &cellMaterialField =
      *cellGraph->field<Definitions::Material>(cellMaterialFieldId);
  auto cellRegionSize = cellGraph->cellRegion(cellId).size();
  auto cellCenter = cellGraph->cellCenterPosition(cellId);
  auto neighborCells = cellGraph->cellNeighbors(cellId);
  // First of all, we need to construct the set of stencil points
  std::vector<Point3d> points;
  points.emplace_back(cellCenter);
  size_t sameLevelNeighbors = 0;
  size_t centerLevel = cellGraph->node(cellId).level();
  for (auto n : neighborCells)
    if (n.id >= 0 || cellMaterialField[n.id] != Definitions::Material::SOLID) {
      if (cellGraph->node(n.id).level() == centerLevel)
        sameLevelNeighbors++;
      points.emplace_back(cellGraph->cellCenterPosition(n.id));
    } else {
      sameLevelNeighbors++;
      switch (n.side) {
      case Definitions::Side::LEFT:
        points.emplace_back(Point3d(cellCenter.x() - cellRegionSize.x(),
                                    cellCenter.y(), cellCenter.z()));
        break;
      case Definitions::Side::RIGHT:
        points.emplace_back(Point3d(cellCenter.x() + cellRegionSize.x(),
                                    cellCenter.y(), cellCenter.z()));
        break;
      case Definitions::Side::BOTTOM:
        points.emplace_back(Point3d(cellCenter.x(),
                                    cellCenter.y() - cellRegionSize.y(),
                                    cellCenter.z()));
        break;
      case Definitions::Side::TOP:
        points.emplace_back(Point3d(cellCenter.x(),
                                    cellCenter.y() + cellRegionSize.y(),
                                    cellCenter.z()));
        break;
      case Definitions::Side::BACK:
        points.emplace_back(Point3d(cellCenter.x(), cellCenter.y(),
                                    cellCenter.z() - cellRegionSize.z()));
        break;
      case Definitions::Side::FRONT:
        points.emplace_back(Point3d(cellCenter.x(), cellCenter.y(),
                                    cellCenter.z() + cellRegionSize.z()));
        break;
      default:
        THROW(false,
              "AdaptiveCenteredPic2::solvePressure invalid neighbor side");
      }
    }
  // Now the weights
  // if we have a regular stencil we can use direct computed weights
  //  if (points.size() == 7 && sameLevelNeighbors == 6)
  //    return rbf.laplacianWeights(distance(points[0], points[1]));
  return rbf->laplacianWeights(points);
}

template <int D>
void SimulationSteps::transferCellFieldToFaceField(
    StructureInterface<D> *cellGraph, size_t cellFieldId, size_t faceFieldId,
    const std::function<bool(size_t)> &useCell,
    const std::function<bool(size_t)> &useFace) {
  auto &cellField = *cellGraph->template field<double>(cellFieldId);
  auto &faceField = *cellGraph->template field<double>(faceFieldId);
  cellGraph->iterateFaces([&](size_t faceId) {
    if (!useFace(faceId))
      return;
    auto cells = cellGraph->faceCells(faceId);
    std::vector<double> values;
    std::vector<Point<double, D>> points;
    std::set<size_t> pointSet;
    for (auto cell : cells) {
      auto cellNeighbors = cellGraph->cellNeighbors(cell.id);
      for (auto neighbor : cellNeighbors)
        if (neighbor.id >= 0 && useCell(neighbor.id))
          pointSet.insert(neighbor.id);
      for (auto pointId : pointSet) {
        points.emplace_back(cellGraph->cellCenterPosition(pointId));
        values.emplace_back(cellField[pointId]);
      }
    }
    MLS<D> interpolant;
    double value = interpolant.interpolateAt(
        cellGraph->faceCenterPosition(faceId), points, values);
    faceField[faceId] = value;
  });
}
template <int D>
void SimulationSteps::computePressureGradient(
    StructureInterface<D> *cellGraph, size_t pressureCellFieldId,
    size_t pressureFaceFieldId,
    const std::vector<size_t> &pressureGradientFaceFieldIds,
    size_t faceMaterialFieldId, DifferentialRBF<D> *rbf) {
  auto &faceMaterialField =
      *cellGraph->template field<Definitions::Material>(faceMaterialFieldId);
  auto &pressureField = *cellGraph->template field<double>(pressureCellFieldId);
  auto &facePressureField =
      *cellGraph->template field<double>(pressureFaceFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        cellGraph->template field<double>(pressureGradientFaceFieldIds[d]);
  cellGraph->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto faceCenter = cellGraph->faceCenterPosition(faceId);
    auto faceCells = cellGraph->faceCells(faceId);
    for (auto cell : faceCells) {
      points.emplace_back(cellGraph->cellCenterPosition(cell.id));
      pressureValues.emplace_back(pressureField[cell.id]);
    }
    points.emplace_back(faceCenter);
    pressureValues.emplace_back(facePressureField[faceId]);

    if (cellGraph->faceOrientation(faceId) ==
        Definitions::Orientation::VERTICAL)
      (*gradientFields[0])[faceId] =
          rbf->gradientAt(faceCenter, 0, points, pressureValues);
    else if (cellGraph->faceOrientation(faceId) ==
             Definitions::Orientation::HORIZONTAL)
      (*gradientFields[1])[faceId] =
          rbf->gradientAt(faceCenter, 1, points, pressureValues);
    else if (cellGraph->faceOrientation(faceId) ==
             Definitions::Orientation::DEPTH)
      (*gradientFields[2])[faceId] =
          rbf->gradientAt(faceCenter, 2, points, pressureValues);
  });
}

// BOTTLENECK
template <int D>
void SimulationSteps::computePressureGradient(
    StructureInterface<D> *cellGraph, size_t faceMaterialFieldId,
    size_t pressureCellFieldId, size_t cellMaterialFieldId,
    const std::vector<size_t> &pressureGradientFaceFieldIds,
    DifferentialRBF<D> *_rbf) {
  auto &faceMaterialField =
      *cellGraph->template field<Definitions::Material>(faceMaterialFieldId);
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  auto &pressureField = *cellGraph->template field<double>(pressureCellFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        cellGraph->template field<double>(pressureGradientFaceFieldIds[d]);
  iterateGridFaces(cellGraph, [&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto faceCenter = cellGraph->faceCenterPosition(faceId);
    auto faceCells = cellGraph->faceCells(faceId);
    // FIXME: what about domain faces???
    if (cellGraph->cellRegion(faceCells[0].id).size() ==
        cellGraph->cellRegion(faceCells[1].id).size()) {
      for (auto cell : faceCells) {
        points.emplace_back(cellGraph->cellCenterPosition(cell.id));
        pressureValues.emplace_back(pressureField[cell.id]);
      }
    } else {
      std::set<size_t> cellsSet;
      cellsSet.insert(faceCells[0].id);
      cellsSet.insert(faceCells[1].id);
      for (auto cell : faceCells) {
        auto neighborCells = cellGraph->cellNeighbors(cell.id);
        for (auto nCell : neighborCells) {
          if (nCell.id >= 0)
            cellsSet.insert(nCell.id);
        }
      }
      for (auto cell : cellsSet) {
        if (cellMaterialField[cell] != Definitions::Material::SOLID) {
          points.emplace_back(cellGraph->cellCenterPosition(cell));
          pressureValues.emplace_back(pressureField[cell]);
        }
      }
    }

    DifferentialRBF<D> rbf{*_rbf->kernel()};

    if (cellGraph->faceOrientation(faceId) ==
        Definitions::Orientation::VERTICAL)
      (*gradientFields[0])[faceId] =
          rbf.gradientAt(faceCenter, 0, points, pressureValues);
    else if (cellGraph->faceOrientation(faceId) ==
             Definitions::Orientation::HORIZONTAL)
      (*gradientFields[1])[faceId] =
          rbf.gradientAt(faceCenter, 1, points, pressureValues);
    else if (cellGraph->faceOrientation(faceId) ==
             Definitions::Orientation::DEPTH)
      (*gradientFields[2])[faceId] =
          rbf.gradientAt(faceCenter, 2, points, pressureValues);
  });
}

template <int D>
void SimulationSteps::computePressureGradientOnMidPointFaces(
    StructureInterface<D> *structure, size_t cellMaterialFieldId,
    size_t cellPressureFieldId, size_t faceMaterialFieldId,
    const std::vector<size_t> &facePressureGradientFieldIds,
    DifferentialRBF<D> *rbf) {
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);
  auto &pressureField = *structure->template field<double>(cellPressureFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        structure->template field<double>(facePressureGradientFieldIds[d]);
  structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto faceCenter = structure->faceCenterPosition(faceId);
    auto faceCells = structure->faceCells(faceId);
    std::set<size_t> cells;
    for (auto cell : faceCells) {
      cells.insert(cell.id);
      if (cellMaterialField[cell.id] == Definitions::Material::FLUID) {
        auto neighbors = structure->cellNeighbors(cell.id);
        for (auto n : neighbors)
          if (n.id >= 0)
            cells.insert(n.id);
      }
    }
    for (auto c : cells) {
      if (cellMaterialField[c] != Definitions::Material::SOLID) {
        points.emplace_back(structure->cellCenterPosition(c));
        pressureValues.emplace_back(pressureField[c]);
      }
    }
    for (int d = 0; d < D; d++)
      (*gradientFields[d])[faceId] =
          rbf->gradientAt(faceCenter, d, points, pressureValues);
  });
}

template <int D>
void SimulationSteps::computePressureGradient(
    StructureInterface<D> *cellGraph, size_t pressureCellFieldId,
    const std::vector<size_t> &pressureGradientCellFieldIds,
    size_t cellMaterialFieldId, DifferentialRBF<D> *rbf) {
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  auto &pressureField = *cellGraph->template field<double>(pressureCellFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        cellGraph->template field<double>(pressureGradientCellFieldIds[d]);
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto cellCenter = cellGraph->cellCenterPosition(cellId);

    auto neighbors = cellGraph->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0)
        continue;
      points.emplace_back(cellGraph->cellCenterPosition(cell.id));
      pressureValues.emplace_back(pressureField[cell.id]);
    }
    points.emplace_back(cellCenter);
    pressureValues.emplace_back(pressureField[cellId]);

    for (size_t d = 0; d < D; d++)
      (*gradientFields[d])[cellId] =
          rbf->gradientAt(cellCenter, d, points, pressureValues);
  });
}

template <int D>
void SimulationSteps::computePressureGradient(
    StructureInterface<D> *cellGraph,
    const ParticleSystem<double, D> *particles, size_t pressureCellFieldId,
    const std::vector<size_t> &pressureGradientCellFieldIds,
    size_t cellMaterialFieldId, size_t cellParticleCenterFieldId,
    DifferentialRBF<D> *rbf) {
  auto &cellMaterialField =
      *cellGraph->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *cellGraph->template field<Point<double, D>>(cellParticleCenterFieldId);
  auto &pressureField = *cellGraph->template field<double>(pressureCellFieldId);
  auto &cellPressureField =
      *cellGraph->template field<double>(pressureCellFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        cellGraph->template field<double>(pressureGradientCellFieldIds[d]);
  std::set<size_t> airSurfaceCells;
  // Gradients are computed only on fluid cells
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto cellCenter = cellParticleCenterField[cellId];
    auto neighbors = cellGraph->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0 ||
          cellMaterialField[cell.id] == Definitions::Material::SOLID)
        continue;
      points.emplace_back(cellParticleCenterField[cell.id]);
      pressureValues.emplace_back(pressureField[cell.id]);
      if (cellMaterialField[cell.id] == Definitions::Material::AIR)
        airSurfaceCells.insert(cell.id);
    }
    points.emplace_back(cellCenter);
    pressureValues.emplace_back(cellPressureField[cellId]);
    for (size_t d = 0; d < D; d++)
      (*gradientFields[d])[cellId] =
          rbf->gradientAt(cellCenter, d, points, pressureValues);
  });
  // now compute pressure gradient on close air as well
  for (size_t cellId : airSurfaceCells) {
    std::vector<Point<double, D>> points;
    std::vector<double> pressureValues;
    auto cellCenter = cellParticleCenterField[cellId];
    auto neighbors = cellGraph->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0 ||
          cellMaterialField[cell.id] == Definitions::Material::SOLID)
        continue;
      points.emplace_back(cellParticleCenterField[cell.id]);
      pressureValues.emplace_back(pressureField[cell.id]);
    }
    points.emplace_back(cellCenter);
    pressureValues.emplace_back(cellPressureField[cellId]);
    std::cerr << "PRESSURE GRADIENT ON AIR CELL " << cellId << std::endl;
    for (size_t d = 0; d < D; d++)
      (*gradientFields[d])[cellId] =
          rbf->gradientAt(cellCenter, d, points, pressureValues);
  }
}

template <int D>
void SimulationSteps::computePressureGradient(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialFieldId, size_t pressureCellFieldId,
    size_t cellParticleCenterFieldId,
    const std::vector<size_t> &particlePressureGradientFieldIds,
    DifferentialRBF<D> *rbf) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *structure->template field<Point<double, D>>(cellParticleCenterFieldId);
  auto &pressureField = *structure->template field<double>(pressureCellFieldId);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    std::vector<Point<double, D>> points(1, cellParticleCenterField[cellId]);
    std::vector<double> values(1, pressureField[cellId]);
    auto neighbors = structure->cellNeighbors(cellId);
    for (auto cell : neighbors) {
      if (cell.id < 0 ||
          cellMaterialField[cell.id] == Definitions::Material::SOLID)
        continue;
      if (cellMaterialField[cell.id] == Definitions::Material::FLUID)
        points.emplace_back(cellParticleCenterField[cell.id]);
      else {
        Point<double, D> zeroPressurePoint =
            (cellParticleCenterField[cellId] +
             cellParticleCenterField[cell.id]) /
            2.;
        switch (cell.side) {
        case Definitions::Side::LEFT:
        case Definitions::Side::RIGHT:
          zeroPressurePoint[1] = cellParticleCenterField[cellId][1];
          break;
        case Definitions::Side::BOTTOM:
        case Definitions::Side::TOP:
          zeroPressurePoint[0] = cellParticleCenterField[cellId][0];
          break;
        case Definitions::Side::BACK:
        case Definitions::Side::FRONT:
          zeroPressurePoint[2] = cellParticleCenterField[cellId][2];
          break;
        }
        points.emplace_back(zeroPressurePoint);
      }
      values.emplace_back(pressureField[cell.id]);
    }
    particles->iterateParticles(
        structure->cellRegion(cellId), [&](size_t id, Point<double, D> p) {
          for (int d = 0; d < D; d++)
            particles->setScalarProperty(particlePressureGradientFieldIds[d],
                                         id,
                                         rbf->gradientAt(p, d, points, values));
        });
  });
}

template <int D>
void SimulationSteps::correctVelocities(
    ParticleSystem<double, D> *particles,
    const std::vector<size_t> &particlePressureGradientFieldIds,
    const std::vector<size_t> &particleVelocityFieldIds, double dt) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    for (int d = 0; d < D; d++) {
      double v = particles->getScalarProperty(particleVelocityFieldIds[d], id);
      double g =
          particles->getScalarProperty(particlePressureGradientFieldIds[d], id);
      particles->setScalarProperty(particleVelocityFieldIds[d], id, v - dt * g);
    }
  });
}

template <int D>
void SimulationSteps::correctVelocities(
    StructureInterface<D> *structure,
    const std::vector<size_t> &pressureGradientFaceFieldIds,
    const std::vector<size_t> &velocityFaceFieldIds, size_t faceMaterialFieldId,
    double dt) {
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);

  std::vector<typename StructureInterface<D>::template Field<double> *>
      pressureGradientFields(D, nullptr), velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++) {
    pressureGradientFields[d] =
        structure->template field<double>(pressureGradientFaceFieldIds[d]);
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  }
  // Correct velocities through pressure gradient
  structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    if (structure->faceOrientation(faceId) ==
        Definitions::Orientation::VERTICAL)
      (*velocityFields[0])[faceId] -= dt * (*pressureGradientFields[0])[faceId];
    else if (structure->faceOrientation(faceId) ==
             Definitions::Orientation::HORIZONTAL)
      (*velocityFields[1])[faceId] -= dt * (*pressureGradientFields[1])[faceId];
    else
      (*velocityFields[2])[faceId] -= dt * (*pressureGradientFields[2])[faceId];
  });
}

template <int D>
void SimulationSteps::correctVelocitiesOnMidPointFaces(
    StructureInterface<D> *structure,
    const std::vector<size_t> &pressureGradientFaceFieldIds,
    const std::vector<size_t> &velocityFaceFieldIds, size_t faceMaterialFieldId,
    double dt) {
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);

  std::vector<typename StructureInterface<D>::template Field<double> *>
      pressureGradientFields(D, nullptr), velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++) {
    pressureGradientFields[d] =
        structure->template field<double>(pressureGradientFaceFieldIds[d]);
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  }
  // Correct velocities through pressure gradient
  structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    for (int d = 0; d < D; d++)
      (*velocityFields[d])[faceId] -= dt * (*pressureGradientFields[d])[faceId];
  });
}

template <int D>
void SimulationSteps::transferCellVelocitiesToParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialFieldId, const std::vector<size_t> &cellVelocityFieldIds,
    const std::vector<size_t> &particleVelocityFieldIds,
    Interpolant<Point<double, D>> *interpolator) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      cellVelocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    cellVelocityFields[d] =
        structure->template field<double>(cellVelocityFieldIds[d]);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      auto neighbors = structure->cellNeighbors(cellId);
      std::vector<Point<double, D>> points(
          1, structure->cellCenterPosition(cellId));
      std::vector<double> values[D];
      for (int d = 0; d < D; d++)
        values[d].emplace_back((*cellVelocityFields[d])[cellId]);
      auto cellRegion = structure->cellRegion(cellId);
      for (auto neighbor : neighbors) {
        if (neighbor.id >= 0 &&
            cellMaterialField[neighbor.id] == Definitions::Material::FLUID) {
          points.emplace_back(structure->cellCenterPosition(neighbor.id));
          for (int d = 0; d < D; d++)
            values[d].emplace_back((*cellVelocityFields[d])[neighbor.id]);
        }
      }
      std::vector<size_t> particleIds;
      particles->iterateParticles(
          cellRegion,
          [&](size_t id, Point<double, D> p) { particleIds.emplace_back(id); });
      for (int d = 0; d < D; d++) {
        interpolator->computeWeights(points, values[d]);
        for (auto id : particleIds) {
          particles->setScalarProperty(
              particleVelocityFieldIds[d], id,
              interpolator->evaluate((*particles)[id]));
        }
      }
    }
  });
}

template <int D>
void SimulationSteps::correctVelocities(
    StructureInterface<D> *structure, size_t pressureGradientCellFieldId,
    ParticleSystem<double, D> *particles, size_t propertyId,
    const std::function<bool(size_t)> &useElement, double dt) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    double value =
        structure->sampleCellField(pressureGradientCellFieldId, p, useElement);
    particles->setScalarProperty(propertyId, id,
                                 particles->getScalarProperty(propertyId, id) -
                                     dt * value);
  });
}

template <int D>
void SimulationSteps::correctVelocities(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    const std::vector<size_t> &pressureGradientCellFieldIds,
    size_t surfaceMaskCellFieldId, size_t cellMaterialFieldId,
    size_t cellParticleCenterFieldId,
    const std::vector<size_t> &particleVelocityIds, double dt,
    Interpolant<Point<double, D>> *interpolator) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellParticleCenterField =
      *structure->template field<Point<double, D>>(cellParticleCenterFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      gradientFields(D, nullptr);

  for (size_t d = 0; d < D; d++)
    gradientFields[d] =
        structure->template field<double>(pressureGradientCellFieldIds[d]);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      std::vector<Point<double, D>> points(1, cellParticleCenterField[cellId]);
      std::vector<double> values[D];
      for (int d = 0; d < D; d++)
        values[d].emplace_back((*gradientFields[d])[cellId]);
      // auto neighbors = structure->cellNeighbors(cellId);
      // for (auto neighbor : neighbors) {
      //   if (neighbor.id <= 0)
      //     continue;
      structure->iterateCellRing(cellId, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::FLUID) {
          // if (cellMaterialField[nid] != Definitions::Material::SOLID) {
          points.emplace_back(cellParticleCenterField[nid]);
          for (int d = 0; d < D; d++)
            values[d].emplace_back((*gradientFields[d])[nid]);
        }
      });

      // if (cellId == 399 || cellId == 91 || cellId == 235 || cellId ==
      // 59 || cellId == 229) {
      std::cerr << "#####################################" << std::endl;
      std::cerr << "#####################################" << std::endl;
      std::cerr << "c" << cellId << std::endl;
      for (int i = 0; i < points.size(); i++) {
        std::cerr << points[i] << ":  \t";
        for (int d = 0; d < D; d++)
          std::cerr << values[d][i] << ' ';
        std::cerr << std::endl;
      }
      // }
      RBFInterpolant<Point<double, D>> rbf[D];
      for (int d = 0; d < D; d++)
        rbf[d].computeWeights(points, values[d]);
      particles->iterateParticles(
          structure->cellRegion(cellId), [&](size_t id, Point<double, D> p) {
            std::cerr << "Target: " << p << std::endl;
            double correction[D];
            for (int d = 0; d < D; d++) {
              double value = interpolator->interpolateAt(p, points, values[d]);
              // double value = rbf[d].evaluate(p);
              particles->setScalarProperty(
                  particleVelocityIds[d], id,
                  particles->getScalarProperty(particleVelocityIds[d], id) -
                      dt * value);
              correction[d] = value;
            }
            std::cerr << "Correction: " << correction[0] << " "
                      << correction[1];
            std::cerr << std::endl;
          });
    }
  });
}

template <>
inline void SimulationSteps::correctVelocitiesOnSurface(
    StructureInterface<2> *structure, ParticleSystem<double, 2> *particles,
    const std::vector<size_t> &pressureGradientCellFieldIds,
    size_t surfaceMaskCellFieldId, size_t cellMaterialFieldId,
    size_t cellParticleCenterFieldId,
    const std::vector<size_t> &particleVelocityIds, double dt,
    Interpolant<Point<double, 2>> *interpolator) {
  int D = 2;
  // for each cell, computes its stencil and updates its particles
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &cellParticleCenterField =
      *cellGraph->template field<Point<double, 2>>(cellParticleCenterFieldId);
  std::vector<typename StructureInterface<2>::template Field<double> *>
      gradientFields(D, nullptr);
  for (int d = 0; d < D; d++)
    gradientFields[d] =
        cellGraph->template field<double>(pressureGradientCellFieldIds[d]);
  std::vector<double> minV(D, INFINITY), maxV(D, -INFINITY);
  cellGraph->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID &&
        surfaceField[cellId]) {
      auto centerCell = cellParticleCenterField[cellId];
      auto cellRegion = cellGraph->cellRegion(cellId);
      auto neighbors = cellGraph->cellNeighbors(cellId);
      std::vector<bool> useVertex(D == 2 ? 4 : 8, false);
      std::vector<std::vector<double>> vertexValue(
          D, std::vector<double>(D == 2 ? 4 : 8, 0.));
      //   0    1     2     3     4    5
      // LEFT RIGHT BOTTOM TOP  BACK FRONT
      std::vector<int> neighborIds(2 * D, -1);
      THROW(neighbors.size() == 2 * D,
            "SimulationSteps::correctVelocitiesOnSurface surface with wrong "
            "number of neighbors");
      for (auto neighbor : neighbors) {
        switch (neighbor.side) {
        case Definitions::Side::LEFT:
          neighborIds[0] = neighbor.id;
          break;
        case Definitions::Side::RIGHT:
          neighborIds[1] = neighbor.id;
          break;
        case Definitions::Side::BOTTOM:
          neighborIds[2] = neighbor.id;
          break;
        case Definitions::Side::TOP:
          neighborIds[3] = neighbor.id;
          break;
        case Definitions::Side::BACK:
          neighborIds[4] = neighbor.id;
          break;
        case Definitions::Side::FRONT:
          neighborIds[5] = neighbor.id;
          break;
        default:
          THROW(false,
                "SimulationSteps::correctVelocities invalid neighbor side");
        }
      }
      for (size_t i = 0; i < 2 * D; i++) {
        if (neighborIds[i] < 0 ||
            cellMaterialField[neighborIds[i]] != Definitions::Material::FLUID) {
          switch (i) {
            //   0    1     2     3     4    5
            // LEFT RIGHT BOTTOM TOP  BACK FRONT
          case 0: // Definitions::Side::LEFT:
            if (D == 2) {
              // interpolate with bottom neighbor
              // free dimensions
              size_t dimensions[1] = {1};
              if (neighborIds[2] >= 0 && cellMaterialField[neighborIds[2]] ==
                                             Definitions::Material::FLUID) {
                useVertex[0] = true;
                for (auto d : dimensions) {
                  double w1 = centerCell[d] - cellRegion.lower()[d];
                  double w2 = cellRegion.lower()[d] -
                              cellParticleCenterField[neighborIds[2]][d];
                  vertexValue[d][0] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[2]]) /
                      (w1 + w2);
                }
              }
              if (neighborIds[3] >= 0 && cellMaterialField[neighborIds[3]] ==
                                             Definitions::Material::FLUID) {
                useVertex[3] = true;
                for (auto d : dimensions) {
                  double w1 = cellRegion.upper()[d] - centerCell[d];
                  double w2 = cellParticleCenterField[neighborIds[3]][d] -
                              cellRegion.upper()[d];
                  vertexValue[d][3] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[3]]) /
                      (w1 + w2);
                }
              }
            } else {
              NOT_IMPLEMENTED();
              useVertex[0] = useVertex[2] = useVertex[4] = useVertex[6] = true;
            }
            break;
            //   0    1     2     3     4    5
            // LEFT RIGHT BOTTOM TOP  BACK FRONT
          case 1: // Definitions::Side::RIGHT:
            if (D == 2) {
              // interpolate with bottom neighbor
              // free dimensions
              size_t dimensions[1] = {1};
              if (neighborIds[2] >= 0 && cellMaterialField[neighborIds[2]] ==
                                             Definitions::Material::FLUID) {
                useVertex[1] = true;
                for (auto d : dimensions) {
                  double w1 = centerCell[d] - cellRegion.lower()[d];
                  double w2 = cellRegion.lower()[d] -
                              cellParticleCenterField[neighborIds[2]][d];
                  vertexValue[d][1] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[2]]) /
                      (w1 + w2);
                }
              }
              if (neighborIds[3] >= 0 && cellMaterialField[neighborIds[3]] ==
                                             Definitions::Material::FLUID) {
                useVertex[2] = true;
                for (auto d : dimensions) {
                  double w1 = cellRegion.upper()[d] - centerCell[d];
                  double w2 = cellParticleCenterField[neighborIds[3]][d] -
                              cellRegion.upper()[d];
                  vertexValue[d][2] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[3]]) /
                      (w1 + w2);
                }
              }
            } else {
              NOT_IMPLEMENTED();
              useVertex[0] = useVertex[2] = useVertex[4] = useVertex[6] = true;
            }
            break;
            //   0    1     2     3     4    5
            // LEFT RIGHT BOTTOM TOP  BACK FRONT
          case 2: // Definitions::Side::BOTTOM:
            if (D == 2) {
              // interpolate with bottom neighbor
              // free dimensions
              size_t dimensions[1] = {0};
              if (neighborIds[0] >= 0 && cellMaterialField[neighborIds[0]] ==
                                             Definitions::Material::FLUID) {
                useVertex[0] = true;
                for (auto d : dimensions) {
                  double w1 = centerCell[d] - cellRegion.lower()[d];
                  double w2 = cellRegion.lower()[d] -
                              cellParticleCenterField[neighborIds[0]][d];
                  vertexValue[d][0] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[0]]) /
                      (w1 + w2);
                }
              }
              if (neighborIds[1] >= 0 && cellMaterialField[neighborIds[1]] ==
                                             Definitions::Material::FLUID) {
                useVertex[1] = true;
                for (auto d : dimensions) {
                  double w1 = cellRegion.upper()[d] - centerCell[d];
                  double w2 = cellParticleCenterField[neighborIds[1]][d] -
                              cellRegion.upper()[d];
                  vertexValue[d][1] =
                      (w1 * (*gradientFields[d])[cellId] +
                       w2 * (*gradientFields[d])[neighborIds[1]]) /
                      (w1 + w2);
                }
              }
            } else {
              NOT_IMPLEMENTED();
              useVertex[0] = useVertex[2] = useVertex[4] = useVertex[6] = true;
            }
            break;
            //   0    1     2     3     4    5
            // LEFT RIGHT BOTTOM TOP  BACK FRONT
          case 3: // Definitions::Side::TOP:            if (D == 2) {
            // interpolate with bottom neighbor
            // free dimensions
            size_t dimensions[1] = {0};
            if (neighborIds[0] >= 0 && cellMaterialField[neighborIds[0]] ==
                                           Definitions::Material::FLUID) {
              useVertex[3] = true;
              for (auto d : dimensions) {
                double w1 = centerCell[d] - cellRegion.lower()[d];
                double w2 = cellRegion.lower()[d] -
                            cellParticleCenterField[neighborIds[0]][d];
                vertexValue[d][3] =
                    (w1 * (*gradientFields[d])[cellId] +
                     w2 * (*gradientFields[d])[neighborIds[0]]) /
                    (w1 + w2);
              }
            }
            if (neighborIds[1] >= 0 && cellMaterialField[neighborIds[1]] ==
                                           Definitions::Material::FLUID) {
              useVertex[2] = true;
              for (auto d : dimensions) {
                double w1 = cellRegion.upper()[d] - centerCell[d];
                double w2 = cellParticleCenterField[neighborIds[1]][d] -
                            cellRegion.upper()[d];
                vertexValue[d][2] =
                    (w1 * (*gradientFields[d])[cellId] +
                     w2 * (*gradientFields[d])[neighborIds[1]]) /
                    (w1 + w2);
              }
            } else {
              NOT_IMPLEMENTED();
              useVertex[0] = useVertex[2] = useVertex[4] = useVertex[6] = true;
            }
            break;
          }
        }
      }
      std::vector<Point<double, 2>> points(1, centerCell);
      std::vector<double> values[2];
      for (size_t d = 0; d < D; d++)
        values[d].emplace_back((*gradientFields[d])[cellId]);
      for (auto nid : neighborIds) {
        if (nid >= 0 &&
            cellMaterialField[nid] == Definitions::Material::FLUID) {
          points.emplace_back(cellParticleCenterField[nid]);
          for (size_t d = 0; d < D; d++)
            values[d].emplace_back((*gradientFields[d])[nid]);
        }
      }
      for (size_t i = 0; i < 4; i++)
        if (useVertex[i]) {
          points.emplace_back(
              structure->vertexPosition(cellGraph->node(cellId).vertex(i)));
          for (size_t d = 0; d < D; d++)
            values[d].emplace_back(vertexValue[d][i]);
        }

      if (cellId == 118 || cellId == 34 || cellId == 153 || cellId == 59 ||
          cellId == 229) {
        std::cerr << "Surface c" << cellId << std::endl;
        std::cerr << "#####################################" << std::endl;
        for (int i = 0; i < points.size(); i++) {
          std::cerr << points[i] << ":  \t";
          for (int d = 0; d < D; d++)
            std::cerr << values[d][i] << ' ';
          std::cerr << std::endl;
        }
      }
      particles->iterateParticles(
          cellGraph->cellRegion(cellId), [&](size_t id, Point<double, 2> p) {
            if (cellId == 399 || cellId == 91 || cellId == 156 ||
                cellId == 59 || cellId == 229) {
              std::cerr << "Particle: " << p << std::endl;
              std::cerr << "Before: ";
              for (int d = 0; d < D; d++)
                std::cerr << particles->getScalarProperty(
                                 particleVelocityIds[d], id)
                          << " ";
              std::cerr << std::endl;
            }
            for (size_t d = 0; d < D; d++)
              particles->setScalarProperty(
                  particleVelocityIds[d], id,
                  particles->getScalarProperty(particleVelocityIds[d], id) -
                      dt * interpolator->interpolateAt(p, points, values[d]));
            if (cellId == 399 || cellId == 91 || cellId == 156 ||
                cellId == 59 || cellId == 229) {
              std::cerr << "After: ";
              for (int d = 0; d < D; d++)
                std::cerr << particles->getScalarProperty(
                                 particleVelocityIds[d], id)
                          << " ";
              std::cerr << std::endl << std::endl;
            }
          });
    }
  });
}

template <int D>
void SimulationSteps::reseedParticles(const StructureInterface<D> *structure,
                                      size_t cellMaterialFieldId,
                                      size_t surfaceMaskCellFieldId,
                                      ParticleSystem<double, D> *particles,
                                      size_t particleGroupPropertyId) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  static BoxSampler sampler;
  std::vector<size_t> toDelete;
  std::vector<Point<double, D>> toCreate;
  std::vector<int> particleGroups;
  size_t minNumber = D == 2 ? 3 : 5;
  size_t maxNumber = D == 2 ? 12 : 15;
  size_t reseedNumber = D == 2 ? 4 : 8;
  structure->iterateCells([&](size_t id) {
    // allow reseed on surface cells touching walls
    if (cellMaterialField[id] != Definitions::Material::FLUID)
      return;
    bool airSurface = false;
    bool solidWall = false;
    structure->iterateCellRing(id, [&](size_t nid) {
      if (cellMaterialField[nid] == Definitions::Material::AIR)
        airSurface = true;
      if (cellMaterialField[nid] == Definitions::Material::SOLID)
        solidWall = true;
    });
    if (surfaceField[id]) {
      if (airSurface)
        return;
    }
    // if (solidWall)
    //  return;
    auto region = structure->cellRegion(id);
    std::vector<size_t> cellParticles;
    auto cellLevel = [&](size_t cellId) -> size_t {
      if (D == 2)
        return dynamic_cast<const CellGraph2 *>(structure)->node(id).level();
      return dynamic_cast<const CellGraph3 *>(structure)->node(id).level();
    };
    size_t level = cellLevel(id);
    particles->iterateParticlesZ(region, level,
                                 [&](size_t i, Point<double, D> p) {
                                   cellParticles.emplace_back(i);
                                   UNUSED_VARIABLE(p);
                                 });
    if (cellParticles.size() > maxNumber || cellParticles.size() < minNumber) {
      int newGroup = 0;
      if (cellParticles.size())
        newGroup = particles->getScalarProperty(particleGroupPropertyId,
                                                cellParticles[0]);
      for (auto p : cellParticles)
        toDelete.emplace_back(p);
      for (size_t i = 0; i < reseedNumber; i++) {
        toCreate.emplace_back(sampler.sample(region));
        particleGroups.emplace_back(newGroup);
      }
    }
  });
  // test: force ressed on walls
  // structure->iterateCells([&](size_t id) {
  //   if (cellMaterialField[id] != Definitions::Material::FLUID ||
  //       surfaceField[id])
  //     return;
  //   bool solidWall = false;
  //   structure->iterateCellRing(id, [&](size_t nid) {
  //     if (cellMaterialField[nid] == Definitions::Material::SOLID)
  //       solidWall = true;
  //   });
  //   if (solidWall) {
  //     auto region = structure->cellRegion(id);
  //     particles->iterateParticles(region, [&](size_t i, Point<double, D> p)
  //     {
  //       toDelete.emplace_back(i);
  //     });
  //     for (size_t i = 0; i < reseedNumber; i++)
  //       toCreate.emplace_back(sampler.sample(region));
  //   }
  // });
  for (auto p : toDelete)
    particles->remove(p);
  for (int i = 0; i < toCreate.size(); ++i) {
    int id = particles->addParticle(toCreate[i]);
    particles->setScalarProperty(particleGroupPropertyId, id,
                                 particleGroups[i]);
  }

#ifdef FUROO_PROFILE
  CodeProfiler::count("ParticlesCreated", toCreate.size());
  CodeProfiler::count("ParticlesDeleted", toDelete.size());
#endif
}

template <int D>
void SimulationSteps::reseedParticles(
    const StructureInterface<D> *structure, size_t cellMaterialFieldId,
    size_t surfaceMaskCellFieldId, ParticleSystem<double, D> *particles,
    const std::vector<size_t> &particleVelocityFieldIds) {
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  static BoxSampler sampler;
  std::vector<size_t> toDelete;
  std::vector<Point<double, D>> toCreate;
  size_t minNumber = D == 2 ? 3 : 5;
  size_t maxNumber = D == 2 ? 12 : 15;
  size_t reseedNumber = D == 2 ? 4 : 8;
  structure->iterateCells([&](size_t id) {
    // allow reseed on surface cells touching walls
    if (cellMaterialField[id] != Definitions::Material::FLUID)
      return;
    if (surfaceField[id]) {
      // return;
      bool airSurface = false;
      structure->iterateCellRing(id, [&](size_t nid) {
        if (cellMaterialField[nid] == Definitions::Material::AIR)
          airSurface = true;
      });
      if (airSurface)
        return;
    }
    // bool notFluid = false;
    // structure->iterateCellRing(id, [&](size_t nid) {
    //   if (cellMaterialField[nid] == Definitions::Material::SOLID)
    //     notFluid = true;
    // });
    // if (notFluid)
    //   return;
    auto region = structure->cellRegion(id);
    std::vector<size_t> cellParticles;
    auto cellLevel = [&](size_t cellId) -> size_t {
      if (D == 2)
        return dynamic_cast<const CellGraph2 *>(structure)->node(id).level();
      return dynamic_cast<const CellGraph3 *>(structure)->node(id).level();
    };
    size_t level = cellLevel(id);
    particles->iterateParticlesZ(region, level,
                                 [&](size_t i, Point<double, D> p) {
                                   cellParticles.emplace_back(i);
                                   UNUSED_VARIABLE(p);
                                 });
    if (cellParticles.size() > maxNumber || cellParticles.size() < minNumber) {
      for (auto p : cellParticles)
        toDelete.emplace_back(p);
      for (size_t i = 0; i < reseedNumber; i++)
        toCreate.emplace_back(sampler.sample(region));
    }
  });
  std::vector<double> newVelocities[D];
  for (auto p : toCreate)
    for (size_t d = 0; d < D; d++)
      newVelocities[d].emplace_back(
          particles->sampleScalarProperty(particleVelocityFieldIds[d], p, 8));
  for (auto p : toDelete)
    particles->remove(p);
  for (size_t i = 0; i < toCreate.size(); i++) {
    size_t id = particles->addParticle(toCreate[i]);
    for (size_t d = 0; d < D; d++)
      particles->setScalarProperty(particleVelocityFieldIds[d], id,
                                   newVelocities[d][i]);
  }

#ifdef FUROO_PROFILE
  CodeProfiler::count("ParticlesCreated", toCreate.size());
  CodeProfiler::count("ParticlesDeleted", toDelete.size());
#endif
}

template <int D>
void SimulationSteps::transferFromGridToParticles(
    const StructureInterface<D> *grid, size_t fieldId,
    ParticleSystem<double, D> *particles, size_t propertyId,
    const std::function<bool(size_t)> &useElement,
    const std::function<bool(size_t)> &useCell,
    Definitions::BoundaryVelocity boundaryVelocity) {
  particles->iterateParticles([&](size_t id, Point<double, D> p) {
    double value = 0.;
    switch (grid->fieldLocation(fieldId)) {
    case Definitions::MeshLocation::FACE_CENTER:
      NOT_IMPLEMENTED();
      break;
    case Definitions::MeshLocation::VERTICAL_FACE_CENTER:
      value =
          grid->sampleFaceField(fieldId, p, Definitions::Orientation::VERTICAL,
                                useElement, useCell, boundaryVelocity);
      break;
    case Definitions::MeshLocation::HORIZONTAL_FACE_CENTER:
      value = grid->sampleFaceField(fieldId, p,
                                    Definitions::Orientation::HORIZONTAL,
                                    useElement, useCell, boundaryVelocity);
      break;
    case Definitions::MeshLocation::DEPTH_FACE_CENTER:
      value = grid->sampleFaceField(fieldId, p, Definitions::Orientation::DEPTH,
                                    useElement, useCell, boundaryVelocity);
      break;
    case Definitions::MeshLocation::VERTEX_CENTER:
      value = grid->sampleVertexField(fieldId, p);
      break;
    case Definitions::MeshLocation::CELL_CENTER:
      value = grid->sampleCellField(fieldId, p, useElement);
      break;
    default:
      THROW(false, "SimulationSteps::transferFromGridToParticles invalid field "
                   "location");
    }
    particles->setScalarProperty(propertyId, id, value);
  });
}

// BOTTLENECK
template <int D>
void SimulationSteps::transferFromGridToParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialTypeFieldId, size_t faceMaterialTypeFieldId,
    size_t cellSurfaceMaskFieldId,
    const std::vector<size_t> &faceVelocityFieldIds,
    const std::vector<size_t> &particleVelocityFieldIds,
    Interpolant<Point<double, D>> *_interpolator) {
  auto &cellMaterialField = *structure->template field<Definitions::Material>(
      cellMaterialTypeFieldId);
  auto &faceMaterialField = *structure->template field<Definitions::Material>(
      faceMaterialTypeFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);

  const_cast<ParticleSystem<double, D> *>(particles)
      ->update(); // invoked here to avoid race conditions with OpenMP
#pragma omp parallel
  {
    iterateGridCells(structure, [&](size_t cellId) {
      QuinticKernel<Point<double, D>> kernel;
      RBFInterpolant<Point<double, D>> interpolator(kernel);
      if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
        std::set<size_t> faceIds;

        //   0    1    2     3     4     5
        // LEFT RIGHT TOP BOTTOM FRONT BACK
        // std::vector<bool> notFluid(D * 2, false);
        auto cellFaces = structure->cellFaces(cellId);
        for (auto f : cellFaces)
          faceIds.insert(f);
        // auto neighbors = structure->cellNeighbors(cellId);
        structure->iterateNeighborCells(cellId, 1, [&](size_t nid) {
          // for (auto n : neighbors) {
          if (cellMaterialField[nid] != Definitions::Material::FLUID) {
            // notFluid[static_cast<int>(n.side)] = true;
            // continue;
            return;
          }
          auto _cellFaces = structure->cellFaces(nid);
          for (auto f : _cellFaces)
            if (faceMaterialField[f] == Definitions::Material::FLUID)
              faceIds.insert(f);
          // }
        });
        // setup stencil
        std::vector<Point<double, D>> points[D];
        std::vector<double> values[D];
        for (auto f : faceIds) {
          auto fOrientation = structure->faceOrientation(f);
          int dimension = 0;
          if (fOrientation == Definitions::Orientation::HORIZONTAL)
            dimension = 1;
          else if (fOrientation == Definitions::Orientation::DEPTH)
            dimension = 2;
          points[dimension].emplace_back(structure->faceCenterPosition(f));
          if (faceMaterialField[f] == Definitions::Material::SOLID) {
            // If face is a solid, replace its velocity to satisfy free slip
            // condition: average from the other faces that are in opposite
            // orientations
            values[dimension].emplace_back(0);
          } else {
            values[dimension].emplace_back((*velocityFields[dimension])[f]);
          }
        }
        auto cellCenter = structure->cellCenterPosition(cellId);
        size_t level = 0;
        if (D == 2)
          level = dynamic_cast<CellGraph2 *>(structure)->node(cellId).level();
        else
          level = dynamic_cast<CellGraph3 *>(structure)->node(cellId).level();
        // for (int d = 0; d < D; d++)
        std::vector<Point<double, D>> particlesPosition;
        std::vector<size_t> particlesId;
        particles->iterateParticlesZ(structure->cellRegion(cellId), level,
                                     [&](size_t id, Point<double, D> p) {
                                       particlesId.emplace_back(id);
                                       particlesPosition.emplace_back(p);
                                     });
        for (int d = 0; d < D; d++) {
          std::vector<double> interpValues = interpolator.interpolateAt(
              particlesPosition, points[d], values[d]);
          for (int pi = 0; pi < particlesPosition.size(); pi++) {
            auto id = particlesId[pi];
            particles->setScalarProperty(particleVelocityFieldIds[d], id,
                                         interpValues[pi]);
          }
        }
      }
    });
  }
}

template <int D>
double SimulationSteps::cfl(const StructureInterface<D> *structure,
                            const std::vector<size_t> &fieldIds, double dx) {
  NOT_IMPLEMENTED();
  double maxVelocity = 0.;
  for (size_t d = 0; d < D; d++)
    std::max(maxVelocity,
             std::fabs(structure->maxScalarFieldValue(fieldIds[d])));
  if (maxVelocity == 0.)
    return INFINITY;
  return dx / maxVelocity;
  return 0;
}

template <>
inline void SimulationSteps::handleInjectors(
    StructureInterface<2> *structure, size_t cellMaterialFieldId,
    size_t surfaceMaskCellFieldId, size_t surfaceLevel,
    ParticleSystem<double, 2> *particles,
    const std::vector<furoo::ParticleInjector<2>> &injectors) {
  auto *cellGraph = dynamic_cast<CellGraph2 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  // queue for every air cell that intersects any injector
  std::queue<size_t> fluidCandidates;
  cellGraph->iterateCells([&](size_t cellId) {
    auto cellRegion = cellGraph->cellRegion(cellId);
    if (cellMaterialField[cellId] != Definitions::Material::AIR)
      return;
    for (auto &injector : injectors)
      if (cellRegion.intersect(injector.region())) {
        fluidCandidates.push(cellId);
        break;
      }
  });
  std::queue<size_t> surfaceCandidates;
  // detect new fluid cells and refine them if necessary
  while (!fluidCandidates.empty()) {
    size_t cellId = fluidCandidates.front();
    fluidCandidates.pop();
    auto cellRegion = cellGraph->cellRegion(cellId);
    if (particles->containsParticle(cellRegion)) {
      if (cellGraph->node(cellId).level() < surfaceLevel) {
        auto children = cellGraph->refine(cellId);
        for (auto child : children)
          fluidCandidates.push(child);
        continue;
      }
      cellMaterialField[cellId] = Definitions::Material::FLUID;
      surfaceCandidates.push(cellId);
    }
  }
  // update surface
  while (!surfaceCandidates.empty()) {
    size_t cellId = surfaceCandidates.front();
    surfaceCandidates.pop();
    size_t neighborCount = 0;
    cellGraph->iterateCellRing(cellId, [&](size_t neighborId) {
      neighborCount++;
      if (cellMaterialField[neighborId] != Definitions::Material::FLUID)
        surfaceField[cellId] = 1;
    });
  }
}

template <>
inline void SimulationSteps::handleInjectors(
    StructureInterface<3> *structure, size_t cellMaterialFieldId,
    size_t surfaceMaskCellFieldId, size_t surfaceLevel,
    ParticleSystem<double, 3> *particles,
    const std::vector<furoo::ParticleInjector<3>> &injectors) {
  auto *cellGraph = dynamic_cast<CellGraph3 *>(structure);
  auto &surfaceField =
      *cellGraph->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  // queue for every air cell that intersects any injector
  std::queue<size_t> fluidCandidates;
  cellGraph->iterateCells([&](size_t cellId) {
    auto cellRegion = cellGraph->cellRegion(cellId);
    if (cellMaterialField[cellId] != Definitions::Material::AIR)
      return;
    for (auto &injector : injectors)
      if (cellRegion.intersect(injector.region())) {
        fluidCandidates.push(cellId);
        break;
      }
  });
  std::queue<size_t> surfaceCandidates;
  // detect new fluid cells and refine them if necessary
  while (!fluidCandidates.empty()) {
    size_t cellId = fluidCandidates.front();
    fluidCandidates.pop();
    auto cellRegion = cellGraph->cellRegion(cellId);
    auto level = cellGraph->node(cellId).level();
    if (particles->containsParticleZ(cellRegion, level)) {
      if (cellGraph->node(cellId).level() < surfaceLevel) {
        auto children = cellGraph->refine(cellId);
        for (auto child : children)
          fluidCandidates.push(child);
        continue;
      }
      cellMaterialField[cellId] = Definitions::Material::FLUID;
      surfaceCandidates.push(cellId);
    }
  }
  // update surface
  while (!surfaceCandidates.empty()) {
    size_t cellId = surfaceCandidates.front();
    surfaceCandidates.pop();
    size_t neighborCount = 0;
    cellGraph->iterateCellRing(cellId, [&](size_t neighborId) {
      neighborCount++;
      if (cellMaterialField[neighborId] != Definitions::Material::FLUID)
        surfaceField[cellId] = 1;
    });
  }
}

template <int D>
void SimulationSteps::propagateVelocities(
    StructureInterface<D> *structure,
    const std::vector<size_t> &velocityFaceFieldIds, size_t cellMaterialFieldId,
    size_t faceMaterialFieldId, size_t surfaceMaskCellFieldId) {
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(velocityFaceFieldIds[d]);
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &faceMaterialField =
      *structure->template field<Definitions::Material>(faceMaterialFieldId);
  structure->iterateCells([&](size_t cellId) {
    if (!surfaceField[cellId])
      return;
    //   0    1     2     3    4    5
    // LEFT RIGHT BOTTOM TOP BACK FRONT
    size_t faceIds[/*D*/ 3 * 2];
    auto faces = structure->cellFaces(cellId);
    THROW(faces.size() == D * 2, "SimulationSteps<>::propagateVelocities "
                                 "neighbor with wrong number of faces");
    for (auto faceId : faces) {
      auto cells = structure->faceCells(faceId);
      // neighbor side
      auto side = cells[0].side;
      if (cells[0].id != cellId)
        side = Definitions::oppositeSide(side);
      // side is the side our neighbor is in respect to the face
      // so its opposite indicates which face we are dealing with
      side = Definitions::oppositeSide(side);
      switch (side) {
      case Definitions::Side::LEFT:
        faceIds[0] = static_cast<size_t>(faceId);
        break;
      case Definitions::Side::RIGHT:
        faceIds[1] = static_cast<size_t>(faceId);
        break;
      case Definitions::Side::BOTTOM:
        faceIds[2] = static_cast<size_t>(faceId);
        break;
      case Definitions::Side::TOP:
        faceIds[3] = static_cast<size_t>(faceId);
        break;
      case Definitions::Side::BACK:
        faceIds[4] = static_cast<size_t>(faceId);
        break;
      case Definitions::Side::FRONT:
        faceIds[5] = static_cast<size_t>(faceId);
        break;
      default:
        break;
      }
    }

    auto neighbors = structure->cellNeighbors(cellId);
    for (auto neighbor : neighbors) {
      if (neighbor.id < 0 ||
          cellMaterialField[neighbor.id] != Definitions::Material::AIR)
        continue;
      // iterate neighbor faces
      auto neighborFaces = structure->cellFaces(neighbor.id);
      THROW(neighborFaces.size() == D * 2,
            "SimulationSteps<>::propagateVelocities "
            "neighbor with wrong number of faces");
      //   0    1     2     3    4    5
      // LEFT RIGHT BOTTOM TOP BACK FRONT
      size_t neighborFaceIds[/*D*/ 3 * 2];
      size_t components[6] = {0, 0, 1, 1, 2, 2};
      for (auto neighborFaceId : neighborFaces) {
        auto cells = structure->faceCells(neighborFaceId);
        // neighbor side
        auto side = cells[0].side;
        if (cells[0].id != neighbor.id)
          side = Definitions::oppositeSide(side);
        // side is the side our neighbor is in respect to the face
        // so its opposite indicates which face we are dealing with
        side = Definitions::oppositeSide(side);
        switch (side) {
        case Definitions::Side::LEFT:
          neighborFaceIds[0] = static_cast<size_t>(neighborFaceId);
          break;
        case Definitions::Side::RIGHT:
          neighborFaceIds[1] = static_cast<size_t>(neighborFaceId);
          break;
        case Definitions::Side::BOTTOM:
          neighborFaceIds[2] = static_cast<size_t>(neighborFaceId);
          break;
        case Definitions::Side::TOP:
          neighborFaceIds[3] = static_cast<size_t>(neighborFaceId);
          break;
        case Definitions::Side::BACK:
          neighborFaceIds[4] = static_cast<size_t>(neighborFaceId);
          break;
        case Definitions::Side::FRONT:
          neighborFaceIds[5] = static_cast<size_t>(neighborFaceId);
          break;
        default:
          break;
        }
      }
      switch (neighbor.side) {
      case Definitions::Side::LEFT:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 1 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      case Definitions::Side::RIGHT:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 0 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      case Definitions::Side::BOTTOM:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 3 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      case Definitions::Side::TOP:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 2 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      case Definitions::Side::BACK:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 5 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      case Definitions::Side::FRONT:
        for (size_t i = 0; i < D * 2; i++)
          if (i != 4 && faceMaterialField[neighborFaceIds[i]] ==
                            Definitions::Material::AIR)
            (*velocityFields[components[i]])[neighborFaceIds[i]] =
                (*velocityFields[components[i]])[faceIds[i]];
        break;
      default:
        break;
      }
    }
  });
}

template <int D>
void SimulationSteps::computeCellParticlesCentroid(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t surfaceMaskCellFieldId, size_t cellMaterialFieldId,
    size_t cellParticleCenterFieldId) {
  auto &surfaceField =
      *structure->template field<unsigned char>(surfaceMaskCellFieldId);
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  auto &cellParticleCenterField =
      *structure->template field<Point<double, D>>(cellParticleCenterFieldId);

  structure->iterateCells([&](size_t cellId) {
    cellParticleCenterField[cellId] = structure->cellCenterPosition(cellId);
    return;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID &&
        surfaceField[cellId]) {
      bool emptyCell = true;
      size_t count = 0;
      Point<double, D> centroid(0.);
      //   particles->iterateParticles(structure->cellRegion(cellId),
      //                               [&](size_t id, Point<double, D>
      //                               position)
      //                               {
      //                                 count++;
      //                                 emptyCell = false;
      //                                 centroid += position;
      //                               });
      //   THROW(!emptyCell,
      //         "SimulationSteps::computeCellRepresentativeParticle empty
      //         cell");
      //   cellParticleCenterField[cellId] = centroid /
      //   static_cast<double>(count);

      // check if cell touches solid
      auto neighbors = structure->cellNeighbors(cellId);
      auto cellRegion = structure->cellRegion(cellId);
      std::vector<Point<double, D>> points;
      std::vector<size_t> particleIds;
      bool interfaceWithAir = false;
      for (auto neighbor : neighbors) {
        if (neighbor.id >= 0 &&
            cellMaterialField[neighbor.id] != Definitions::Material::FLUID) {
          interfaceWithAir = true;
          break;
        }
      }
      std::cerr << "******************************************\n";
      std::cerr << "CENTROID ON CELL " << cellId << std::endl;
      if (interfaceWithAir)
        particles->iterateClosestParticles(cellParticleCenterField[cellId], 5,
                                           [&](size_t id, Point<double, D> p) {
                                             count++;
                                             emptyCell = false;
                                             centroid += p;
                                             std::cerr << p << std::endl;
                                           });
      else
        particles->iterateParticles(cellRegion,
                                    [&](size_t id, Point<double, D> p) {
                                      count++;
                                      emptyCell = false;
                                      centroid += p;
                                      std::cerr << p << std::endl;
                                    });
      centroid = centroid / static_cast<double>(count);
      cellParticleCenterField[cellId] = cellRegion.clamp(centroid);
      std::cerr << cellParticleCenterField[cellId] << std::endl;
      std::cerr << "*******************************************\n";
    }
  });
}

template <int D>
void SimulationSteps::transferMidPointFaceVelocitiesToParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialTypeFieldId,
    const std::vector<size_t> &faceVelocityFieldIds,
    const std::vector<size_t> &particleVelocityFieldIds,
    Interpolant<Point<double, D>> *interpolator) {
  auto &cellMaterialField = *structure->template field<Definitions::Material>(
      cellMaterialTypeFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      // setup stencil
      std::set<size_t> faces;
      auto cellFaces = structure->cellFaces(cellId);
      for (auto f : cellFaces)
        faces.insert(f);
      auto neighbors = structure->cellNeighbors(cellId);
      for (auto n : neighbors) {
        if (n.id < 0 || cellMaterialField[n.id] != Definitions::Material::FLUID)
          continue;
        auto nFaces = structure->cellFaces(n.id);
        for (auto f : nFaces)
          faces.insert(f);
      }
      std::vector<Point<double, D>> points;
      for (auto f : faces)
        points.emplace_back(structure->faceCenterPosition(f));
      std::vector<Point<double, D>> values(points.size());
      for (int d = 0; d < D; d++) {
        size_t k = 0;
        for (auto f : faces)
          values[k++][d] = (*velocityFields[d])[f];
      }
      particles->iterateParticles(
          structure->cellRegion(cellId), [&](size_t id, Point<double, D> p) {
            auto value =
                interpolator->interpolateAt((*particles)[id], points, values);
            for (int d = 0; d < D; d++)
              particles->setScalarProperty(particleVelocityFieldIds[d], id,
                                           value[d]);
          });
    }
  });
}

template <int D>
void SimulationSteps::transferMidPointFaceVelocitiesToParticles(
    StructureInterface<D> *structure, ParticleSystem<double, D> *particles,
    size_t cellMaterialTypeFieldId, size_t faceMaterialTypeFieldId,
    const std::vector<size_t> &faceVelocityFieldIds,
    const std::vector<size_t> &particleVelocityFieldIds,
    Definitions::BoundaryVelocity boundaryVelocityCondition,
    Interpolant<Point<double, D>> *interpolator) {
  auto &cellMaterialField = *structure->template field<Definitions::Material>(
      cellMaterialTypeFieldId);
  auto &faceMaterialField = *structure->template field<Definitions::Material>(
      faceMaterialTypeFieldId);
  std::vector<typename StructureInterface<D>::template Field<double> *>
      velocityFields(D, nullptr);
  for (size_t d = 0; d < D; d++)
    velocityFields[d] =
        structure->template field<double>(faceVelocityFieldIds[d]);

  if (boundaryVelocityCondition == Definitions::BoundaryVelocity::NO_SLIP) {
    structure->iterateCells([&](size_t cellId) {
      if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
        // setup stencil
        std::set<size_t> faces;
        auto cellFaces = structure->cellFaces(cellId);
        for (auto f : cellFaces)
          faces.insert(f);
        auto neighbors = structure->cellNeighbors(cellId);
        for (auto n : neighbors) {
          if (n.id < 0 ||
              cellMaterialField[n.id] != Definitions::Material::FLUID)
            continue;
          auto nFaces = structure->cellFaces(n.id);
          for (auto f : nFaces)
            faces.insert(f);
        }
        std::vector<Point<double, D>> points;
        for (auto f : faces)
          points.emplace_back(structure->faceCenterPosition(f));
        std::vector<Point<double, D>> values(points.size());
        for (int d = 0; d < D; d++) {
          size_t k = 0;
          for (auto f : faces)
            values[k++][d] = (*velocityFields[d])[f];
        }
        particles->iterateParticles(
            structure->cellRegion(cellId), [&](size_t id, Point<double, D> p) {
              auto value =
                  interpolator->interpolateAt((*particles)[id], points, values);
              for (int d = 0; d < D; d++)
                particles->setScalarProperty(particleVelocityFieldIds[d], id,
                                             value[d]);
            });
      }
    });
  } else if (boundaryVelocityCondition ==
             Definitions::BoundaryVelocity::FREE_SLIP) {
    structure->iterateCells([&](size_t cellId) {
      if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
        // setup stencil
        std::set<size_t> faces;
        auto cellFaces = structure->cellFaces(cellId);
        for (auto f : cellFaces) {
          faces.insert(f);
          if (faceMaterialField[f] == Definitions::Material::SOLID) {
            // If face is a solid, replace its velocity to satisfy free slip
            // condition: average from the other faces that are in opposite
            // orientations
            auto fOrientation = structure->faceOrientation(f);
            int normalId = 0;
            if (fOrientation == Definitions::Orientation::HORIZONTAL)
              normalId = 1;
            else if (fOrientation == Definitions::Orientation::DEPTH)
              normalId = 2;
            Vector<double, D> velocity;
            int count = 0;
            for (auto otherF : cellFaces) {
              if (structure->faceOrientation(otherF) != fOrientation) {
                count++;
                for (int d = 0; d < D; d++) {
                  if (d != normalId)
                    velocity[d] += (*velocityFields[d])[otherF];
                }
              }
            }
            for (int d = 0; d < D; d++)
              (*velocityFields[d])[f] = velocity[d] / count;
          }
        }
        auto neighbors = structure->cellNeighbors(cellId);
        for (auto n : neighbors) {
          if (n.id < 0 ||
              cellMaterialField[n.id] != Definitions::Material::FLUID)
            continue;
          auto nFaces = structure->cellFaces(n.id);
          for (auto f : nFaces)
            faces.insert(f);
        }
        std::vector<Point<double, D>> points;
        for (auto f : faces)
          points.emplace_back(structure->faceCenterPosition(f));
        std::vector<Point<double, D>> values(points.size());
        for (int d = 0; d < D; d++) {
          size_t k = 0;
          for (auto f : faces)
            values[k++][d] = (*velocityFields[d])[f];
        }
        particles->iterateParticles(
            structure->cellRegion(cellId), [&](size_t id, Point<double, D> p) {
              auto value =
                  interpolator->interpolateAt((*particles)[id], points, values);
              for (int d = 0; d < D; d++)
                particles->setScalarProperty(particleVelocityFieldIds[d], id,
                                             value[d]);
            });
      }
    });
  }
}

template <>
inline void SimulationSteps::initImplicitScene(
    const std::vector<int> &fluidBoxes,
    const std::vector<std::shared_ptr<SolidInterface<2>>> &solids,
    size_t maxLevel, size_t minLevel, size_t cellDistanceFieldId,
    StructureInterface<2> *structure, ParticleSystem<double, 2> *particles,
    size_t cellMaterialFieldId, size_t cellSurfaceMaskFieldId, bool graded) {
  int D = 2;
  std::vector<BBox<double, 2>> fBoxes;
  size_t n = (fluidBoxes.size() - D) / (D * 2);
  Vector<double, 2> step;
  for (int d = 0; d < D; d++)
    step[d] = 1. / (1 << fluidBoxes[d]);
  for (size_t i = 0; i < n; i++) {
    Point<double, 2> pmin, pmax;
    for (int d = 0; d < D; d++) {
      pmin[d] = fluidBoxes[D + i * D * 2 + d] * step[d];
      pmax[d] = fluidBoxes[D + i * D * 2 + D + d] * step[d];
    }
    fBoxes.emplace_back(BBox<double, 2>(pmin, pmax));
    std::cerr << "setup fluid box " << BBox<double, 2>(pmin, pmax) << std::endl;
  }
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(structure);
  graph->refine(1, [&](CellGraph2::Node node) -> bool {
    if (node.level() >= maxLevel)
      return false;
    for (auto box : fBoxes)
      if (node.region().intersect(box))
        return true;
    return false;
  });
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  graph->iterateCells([&](size_t cellId) {
    cellMaterialField[cellId] = Definitions::Material::AIR;
    for (auto box : fBoxes)
      if (box.contains(structure->cellCenterPosition(cellId)))
        cellMaterialField[cellId] = Definitions::Material::FLUID;
  });

  auto &cellSurfaceMaskField =
      *structure->template field<unsigned char>(cellSurfaceMaskFieldId);
  // mark surface
  std::queue<size_t> surfaceCells;
  graph->iterateCells([&](size_t cellId) {
    cellSurfaceMaskField[cellId] = 0;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      bool isSurface = false;
      graph->iterateCellRing(cellId, [&](size_t nid) {
        if (cellMaterialField[nid] != Definitions::Material::FLUID)
          isSurface = true;
      });
      if (isSurface) {
        cellSurfaceMaskField[cellId] = 1;
        surfaceCells.push(cellId);
      }
    }
  });
  markSolids(structure, cellMaterialFieldId, solids);
  std::cerr << "coarsing\n";
  coarseGraph(structure, cellSurfaceMaskFieldId, cellMaterialFieldId,
              [&](Definitions::Material nodeMaterialType, size_t level) {
                return !(nodeMaterialType == Definitions::Material::FLUID &&
                         level <= minLevel);
              });
  std::cerr << "refining near surface\n";
  try {
    auto newSurface = refineNearSurface(structure, cellSurfaceMaskFieldId,
                                        cellDistanceFieldId,
                                        cellMaterialFieldId, surfaceCells);
  } catch (const char *e) {
    std::cerr << e << std::endl;
    THROW(false, "SimulationSteps::updateTree");
  }
  if (graded)
    makeGraded(structure, cellMaterialFieldId);

  // Iterate over fluid cells to inject particles over them
  Injector2 injector(particles, 1.);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      auto region = structure->cellRegion(cellId);
      injector.setSpacing(region.size()[0] / 2);
      injector.setupBoxShape(region);
    }
  });
}

template <>
inline void SimulationSteps::initImplicitScene(
    const std::vector<int> &fluidBoxes,
    const std::vector<std::shared_ptr<SolidInterface<3>>> &solids,
    size_t maxLevel, size_t minLevel, size_t cellDistanceFieldId,
    StructureInterface<3> *structure, ParticleSystem<double, 3> *particles,
    size_t cellMaterialFieldId, size_t cellSurfaceMaskFieldId, bool graded) {
  int D = 3;
  std::vector<BBox<double, 3>> fBoxes;
  size_t n = (fluidBoxes.size() - D) / (D * 2);
  Vector<double, 3> step;
  for (int d = 0; d < D; d++)
    step[d] = 1. / (1 << fluidBoxes[d]);
  for (size_t i = 0; i < n; i++) {
    Point<double, 3> pmin, pmax;
    for (int d = 0; d < D; d++) {
      pmin[d] = fluidBoxes[D + i * D * 2 + d] * step[d];
      pmax[d] = fluidBoxes[D + i * D * 2 + D + d] * step[d];
    }
    fBoxes.emplace_back(BBox<double, 3>(pmin, pmax));
    std::cerr << "setup fluid box " << BBox<double, 3>(pmin, pmax) << std::endl;
  }
  CellGraph3 *graph = dynamic_cast<CellGraph3 *>(structure);
  graph->refine(1, [&](CellGraph3::Node node) -> bool {
    if (node.level() >= maxLevel)
      return false;
    for (auto box : fBoxes)
      if (node.region().intersect(box))
        return true;
    return false;
  });
  auto &cellMaterialField =
      *structure->template field<Definitions::Material>(cellMaterialFieldId);
  graph->iterateCells([&](size_t cellId) {
    cellMaterialField[cellId] = Definitions::Material::AIR;
    for (auto box : fBoxes)
      if (box.contains(structure->cellCenterPosition(cellId)))
        cellMaterialField[cellId] = Definitions::Material::FLUID;
  });

  auto &cellSurfaceMaskField =
      *structure->template field<unsigned char>(cellSurfaceMaskFieldId);
  // mark surface
  std::queue<size_t> surfaceCells;
  graph->iterateCells([&](size_t cellId) {
    cellSurfaceMaskField[cellId] = 0;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      bool isSurface = false;
      graph->iterateCellRing(cellId, [&](size_t nid) {
        if (cellMaterialField[nid] != Definitions::Material::FLUID)
          isSurface = true;
      });
      if (isSurface) {
        cellSurfaceMaskField[cellId] = 1;
        surfaceCells.push(cellId);
      }
    }
  });
  markSolids(structure, cellMaterialFieldId, solids);
  std::cerr << "coarsing\n";
  coarseGraph(structure, cellSurfaceMaskFieldId, cellMaterialFieldId,
              [&](Definitions::Material nodeMaterialType, size_t level) {
                return !(nodeMaterialType == Definitions::Material::FLUID &&
                         level <= minLevel);
              });
  std::cerr << "refining near surface\n";
  try {
    auto newSurface = refineNearSurface(structure, cellSurfaceMaskFieldId,
                                        cellDistanceFieldId,
                                        cellMaterialFieldId, surfaceCells);
  } catch (const char *e) {
    std::cerr << e << std::endl;
    THROW(false, "SimulationSteps::updateTree");
  }
  if (graded)
    makeGraded(structure, cellMaterialFieldId);

  // Iterate over fluid cells to inject particles over them
  Injector3 injector(particles, 1.);
  structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      auto region = structure->cellRegion(cellId);
      injector.setSpacing(region.size()[0] / 2);
      injector.setupBoxShape(region);
    }
  });
}
