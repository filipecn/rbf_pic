#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <furoo.h>

namespace furoo {

AdaptiveCenteredPic2::AdaptiveCenteredPic2(BBox2d region, size_t minLevel,
                                           size_t surfaceLevel) {
  this->_structure.reset(new CellGraph2(region));
  this->_particles.reset(new ParticleSystem2d(
      new PointZGrid2d(16),
      //                           new ShepardInterpolant<2>()));
      new MLS<2>(Definitions::PolynomialType::QUADRATIC)));
  this->_integrator.reset(new EulerTemporalIntegrator2());
  this->_minLevel = minLevel;
  this->_surfaceLevel = surfaceLevel;
  this->setRegion(region);

  // FACE BOUNDARY TYPE 0
  faceBoundaryTypeFieldId = _structure->addFaceField<Definitions::Boundary>(
      Definitions::Boundary::NONE);
  _structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId)
      ->increaseSize(_structure->faceCount());

  // FACE MATERIAL TYPE 1
  faceMaterialTypeFieldId = _structure->addFaceField<Definitions::Material>(
      Definitions::Material::AIR);
  _structure->field<Definitions::Material>(faceMaterialTypeFieldId)
      ->increaseSize(_structure->faceCount());

  // CELL MATERIAL TYPE 2
  cellMaterialTypeFieldId = _structure->addCellField<Definitions::Material>(
      Definitions::Material::AIR);
  _structure->field<Definitions::Material>(cellMaterialTypeFieldId)
      ->increaseSize(_structure->cellCount() + 1);

  // CELL SURFACE MASK 3
  surfaceMaskFieldId = _structure->addCellField<unsigned char>(0);
  _structure->field<unsigned char>(surfaceMaskFieldId)
      ->increaseSize(_structure->cellCount() + 1);

  // CELL PRESSURE FIELD 4
  pressureFieldId = _structure->addCellField<double>(0.);
  _structure->field<double>(pressureFieldId)
      ->increaseSize(_structure->cellCount() + 1);

  // CELL DIVERGENCE FIELD 5
  divergenceFieldId = _structure->addCellField<double>(0.);
  _structure->field<double>(divergenceFieldId)
      ->increaseSize(_structure->cellCount() + 1);

  // CELL GRADIENT FIELD 6 - 7
  for (size_t &cellGradientFieldId : _cellGradientFieldIds) {
    cellGradientFieldId = _structure->addCellField<double>(0.);
    _structure->field<double>(cellGradientFieldId)
        ->increaseSize(_structure->cellCount() + 1);
  }

  // CELL VELOCITY FIELD 8 - 9
  for (size_t &cellVelocityFieldId : _cellVelocityFieldIds) {
    cellVelocityFieldId = _structure->addCellField<double>(0.);
    _structure->field<double>(cellVelocityFieldId)
        ->increaseSize(_structure->cellCount() + 1);
  }
  for (size_t &particleVelocityField : particleVelocityFieldIds)
    particleVelocityField = _particles->addScalarProperty();
  for (auto &_domainBoundary : _domainBoundaries)
    _domainBoundary = Definitions::Boundary::NEUMANN;

  // FACE VELOCITY FIELD 10
  faceVelocityFieldId = _structure->addFaceField<double>(0.);
  _structure->field<double>(faceVelocityFieldId)
      ->increaseSize(_structure->faceCount());

  // FACE PRESSURE FIELD 11
  _facePressureFieldId = _structure->addFaceField<double>(0.);
  _structure->field<double>(_facePressureFieldId)
      ->increaseSize(_structure->faceCount());
}

void AdaptiveCenteredPic2::solvePressure(double dt) {
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  auto &velocityFieldX = *_structure->field<double>(_cellVelocityFieldIds[0]);
  auto &velocityFieldY = *_structure->field<double>(_cellVelocityFieldIds[1]);
  auto &faceVelocityField = *_structure->field<double>(faceVelocityFieldId);
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  auto &facePressureField = *_structure->field<double>(_facePressureFieldId);
  //  auto &faceVelocityField = *_structure->field<double>(faceVelocityFieldId);
  //  CellGraph2 *structure = dynamic_cast<CellGraph2 *>(_structure.get());
  DifferentialRBF<2> rbf;
  rbf.setKernel(new GaussianKernel<Point2d>(0.1));

  size_t fluidCount = 0;
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID) {
      fluidCount++;
    }
  });

  // Map for fluid cell indexes
  std::map<size_t, size_t> indexMap;
  std::function<size_t(size_t)> cellIndex;
  cellIndex = [&](size_t cellId) -> size_t {
    auto it = indexMap.find(cellId);
    if (it == indexMap.end()) {
      size_t nextIndex = indexMap.size();
      indexMap[cellId] = nextIndex;
      return indexMap.size() - 1;
    }
    return it->second;
  };

  // Assemble Divergence RHS
  Eigen::VectorXd rhs;
  std::vector<Eigen::Triplet<double>> triplets;
  rhs = Eigen::VectorXd::Zero(fluidCount /* + 1*/);
  auto &divergenceField = *_structure->field<double>(divergenceFieldId);
  divergenceField.setAll(0.);
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID) {
      return;
    }
    auto cellRegionSize = _structure->cellRegion(cellId).size();
    auto cellCenter = _structure->cellCenterPosition(cellId);
    auto neighborFaces = _structure->cellFaces(cellId);
    auto neighborCells = _structure->cellNeighbors(cellId);
    {
      std::vector<Point2d> pointsX, pointsY;
      std::vector<double> valuesX, valuesY;
      double div = 0.;
      for (auto face : neighborFaces) {
        bool verticalFace = _structure->faceOrientation(face) ==
                            Definitions::Orientation::VERTICAL;
        if (verticalFace) {
          pointsX.emplace_back(_structure->faceCenterPosition(face));
          valuesX.emplace_back(faceVelocityField[face]);
        } else {
          pointsY.emplace_back(_structure->faceCenterPosition(face));
          valuesY.emplace_back(faceVelocityField[face]);
        }
      }
      // X-direction divergent
      pointsX.emplace_back(cellCenter);
      valuesX.emplace_back(velocityFieldX[cellId]);
      THROW(pointsX.size() == valuesX.size(),
            "AdaptiveCenteredPic2::solvePressure");
      div += rbf.gradientAt(cellCenter, 0, pointsX, valuesX);
      // Y-direction divergent
      pointsY.emplace_back(cellCenter);
      valuesY.emplace_back(velocityFieldY[cellId]);
      THROW(pointsY.size() == valuesY.size(),
            "AdaptiveCenteredPic2::solvePressure");
      //      std::cerr << "divergence for cell " << cellId << std::endl;
      div += rbf.gradientAt(cellCenter, 1, pointsY, valuesY);
      //      std::cerr << " = " << div << std::endl;
      divergenceField[cellId] = div;
      rhs[cellIndex(cellId)] = div / dt;
    }

    // Laplacian weights
    {
      std::vector<Point2d> points;
      points.emplace_back(cellCenter);
      for (auto n : neighborCells) {
        if (n.id >= 0) {
          points.emplace_back(_structure->cellCenterPosition(n.id));
        } else {
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
      }
      std::vector<double> weights;
      // if we have a regular stencil we can use direct computed weights
      if (points.size() == 5) {
        char name[20];
        sprintf(name, "%zu-stencil", points.size());
        CodeProfiler::instance()[name]++;
        weights = rbf.laplacianWeights(distance(points[0], points[1]));
      } else
        weights = rbf.laplacianWeights(points);

      // Assemble triplets for sparse matrix to solve Poisson system
      double centerWeight = weights[0];
      double airWeight = 0.;
      for (size_t i = 0; i < neighborCells.size(); i++) {
        if (neighborCells[i].id >= 0 &&
            cellMaterialField[neighborCells[i].id] ==
                Definitions::Material::AIR) {
          airWeight += weights[i + 1];
          continue;
        }
        if (neighborCells[i].id >= 0)
          triplets.emplace_back(Eigen::Triplet<double>(
              cellIndex(cellId), cellIndex(neighborCells[i].id),
              weights[i + 1]));
        else
          centerWeight += weights[i + 1];
      }
      triplets.emplace_back(Eigen::Triplet<double>(
          cellIndex(cellId), cellIndex(cellId), centerWeight));
      //      triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
      //                                                   fluidCount,
      //                                                   airWeight));
    }
  });

  // Assemble sparse matrix and solve
  //  triplets.emplace_back(Eigen::Triplet<double>(indexMap.size(),
  //                                               indexMap.size(),
  //                                               1.));
  Eigen::SparseMatrix<double> A(indexMap.size() /* + 1*/,
                                indexMap.size() /* + 1*/);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver;
  Eigen::VectorXd pressure;
  auto &pressureField = *_structure->field<double>(pressureFieldId);
  A.setFromTriplets(triplets.begin(), triplets.end());
  solver.compute(A);
  pressure = solver.solve(rhs);
  pressureField.setAll(0.);
  _structure->iterateCells([&](size_t cellId) {
    pressureField[cellId] = 0.;
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      pressureField[cellId] = pressure[cellIndex(cellId)];
  });

  // compute pressure gradient in cell center and interpolate to particles
  /*auto &cellPressureGradientXField =
      *_structure->field<double>(_cellGradientFieldIds[0]);
  auto &cellPressureGradientYField =
      *_structure->field<double>(_cellGradientFieldIds[1]);
  _structure->iterateCells([&](size_t cellId) {
    std::vector<double> pressureValues(1, pressureField[cellId]);
    std::vector<Point2d> points(1, _structure->cellCenterPosition(cellId));
    auto cellNeighbors = _structure->cellNeighbors(cellId);
    for (auto neighbor : cellNeighbors)
      if (neighbor.id >= 0) {
        pressureValues.emplace_back(pressureField[neighbor.id]);
        points.emplace_back(_structure->cellCenterPosition(neighbor.id));
      }
    DifferentialRBF<2> rbf;
    cellPressureGradientXField[cellId] =
        rbf.gradientAt(points[0], 0, points, pressureValues);
    cellPressureGradientYField[cellId] =
        rbf.gradientAt(points[0], 1, points, pressureValues);
  });
  // correct velocities directly on particles
  std::function<bool(size_t)> isValid = [&](size_t id) -> bool {
    UNUSED_VARIABLE(id);
    return true; //faceMaterialField[id] == Definitions::Material::FLUID;
  };
  _particles->iterateParticles([&](size_t id, Point2d p) {
    Vector2d v(_particles->getScalarProperty(particleVelocityFieldIds[0], id),
               _particles->getScalarProperty(particleVelocityFieldIds[1], id));
    _particles->setScalarProperty(
        particleVelocityFieldIds[0], id, v.x() - dt *
            _structure->sampleCellField(_cellGradientFieldIds[0], p, isValid));
    _particles->setScalarProperty(
        particleVelocityFieldIds[1], id, v.y() - dt *
            _structure->sampleCellField(_cellGradientFieldIds[1], p, isValid));
  });*/
  //return;
  // Interpolate pressure first, compute gradient after
  _structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;

    auto cells = _structure->faceCells(faceId);
    std::vector<double> pressureValues;
    std::vector<Point2d> points;
    std::set<size_t> pointSet;
    for (auto cell : cells) {
      auto cellNeighbors = _structure->cellNeighbors(cell.id);
      for (auto neighbor : cellNeighbors) {
        if (neighbor.id >= 0 &&
            cellMaterialField[neighbor.id] == Definitions::Material::FLUID) {
          pointSet.insert(neighbor.id);
        }
      }
      for (auto pointId : pointSet) {
        points.emplace_back(_structure->cellCenterPosition(pointId));
        pressureValues.emplace_back(pressureField[pointId]);
      }
    }
    MLS<2> interpolant;
    double facePressure = interpolant.interpolateAt(
        _structure->faceCenterPosition(faceId), points, pressureValues);
    facePressureField[faceId] = facePressure;
  });

  // Correct velocities through pressure gradient
  _structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    std::vector<Point2d> points;
    std::vector<double> pressureValues;
    auto faceCenter = _structure->faceCenterPosition(faceId);
    DifferentialRBF<2> rbf;
    auto faceCells = _structure->faceCells(faceId);
    for (auto cell : faceCells) {
      points.emplace_back(_structure->cellCenterPosition(cell.id));
      pressureValues.emplace_back(pressureField[cell.id]);
    }
    points.emplace_back(faceCenter);
    pressureValues.emplace_back(facePressureField[faceId]);

    if (_structure->faceOrientation(faceId) ==
        Definitions::Orientation::VERTICAL)
      faceVelocityField[faceId] -=
          dt * rbf.gradientAt(faceCenter, 0, points, pressureValues);
    else
      faceVelocityField[faceId] -=
          dt * rbf.gradientAt(faceCenter, 1, points, pressureValues);
  });
}

void AdaptiveCenteredPic2::computeExternalForces(double dt) {
  auto &velocityYField = *_structure->field<double>(_cellVelocityFieldIds[1]);
  auto &faceVelocityField = *_structure->field<double>(faceVelocityFieldId);
  auto &faceMaterialType =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);

  _structure->iterateCells([&](size_t i) { velocityYField[i] -= 9.81 * dt; });

  _structure->iterateFaces([&](size_t id) {
    if (faceMaterialType[id] != Definitions::Material::FLUID)
      return;
    if (_structure->faceOrientation(id) ==
        Definitions::Orientation::HORIZONTAL) {
      faceVelocityField[id] -= 9.81 * dt;
    }
  });
}

void AdaptiveCenteredPic2::defineBoundaryConditions(Definitions::Side side,
                                                    Vector2d vel) {
  // _boundaryValues[side] = vel;
}

void AdaptiveCenteredPic2::applyBoundaryCondition() {
  // std::vector<size_t> boundaryFaces;

  // Look for boundary faces that correspond to side
  // _structure->iterateFaces([&](size_t id) {
  //   auto cells = _structure->faceCells(id);
  //   if( faceBoundaryTypeField(id) == Definitions::Boundary::DIRICHLET ) {
  //     if( cells[0].id >=0 && cells[0].side == Definitions.oppositeSide(side) ){
  //
  //     }
  //   }
  // });

  // Get faces id and set them up with vel value
  // TODO particles should interpolate x and y velocity from boundary

}

void AdaptiveCenteredPic2::transferFromParticlesToGrid() {
  CellGraph2 *structure = dynamic_cast<CellGraph2 *>(_structure.get());
  auto &velocityXField = *_structure->field<double>(_cellVelocityFieldIds[0]);
  auto &velocityYField = *_structure->field<double>(_cellVelocityFieldIds[1]);
  auto &faceVelocityField = *_structure->field<double>(faceVelocityFieldId);
  auto &faceMaterialType =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  //  auto &faceBoundaryTypeField =
  //      *_structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId);

  // Transfer particle velocities to cell CENTER
  // It is needed because differential operators have higher accuracy where the
  // the value is known
  _structure->iterateCells([&](size_t id) {
    velocityXField[id] = velocityYField[id] = 0.;
    if (cellMaterialField[id] == Definitions::Material::FLUID) {
      // todo: interpolate both fields inside particles to avoid making two
      // identical searches
      try {
        velocityXField[id] = _particles->sampleScalarProperty(
            particleVelocityFieldIds[0], _structure->cellCenterPosition(id),
            20);
        velocityYField[id] = _particles->sampleScalarProperty(
            particleVelocityFieldIds[1], _structure->cellCenterPosition(id),
            20);
      } catch (std::string e) {
        std::cerr << e << "transferFromParticlesToGrid" << std::endl;
        //        throw (std::string(
        //            "transferFromParticlesToGrid: particle to cell centers"));
      }
    }
  });
  TimeLogger timer;
  // Transfer to cells FACES
  size_t count = 0;
  _structure->iterateFaces([&](size_t id) {
    if (faceMaterialType[id] != Definitions::Material::FLUID)
      return;
    {
      auto neighbors = structure->faceCells(id);

      std::vector<Point2d> points;
      std::vector<double> values;
      size_t fieldId =
          _structure->faceOrientation(id) == Definitions::Orientation::VERTICAL
              ? particleVelocityFieldIds[0]
              : particleVelocityFieldIds[1];

      for (auto n : neighbors) {
        timer.log("iterating faces");
        if (n.id >= 0 &&
            cellMaterialField[n.id] == Definitions::Material::FLUID) {
          count++;
          _particles->iterateParticles(
              structure->cellRegion(n.id), [&](size_t pid, Point2d p) {
                values.emplace_back(
                    _particles->getScalarProperty(fieldId, pid));
                points.emplace_back(p);
              });
          timer.log("iterating particles");
        }
      }

      MLS<2> mls(Definitions::PolynomialType::QUADRATIC);
      //      ShepardInterpolant2 mls;
      try {
        faceVelocityField[id] = mls.interpolateAt(
            structure->faceCenterPosition(id), points, values);
      } catch (std::string &e) {
        std::cerr << "ParticleToGrid: attempting with bigger neighborhood"
                  << std::endl;
        std::set<int> neighborSet;

        for (auto cell : neighbors) {
          for (auto n : structure->cellNeighbors(cell.id))
            neighborSet.insert(n.id);
        }
        for (auto n : neighborSet) {
          if (n >= 0 && cellMaterialField[n] == Definitions::Material::FLUID)
            _particles->iterateParticles(
                structure->cellRegion(n), [&](size_t pid, Point2d p) {
                  values.emplace_back(
                      _particles->getScalarProperty(fieldId, pid));
                  points.emplace_back(p);
                });
        }
        try {
          faceVelocityField[id] = mls.interpolateAt(
              structure->faceCenterPosition(id), points, values);
        } catch (std::string &e) {
          std::cerr << structure->faceCenterPosition(id) << "   <<< " << id
                    << std::endl;
          THROW(false, "ParticleToGrid: second attempt failed");
        }
        //        std::cerr << _structure->faceCenterPosition(id) << std::endl;
        //        faceVelocityField[id] = _particles->sampleScalarProperty(
        //            particleVelocityFieldIds[1],
        //            _structure->faceCenterPosition(id), 15);
      }
    }
    timer.log("mls");
    return;
    size_t numberOfSamples = 10;
    //        (faceBoundaryTypeField[id] == Definitions::Boundary::DIRICHLET) ?
    //        30
    //                                                                        :
    //                                                                        15;
    if (_structure->faceOrientation(id) ==
        Definitions::Orientation::HORIZONTAL) {
      try {
        faceVelocityField[id] = _particles->sampleScalarProperty(
            particleVelocityFieldIds[1], _structure->faceCenterPosition(id),
            numberOfSamples);
      } catch (std::string e) {
        std::cerr << e << std::endl;
        //        throw (std::string(
        //            "transferFromParticlesToGrid: particle to horizontal face
        //            centers"));
      }
    } else if (_structure->faceOrientation(id) ==
               Definitions::Orientation::VERTICAL) {
      try {

        faceVelocityField[id] = _particles->sampleScalarProperty(
            particleVelocityFieldIds[0], _structure->faceCenterPosition(id),
            numberOfSamples);
      } catch (std::string e) {
        std::cerr << e << std::endl;
        //        throw (std::string(
        //            "transferFromParticlesToGrid: particle to vertival face
        //            centers"));
      }
    }
  });
  std::cerr << count << " searches, " << count * 0.0045 << std::endl;
  timer.report();
}

void AdaptiveCenteredPic2::transferFromGridToParticles() {
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  std::function<bool(size_t)> isValid = [&](size_t id) -> bool {
    return faceMaterialField[id] == Definitions::Material::FLUID;
  };

  _particles->iterateParticles([&](size_t id, Point2d p) {
    // todo: same optimization of transfer from particles to grid
    _particles->setScalarProperty(
        particleVelocityFieldIds[0], id,
        _structure->sampleFaceField(faceVelocityFieldId, p,
                                    Definitions::Orientation::VERTICAL,
                                    isValid));
    //        _structure->sampleCellField(_cellVelocityFieldIds[0], p,
    //        isValid));
    _particles->setScalarProperty(
        particleVelocityFieldIds[1], id,
        _structure->sampleFaceField(faceVelocityFieldId, p,
                                    Definitions::Orientation::HORIZONTAL,
                                    isValid));
    //        _structure->sampleCellField(_cellVelocityFieldIds[1], p,
    //        isValid));
  });
}

double AdaptiveCenteredPic2::cfl() const {
  double maxVelocity = 0;
  //  auto &velocityXField =
  //  *_structure->field<double>(_cellVelocityFieldIds[0]); auto &velocityYField
  //  = *_structure->field<double>(_cellVelocityFieldIds[1]);
  //  _structure->iterateCells([&](size_t id) {
  //    maxVelocity = std::max(maxVelocity, std::fabs(velocityXField[id]));
  //    maxVelocity = std::max(maxVelocity, std::fabs(velocityYField[id]));
  //  });
  _particles->iterateParticles([&](size_t id, Point2d p) {
    UNUSED_VARIABLE(p);
    maxVelocity = std::max(maxVelocity, std::fabs(_particles->getScalarProperty(
                                            particleVelocityFieldIds[0], id)));
    maxVelocity = std::max(maxVelocity, std::fabs(_particles->getScalarProperty(
                                            particleVelocityFieldIds[1], id)));
  });
  if (maxVelocity == 0.)
    return INFINITY;
  return (_structure->domainRegion().size().x() / (1 << _surfaceLevel)) /
         maxVelocity;
}

void AdaptiveCenteredPic2::loadFromFile(const char *filename) {
//  IO::loadCellGraph(filename, dynamic_cast<CellGraph2 *>(_structure.get()),
//                    cellMaterialTypeFieldId);
  // fix fields
  // FACE BOUNDARY TYPE 0
  _structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId)->clear();
  _structure->field<Definitions::Boundary>(faceBoundaryTypeFieldId)
      ->increaseSize(_structure->faceCount());
  // FACE MATERIAL TYPE 1
  _structure->field<Definitions::Material>(faceMaterialTypeFieldId)->clear();
  _structure->field<Definitions::Material>(faceMaterialTypeFieldId)
      ->increaseSize(_structure->faceCount());
  // CELL MATERIAL TYPE 2
  // CELL SURFACE MASK 3
  _structure->field<unsigned char>(surfaceMaskFieldId)->clear();
  _structure->field<unsigned char>(surfaceMaskFieldId)
      ->increaseSize(_structure->cellCount() + 1);
  // CELL PRESSURE FIELD 4
  _structure->field<double>(pressureFieldId)->clear();
  _structure->field<double>(pressureFieldId)
      ->increaseSize(_structure->cellCount() + 1);
  // CELL DIVERGENCE FIELD 5
  _structure->field<double>(divergenceFieldId)->clear();
  _structure->field<double>(divergenceFieldId)
      ->increaseSize(_structure->cellCount() + 1);
  // CELL GRADIENT FIELD 6 - 7
  for (size_t &cellGradientFieldId : _cellGradientFieldIds) {
    _structure->field<double>(cellGradientFieldId)->clear();
    _structure->field<double>(cellGradientFieldId)
        ->increaseSize(_structure->cellCount() + 1);
  }
  // CELL VELOCITY FIELD 8 - 9
  for (size_t &cellVelocityFieldId : _cellVelocityFieldIds) {
    _structure->field<double>(cellVelocityFieldId)->clear();
    _structure->field<double>(cellVelocityFieldId)
        ->increaseSize(_structure->cellCount() + 1);
  }
  // FACE VELOCITY FIELD 10
  _structure->field<double>(faceVelocityFieldId)->clear();
  _structure->field<double>(faceVelocityFieldId)
      ->increaseSize(_structure->faceCount());
  // FACE PRESSURE FIELD 11
  _structure->field<double>(_facePressureFieldId)->clear();
  _structure->field<double>(_facePressureFieldId)
      ->increaseSize(_structure->faceCount());
  auto &surfaceField = *_structure->field<unsigned char>(surfaceMaskFieldId);
  _structure->iterateCells([&](size_t i) { surfaceField[i] = 0; });
  this->trackSurface();
}

} // namespace furoo
