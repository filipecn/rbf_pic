#include "adaptive_pic2.h"

#include <common/random.h>
#include <structures/cell_graph.h>
#include <structures/point_z_grid.h>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <blas/rbf.h>
#include <map>

namespace furoo {

AdaptivePic2::AdaptivePic2(BBox2d region, size_t minLevel, size_t surfaceLevel)
    : PicSimulation2(new CellGraph2(region),
                     new ParticleSystem2d(new PointZGrid2d(16)),
                     new EulerTemporalIntegrator2(),
                     nullptr),
      _minLevel(minLevel),
      _surfaceLevel(surfaceLevel) {
  this->setRegion(region);
}

void AdaptivePic2::resetTree(size_t level) {
  CellGraph2
      *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &surfaceField = *_structure->field<unsigned char>(surfaceMaskFieldId);
  std::function<bool(size_t, size_t)>
      comp = [&](size_t cellA, size_t cellB) -> bool {
    return graph->node(cellA).level()
        > graph->node(cellB).level();
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  q.push(1);
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    auto curLevel = graph->node(curCellId).level();
    if (curLevel >= level)
      continue;
    auto children = graph->refine(curCellId);
    for (auto child : children) {
      if (!child)
        continue;
      q.push(child);
    }
  }
  _structure->iterateCells([&](size_t i) { surfaceField[i] = 0; });
  coarseTree([&](Definitions::Material material, size_t curLevel) {
    UNUSED_VARIABLE(material);
    return curLevel > level;
  });
  markCells();
  refineSurface();
  trackSurface();
}

std::set<size_t> AdaptivePic2::trackSurface() {
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &surfaceField =
      *_structure->field<unsigned char>(this->surfaceMaskFieldId);
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(this->cellMaterialTypeFieldId);
  std::set<size_t> surfaceCells;
  _structure->iterateFaces([&](size_t faceId) {
    auto edge = graph->edge(faceId);
    auto nodes = edge.nodes();
    if (!nodes[0]
        && cellMaterialField[nodes[1]] == Definitions::Material::FLUID) {
      surfaceField[nodes[1]] = 1;
      surfaceCells.insert(nodes[1]);
    } else if (!nodes[1]
        && cellMaterialField[nodes[0]] == Definitions::Material::FLUID) {
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

void AdaptivePic2::refineSurface() {
  auto surfaceCells = trackSurface();
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(this->cellMaterialTypeFieldId);
  auto &surfaceField =
      *_structure->field<unsigned char>(this->surfaceMaskFieldId);
  std::function<bool(size_t, size_t)>
      comp = [&](size_t cellA, size_t cellB) -> bool {
    if (cellMaterialField[cellA] == cellMaterialField[cellB])
      return graph->node(cellA).level() > graph->node(cellB).level();
    return cellMaterialField[cellA] != Definitions::Material::FLUID;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  for (auto sc : surfaceCells)
    q.push(sc);
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    auto curLevel = graph->node(curCellId).level();
    if (cellMaterialField[curCellId] == Definitions::Material::FLUID) {
      if (curLevel >= _surfaceLevel)
        continue;
      auto children = graph->refine(curCellId);
      for (auto child : children) {
        if (!child)
          continue;
        q.push(child);
        size_t count = 0;
        this->_particles->iterateParticles(graph->cellRegion(child),
                                           [&count](size_t id,
                                                    Point2d p) {
                                             UNUSED_VARIABLE(id);
                                             UNUSED_VARIABLE(p);
                                             count++;
                                           });
        cellMaterialField[child] = (count) ? Definitions::Material::FLUID
                                           : Definitions::Material::AIR;
        surfaceField[child] = 0;
      }
    }
  }
}

void AdaptivePic2::refineFluid(size_t fluidLevel) {
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(this->cellMaterialTypeFieldId);
  std::function<bool(size_t, size_t)>
      comp = [&](size_t cellA, size_t cellB) -> bool {
    if (cellMaterialField[cellA] == cellMaterialField[cellB])
      return graph->node(cellA).level() > graph->node(cellB).level();
    return cellMaterialField[cellA] != Definitions::Material::FLUID;
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(comp);
  q.push(1);
  while (!q.empty()) {
    auto curCellId = q.top();
    q.pop();
    auto curLevel = graph->node(curCellId).level();
    if (cellMaterialField[curCellId] == Definitions::Material::FLUID) {
      if (curLevel >= fluidLevel)
        continue;
      auto children = graph->refine(curCellId);
      for (auto child : children) {
        if (!child)
          continue;
        q.push(child);
        size_t count = 0;
        this->_particles->iterateParticles(graph->cellRegion(child),
                                           [&count](size_t id,
                                                    Point2d p) {
                                             UNUSED_VARIABLE(id);
                                             UNUSED_VARIABLE(p);
                                             count++;
                                           });
        cellMaterialField[child] = (count) ? Definitions::Material::FLUID
                                           : Definitions::Material::AIR;
      }
    }
  }
}

void AdaptivePic2::coarseTree(std::function<bool(Definitions::Material,
                                                 size_t)> doCoarse) {
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &surfaceField =
      *_structure->field<unsigned char>(this->surfaceMaskFieldId);
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(this->cellMaterialTypeFieldId);
  std::function<bool(size_t, size_t)>
      comp = [&](size_t cellA, size_t cellB) -> bool {
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
      for (auto s : siblings)
        if (surfaceField[s])
          coarse = false;
      auto nodeMaterialType = cellMaterialField[siblings[0]];
      if (coarse) {
        if (doCoarse(nodeMaterialType, graph->node(curCellId).level())) {
          try {
            size_t newNode = graph->coarse(curCellId);
            cellMaterialField[newNode] = nodeMaterialType;
            if (newNode)
              q.push(newNode);
          } catch (std::string e) {
            std::ostringstream errorMessage;
            errorMessage << "AdaptivePic2::coarseTree failed to coarse node "
                         << curCellId << " with region "
                         << graph->cellRegion(curCellId);
            THROW(false, errorMessage.str());
          }
        }
      }
    }
  }
}

void AdaptivePic2::updateTree() {
  CellGraph2 *graph = dynamic_cast<CellGraph2 *>(this->_structure.get());
  auto &surfaceField =
      *_structure->field<unsigned char>(this->surfaceMaskFieldId);
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(this->cellMaterialTypeFieldId);

  std::queue<size_t> surfaceCells;
  // iterate air cells fixing surface cells
  graph->iterateCells([&](size_t id) {
    if (cellMaterialField[id] == Definitions::Material::AIR) {
      // TODO when we allow the existence of injectors, their region must be refined too!
      // We don't need to look for particles on coarse air cells, as the air near fluid is refined and
      // the cfl don't permit any particles to jump into a coarse air cell
      if (graph->node(id).level() < _surfaceLevel)
        return;
      bool isFluidCell = false;
      _particles->iterateParticles(graph->cellRegion(id),
                                   [&](size_t pid, Point2d p) {
                                     isFluidCell = true;
                                     UNUSED_VARIABLE(pid);
                                     UNUSED_VARIABLE(p);
                                   });
      if (!isFluidCell)
        return;
      cellMaterialField[id] = Definitions::Material::FLUID;
      bool amIASurface = false;
      for (auto n : graph->neighbors(id)) {
        if (n.id > 0 && cellMaterialField[n.id] == Definitions::Material::AIR)
          amIASurface = true;
        if (n.id > 0 && surfaceField[n.id]) {
          bool isSurface = false;
          for (auto nn : graph->neighbors(n.id))
            if (nn.id > 0
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
      _particles->iterateParticles(graph->cellRegion(id),
                                   [&](size_t pid, Point2d p) {
                                     isAirCell = false;
                                     UNUSED_VARIABLE(pid);
                                     UNUSED_VARIABLE(p);
                                   });
      if (isAirCell) {
        surfaceField[id] = 0;
        cellMaterialField[id] = Definitions::Material::AIR;
        for (auto n : graph->neighbors(id))
          if (n.id > 0
              && cellMaterialField[n.id] == Definitions::Material::FLUID) {
            surfaceField[n.id] = 1;
            surfaceCells.push(n.id);
          }
      } else {
        surfaceCells.push(id);
      }
    }
  });
  // engorda
  coarseTree([&](Definitions::Material nodeMaterialType, size_t level) {
    return !(nodeMaterialType == Definitions::Material::FLUID
        && level <= _minLevel);
  });
  // refine near surface
  std::queue<size_t> newSurface;
  {
    while (!surfaceCells.empty()) {
      size_t curCellId = surfaceCells.front();
      surfaceCells.pop();
      auto neighbors = graph->neighbors(curCellId);
      bool newSurfaceFlag = true;
      for (auto n : neighbors)
        if (n.id > 0 && graph->node(n.id).isValid()
//            && Definitions::Material::FLUID == cellMaterialField[n.id]
            && graph->node(n.id).level()
                < graph->node(curCellId).level()) {
          auto parentMaterial = cellMaterialField[n.id];
          auto children = graph->refine(n.id);
          for (auto c : children)
            cellMaterialField[c] = parentMaterial;
          if (graph->node(children[0]).level()
              < graph->node(curCellId).level()) {
            surfaceCells.push(curCellId);
            newSurfaceFlag = false;
          }
        }
      if (newSurfaceFlag)
        newSurface.push(curCellId);
    }
  }
// make graded
  {
    std::queue<size_t> fluidCells;
    graph->iterateCells([&](size_t i) {
                          if (cellMaterialField[i] == Definitions::Material::FLUID)
                            fluidCells.push(i);
                        }
    );
    while (!fluidCells.empty()) {
      size_t curCellId = fluidCells.front();
      fluidCells.pop();
      auto neighbors = graph->neighbors(curCellId);
      for (auto n : neighbors)
        if (n.id > 0 && graph->node(n.id).isValid()) {
          if (Definitions::Material::FLUID == cellMaterialField[n.id]) {
            if (graph->node(n.id).level()
                < graph->node(curCellId).level() - 1) {
              auto children = graph->refine(n.id);
              for (auto c : children) {
                fluidCells.push(c);
                cellMaterialField[c] = Definitions::Material::FLUID;
              }
            }
          }
        }
    }
  }
// emagrece
  /* std::queue<size_t> childrenCells;
   std::queue<size_t> newSurface;
   while (!surfaceCells.
       empty()
       ) {
     size_t curCell = surfaceCells.front();
     surfaceCells.
         pop();
     if (!surfaceField[curCell] || !structure->
             node(curCell)
         .
             isValid()
         ) {
       continue;
     }
     bool goBackToLineSir = false; // curCell must be processed again
     structure->
         iterateNeighborCells(curCell,
                              2, [&](
             size_t neighborId
         ) {
           if (cellMaterialField[neighborId] ==
               Definitions::Material::AIR
               && structure
                   ->
                       node(neighborId)
                   .
                       level()
                   < _surfaceLevel) {
             try {
               auto children = structure->refine(neighborId);
               for (
                 auto c :
                   children)
                 cellMaterialField[c] =
                     Definitions::Material::AIR;
             } catch (
                 std::string e
             ) {
               std::ostringstream errorMessage;
               errorMessage
                   << "AdaptivePice2::updateTree failed to refine air node: "
                   << neighborId << " with region "
                   << structure->
                       cellRegion(curCell);
               THROW(false, errorMessage.str());
             }
             goBackToLineSir = true;
           } else if (cellMaterialField[neighborId] ==
               Definitions::Material::FLUID
               && structure
                   ->
                       node(neighborId)
                   .
                       level()
                   < _surfaceLevel) {
             try {
               auto children = structure->refine(neighborId);
               for (
                 auto c :
                   children) {
                 cellMaterialField[c] =
                     Definitions::Material::FLUID;
                 childrenCells.
                     push(c);
               }
             } catch (
                 std::string e
             ) {
               std::ostringstream errorMessage;
               errorMessage
                   << "AdaptivePice2::updateTree failed to refine fluid node: "
                   << neighborId << " with region "
                   << structure->
                       cellRegion(curCell);
               THROW(false, errorMessage.str());
             }
             goBackToLineSir = true;
           }
         });
     if (goBackToLineSir)
       surfaceCells.
           push(curCell);
     else
       newSurface.
           push(curCell);
   }*/
// fix surface
  while (!newSurface.empty()) {
    size_t curCell = newSurface.front();
    newSurface.pop();
    /* LEFT BOTTOM RIGHT TOP*/
    int neighborsId[4] = {0, 0, 0, 0};
    auto neighbors = graph->cellNeighbors(curCell);
    for (auto n : neighbors) neighborsId[static_cast <size_t>(n.side)] = n.id;
    Definitions::Side sides[4] =
        {Definitions::Side::LEFT, Definitions::Side::BOTTOM,
         Definitions::Side::RIGHT, Definitions::Side::TOP};
    Definitions::Material materials[4];
    bool allFluid = true;
    for (size_t i = 0; i < 4; i++) {
      materials[i] = (neighborsId[static_cast <size_t>(sides[i])] > 0)
                     ? cellMaterialField[neighborsId[static_cast <size_t>(sides[i])]]
                     : Definitions::Material::AIR;
      if (materials[i] != Definitions::Material::FLUID)
        allFluid = false;
    }
    if (allFluid)
      surfaceField[curCell] = 0;
    for (size_t i = 0; i < 4; i++) {
      if (sides[i] != Definitions::oppositeSide(sides[(i + 1) % 4])) {
        if (materials[i] == Definitions::Material::FLUID
            && materials[(i + 1) % 4] != Definitions::Material::FLUID) {
          surfaceField[neighborsId[static_cast <size_t>(sides[i])]] = 1;
        } else if (materials[i] != Definitions::Material::FLUID
            && materials[(i + 1) % 4]
                == Definitions::Material::FLUID)
          surfaceField[neighborsId[static_cast <size_t>(sides[(i + 1) % 4])]] =
              1;
      }
    }
  }



/*
std::cerr << "Re-refining ineer cells...";
while (!childrenCells.empty()) {
  size_t currCell = childrenCells.front();
  size_t currLevel = structure->node(currCell).level();
  childrenCells.pop();
  auto neighbors = structure->cellNeighbors(currCell);
  for (auto n: neighbors) {
    if (n.id >= 0 && structure->node(n.id).level() < currLevel - 1) {
      auto children = structure->refine(n.id);
      std::cerr << "refined ";
      for (auto c:children) {
        cellMaterialField[c] = Definitions::Material::FLUID;
        childrenCells.push(c);
      }
    }
  }
}
std::cerr << "Done\n";*/
}

void AdaptivePic2::buildTree() {
  markCells();
  refineFluid(_surfaceLevel);
  trackSurface();
}

void AdaptivePic2::solvePressure(double dt) {
  DifferentialRBF<2> rbfGaussian;
  DifferentialRBF<2> rbfQuintic;
  rbfGaussian.setKernel(new GaussianKernel<Point2d>(1.));
  rbfQuintic.setKernel(new QuinticKernel<Point2d>());
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  size_t fluidCount = 0;
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      fluidCount++;
  });
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(fluidCount + 1);
  std::vector<Eigen::Triplet<double>> triplets;
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
  auto &divergenceField = *_structure->field<double>(divergenceFieldId);
  divergenceField.setAll(0.);
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] != Definitions::Material::FLUID)
      return;
    auto cellRegionSize = _structure->cellRegion(cellId).size();
    auto center = _structure->cellCenterPosition(cellId);
    auto neighbors = _structure->cellNeighbors(cellId);
    auto faces = _structure->cellFaces(cellId);
    { // DIV
      double div = 0.;
      // compute divergence
      { // X direction
        std::vector<Point2d> points;
        std::vector<double> values;
        for (auto f : faces)
          if (_structure->faceOrientation(f)
              == Definitions::Orientation::HORIZONTAL) {
            values.emplace_back(velocityField[f]);
            points.emplace_back(_structure->faceCenterPosition(f));
          }
        div += rbfQuintic.gradientAt(center, 1, points, values);
      }
      { // Y direction
        std::vector<Point2d> points;
        std::vector<double> values;
        for (auto f : faces)
          if (_structure->faceOrientation(f)
              == Definitions::Orientation::VERTICAL) {
            values.emplace_back(velocityField[f]);
            points.emplace_back(_structure->faceCenterPosition(f));
          }
        div += rbfQuintic.gradientAt(center, 0, points, values);
      }
      divergenceField[cellId] = div;
      rhs[cellIndex(cellId)] = div / dt;
    }
    {
      std::vector<Point2d> points(neighbors.size() + 1);
      points[0] = center;
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0)
          points[i + 1] = _structure->cellCenterPosition(neighbors[i].id);
        else {
          switch (neighbors[i].side) {
            case Definitions::Side::LEFT:
              points[i + 1] = Point2d(center.x() - cellRegionSize.x(),
                                      center.y());
              break;
            case Definitions::Side::RIGHT:
              points[i + 1] = Point2d(center.x() + cellRegionSize.x(),
                                      center.y());
              break;
            case Definitions::Side::BOTTOM:
              points[i + 1] = Point2d(center.x(),
                                      center.y() - cellRegionSize.y());
              break;
            case Definitions::Side::TOP:
              points[i + 1] = Point2d(center.x(),
                                      center.y() + cellRegionSize.y());
              break;
            default: THROW(false, "PicSimulation2: invalid neighbor side");
          }
        }
      }
      std::vector<double> weights = rbfGaussian.laplacianWeights(points);
      double centerWeight = weights[0];
      double airWeight = 0.;
      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0
            && cellMaterialField[neighbors[i].id]
                == Definitions::Material::AIR) {
          airWeight += weights[i + 1];
          continue;
        }
        if (neighbors[i].id >= 0)
          triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
                                                       cellIndex(neighbors[i].id),
                                                       weights[i + 1]));
        else
          centerWeight += weights[i + 1];
      }
      triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
                                                   cellIndex(cellId),
                                                   centerWeight));
      triplets.emplace_back(Eigen::Triplet<double>(cellIndex(cellId),
                                                   fluidCount,
                                                   airWeight));
    }
  });
  triplets.emplace_back(Eigen::Triplet<double>(indexMap.size(),
                                               indexMap.size(),
                                               1.));
  Eigen::SparseMatrix<double> A(indexMap.size() + 1, indexMap.size() + 1);
  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver;
  solver.compute(A);
  Eigen::VectorXd pressure = solver.solve(rhs);
  auto &faceMaterialField =
      *_structure->field<Definitions::Material>(faceMaterialTypeFieldId);
  auto &pressureField = *_structure->field<double>(pressureFieldId);
  pressureField.setAll(0.);
  _structure->iterateCells([&](size_t cellId) {
    if (cellMaterialField[cellId] == Definitions::Material::FLUID)
      pressureField[cellId] = pressure[cellIndex(cellId)];
  });
  auto &gradientField = *_structure->field<double>(gradientFieldId);
  gradientField.setAll(0.);
  CellGraph2
      *structure = dynamic_cast<CellGraph2 *>(this->_structure.get());
  // correct velocities
  _structure->iterateFaces([&](size_t faceId) {
    if (faceMaterialField[faceId] != Definitions::Material::FLUID)
      return;
    auto cells = _structure->faceCells(faceId);
    double grad = 0.;
    if (structure->node(cells[0].id).level()
        == structure->node(cells[1].id).level()) {
      // finite differences on same level neighbors
      if ((cells[0].side == Definitions::Side::RIGHT
          || cells[0].side == Definitions::Side::TOP)
          && (cells[1].side == Definitions::Side::LEFT
              || cells[1].side == Definitions::Side::BOTTOM))
        std::swap(cells[0], cells[1]);
      double pressureValues[2] = {0., 0.};
      for (size_t i = 0; i < 2; i++)
        if (cellMaterialField[cells[i].id] == Definitions::Material::FLUID)
          pressureValues[i] = pressure[cellIndex(cells[i].id)];
      grad = (pressureValues[1]
          - pressureValues[0])
          / _structure->cellRegion(cells[0].id).size().x();
    } else {
      std::set<size_t> neighborCells;
      for (auto c : cells)
        if (c.id >= 0)
          _structure->iterateNeighborCells(c.id,
                                           1,
                                           [&](size_t i) {
                                             neighborCells.insert(i);
                                           });
      std::vector<double> pressureValues;
      std::vector<Point2d> points;
      for (auto &it : neighborCells) {
        points.emplace_back(_structure->cellCenterPosition(it));
        if (cellMaterialField[it] == Definitions::Material::FLUID)
          pressureValues.emplace_back(pressure[cellIndex(it)]);
        else pressureValues.emplace_back(0.);
      }
      switch (_structure->faceOrientation(faceId)) {
        case Definitions::Orientation::VERTICAL:
          grad =
              rbfQuintic.gradientAt(_structure->faceCenterPosition(faceId),
                                    0,
                                    points,
                                    pressureValues);
          break;
        case Definitions::Orientation::HORIZONTAL:
          grad =
              rbfQuintic.gradientAt(_structure->faceCenterPosition(faceId),
                                    1,
                                    points,
                                    pressureValues);
          break;
        default: THROW(false,
                       "PicSimulation2::computePressure wrong face orientation");
      }
    }
    velocityField[faceId] -= grad * dt;
    gradientField[faceId] = grad;
  });
}

double AdaptivePic2::cfl() const {
  double maxVelocity = 0;
  auto &velocityField = *_structure->field<double>(faceVelocityFieldId);
  _structure->iterateFaces([&](size_t id) {
    maxVelocity = std::max(maxVelocity, std::fabs(velocityField[id]));
  });
  if (maxVelocity == 0.)
    return INFINITY;
  return (_structure->domainRegion().size().x() / (1 << _surfaceLevel))
      / maxVelocity;
}

void AdaptivePic2::reseedParticles() {
  auto &cellMaterialField =
      *_structure->field<Definitions::Material>(cellMaterialTypeFieldId);
  auto
      &surfaceField = *_structure->field<unsigned char>(surfaceMaskFieldId);
  BoxSampler sampler;
  std::vector<size_t> toDelete;
  std::vector<Point2d> toCreate;
  _structure->iterateCells([&](size_t id) {
    if (cellMaterialField[id] != Definitions::Material::FLUID
        || surfaceField[id])
      return;
    auto region = _structure->cellRegion(id);
    std::vector<size_t> cellParticles;
    _particles->iterateParticles(region,
                                 [&](size_t i, Point2d p) {
                                   cellParticles.emplace_back(i);
                                   UNUSED_VARIABLE(p);
                                 });
    if (cellParticles.size() > 25 || cellParticles.size() < 3) {
      for (auto p : cellParticles)
        toDelete.emplace_back(p);

//      auto cellLower = _structure->cellRegion(id).lower();
//      auto cellSize = _structure->cellRegion(id).size().x();
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.5 * cellSize,
//                                    cellLower.y() + cellSize / 6));
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.5 * cellSize,
//                                    cellLower.y() + 5 * cellSize / 6));
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.25 * cellSize,
//                                    cellLower.y() + cellSize / 3));
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.25 * cellSize,
//                                    cellLower.y() + 2 * cellSize / 3));
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.75 * cellSize,
//                                    cellLower.y() + cellSize / 3));
//      toCreate.emplace_back(Point2d(cellLower.x() + 0.75 * cellSize,
//                                    cellLower.y() + 2 * cellSize / 3));
      for (size_t i = 0; i < 6; i++)
        toCreate.emplace_back(sampler.sample(region));
    }
  });
  std::cerr << "BEFORE DELETE " << _particles->size() << std::endl;
  for (auto p : toDelete)
    _particles->remove(p);
  std::cerr << "AFTER DELETE " << _particles->size() << std::endl;
  for (auto p : toCreate)
    _particles->addParticle(p);
  std::cerr << "AFTER CREATE " << _particles->size() << std::endl;
}

}
