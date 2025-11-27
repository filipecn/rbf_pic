#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <blas/rbf.h>
#include <cmath>
#include <common/debug.h>
#include <common/definitions.h>
#include <furoo.h>
#include <geometry/point.h>
#include <structures/cell_graph.h>
#include <utility>

/// Lu = f
/// u(x,y) = e^{-x-y}
/// Dirichlet boundary conditions
double min_9_1_1_analytic(furoo::Point2d target) {
  return exp(-target.x() - target.y());
}

double min_9_1_1(furoo::Point2d target) {
  return 2 * exp(-target.x() - target.y());
}

/// Lu = f
/// u(x,y) = cos(x) cos(y) - 1
/// Homogeneous Neumann boundary conditions
double min_9_1_5_analytic(furoo::Point2d target) {
  return cos(target.x()) * cos(target.y()) - 1;
}

double min_9_1_5(furoo::Point2d target) {
  return -2 * cos(target.x()) * cos(target.y());
}

/****************************************************/
/****************************************************/
int main(int argc, char const *argv[]) {
  if (argc < 4) {
    std::cout << argv[0] << "[0]Regular|[1]Random Refinement rbfEps"
              << std::endl;
    return -1;
  }
  int refinement, problemType;
  double rbfEps;
  sscanf(argv[1], "%d", &problemType);
  sscanf(argv[2], "%d", &refinement);
  sscanf(argv[3], "%lf", &rbfEps);

  furoo::CellGraph2 _structure(
      furoo::BBox2d(furoo::Point2d(0, 0), furoo::Point2d(PI, PI)));
  if (problemType == 1) {
    _structure.refine(1);
    _structure.refine(2);
    _structure.refine(1);
    _structure.refine(2);
    _structure.refine(1);
    _structure.refine(3);
    _structure.refine(19);
    _structure.refine(22);
    _structure.refine(25);
    _structure.refine(5);
    _structure.refine(30);
    _structure.refine(33);
    _structure.refine(36);
    _structure.refine(4);
    _structure.refine(43);
    _structure.refine(46);
    _structure.refine(49);
    for (int k = 0; k < refinement; ++k) {
      _structure.refineAll();
    }
  } else {
    _structure.refine(1, [&](const furoo::CellGraph2::Node &node) -> bool {
      return node.level() < refinement;
    });
  }

  size_t boundaryConditionsFieldId;
  size_t cellValueFieldId;
  size_t faceValueFieldId;
  size_t gradientFieldId;
  std::vector<size_t> _boundaryFacesId;

  // Structure initialization
  cellValueFieldId = _structure.addCellField<double>(0.0);
  _structure.field<double>(cellValueFieldId)
      ->increaseSize(_structure.cellCount());
  auto &cellValueField = *_structure.field<double>(cellValueFieldId);

  faceValueFieldId = _structure.addFaceField<double>(0.0);
  _structure.field<double>(faceValueFieldId)
      ->increaseSize(_structure.faceCount());
  auto &faceValueField = *_structure.field<double>(faceValueFieldId);

  gradientFieldId = _structure.addFaceField<double>(0.0);
  _structure.field<double>(gradientFieldId)
      ->increaseSize(_structure.faceCount());
  auto &gradientField = *_structure.field<double>(gradientFieldId);

  furoo::Timer timer;

  // Variables for numeric solution of Poisson equation
  Eigen::VectorXd rhs(
      _structure.cellCount() +
      _structure.neighbors(0).size()); // Function values to be solved
  Eigen::VectorXd solution, analyticSolution(_structure.cellCount());
  std::vector<Eigen::Triplet<double>> triplets; // to assemble matrix
  Eigen::SparseMatrix<double> A(
      _structure.cellCount() + _structure.neighbors(0).size(),
      _structure.cellCount() +
          _structure.neighbors(0).size()); // Number of cells + boundary faces
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
      solver; // Biconjugate Gradient solver
  std::map<std::pair<size_t, furoo::Definitions::Side>, size_t>
      boundaryFaceMap; // Map cell's boundary face to correct rhs place. Used
                       // later when assemblying matrix with respect Dirichet
                       // boundary conditions
  {
    // =====
    // Defining boundaryConditions
    size_t boundaryFaceCounter = _structure.cellCount();
    _structure.iterateFaces([&](size_t faceId) {
      std::vector<furoo::Definitions::Neighbor> cells =
          _structure.faceCells(faceId);
      bool boundary = false;
      if (cells[0].id <= 0 || cells[1].id <= 0 || cells.size() < 2) {
        // The face is boundary
        size_t innerCell;
        furoo::Definitions::Side innerSide;
        if (cells.size() < 2) {
          innerCell = cells[0].id;
          innerSide = cells[0].side;
        } else if (cells[0].id <= 0) {
          innerCell = cells[1].id;
          innerSide = cells[1].side;
        } else if (cells[1].id <= 0) {
          innerCell = cells[0].id;
          innerSide = cells[0].side;
        }
        rhs[boundaryFaceCounter] =
            min_9_1_1_analytic(_structure.faceCenterPosition(faceId));
        boundaryFaceMap[std::make_pair(innerCell, innerSide)] =
            boundaryFaceCounter++;
      }
    });

    // Compute RBF weights for stencil (including face boundaries)
    // and assemble matrix
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
    timer.reset();
    std::cout << "Assembling weights" << '\n';
    _structure.iterateCells([&](size_t cellId) {
      furoo::DifferentialRBF<2> rbf;
      // rbf.setKernel(new furoo::GaussianKernel<furoo::Point2d>(rbfEps));
      rbf.setKernel(new furoo::QuinticKernel<furoo::Point2d>());
      std::vector<furoo::Point2d> points;
      auto neighborCells = _structure.cellNeighbors(cellId);
      auto cellRegionSize = _structure.cellRegion(cellId).size();
      auto cellCenter = _structure.cellCenterPosition(cellId);

      points.emplace_back(cellCenter);
      for (auto n : neighborCells) {
        if (n.id > 0) {
          points.emplace_back(_structure.cellCenterPosition(n.id));
        } else {
          // Face position (Dirichlet boundary)
          switch (n.side) {
          case furoo::Definitions::Side::LEFT:
            points.emplace_back(furoo::Point2d(
                cellCenter.x() - cellRegionSize.x() / 2, cellCenter.y()));
            break;
          case furoo::Definitions::Side::RIGHT:
            points.emplace_back(furoo::Point2d(
                cellCenter.x() + cellRegionSize.x() / 2, cellCenter.y()));
            break;
          case furoo::Definitions::Side::TOP:
            points.emplace_back(furoo::Point2d(
                cellCenter.x(), cellCenter.y() + cellRegionSize.y() / 2));
            break;
          case furoo::Definitions::Side::BOTTOM:
            points.emplace_back(furoo::Point2d(
                cellCenter.x(), cellCenter.y() - cellRegionSize.y() / 2));
            break;
          default:
            THROW(false,
                  "AdaptiveCenteredPic2::solvePressure invalid neighbor side");
          }
        }
      }

      // Assemble the RBF stencil
      // If we have a regular stencil we can use direct computed weights
      std::vector<double> weights;
      if (distance(points[0], points[1]) == distance(points[0], points[2]) &&
          distance(points[0], points[1]) == distance(points[0], points[3]) &&
          distance(points[0], points[1]) == distance(points[0], points[4])) {
        //        if (points.size() == 5) {
        char name[20];
        sprintf(name, "%zu-stencil", points.size());
        weights = rbf.laplacianWeights(distance(points[0], points[1]));
      } else
        weights = rbf.laplacianWeights(points);

      // Assemble triplets for sparse matrix to solve Poisson system
      double centerWeight = weights[0];
      size_t cellQty = _structure.cellCount();
      for (size_t i = 0; i < neighborCells.size(); i++) {
        if (neighborCells[i].id <= 0) {
          triplets.emplace_back(Eigen::Triplet<double>(
              cellIndex(cellId),
              boundaryFaceMap[std::make_pair(
                  cellId,
                  furoo::Definitions::oppositeSide(neighborCells[i].side))],
              weights[i + 1]));
        } else if (neighborCells[i].id > 0) {
          triplets.emplace_back(Eigen::Triplet<double>(
              cellIndex(cellId), cellIndex(neighborCells[i].id),
              weights[i + 1]));
        }
      }

      triplets.emplace_back(Eigen::Triplet<double>(
          cellIndex(cellId), cellIndex(cellId), centerWeight));
      rhs[cellIndex(cellId)] = min_9_1_1(_structure.cellCenterPosition(cellId));
      analyticSolution[cellIndex(cellId)] =
          min_9_1_1_analytic(_structure.cellCenterPosition(cellId));
    });
    std::cout << timer.ellapsedTimeInSeconds() << std::endl;
    // Dirichlet Identity part
    for (int j = _structure.cellCount();
         j < _structure.cellCount() + _structure.neighbors(0).size(); ++j) {
      triplets.emplace_back(j, j, 1);
    }

    // Solve
    timer.reset();
    A.setFromTriplets(triplets.begin(), triplets.end());
    std::cout << "Solving Poisson system" << '\n';
    Eigen::setNbThreads(8);
    solver.compute(A);
       std::cout << A;
       std::cout << rhs;
    solution = solver.solve(rhs);
    std::cout << solution << '\n';
    std::cout << "analyticSolution\n" << analyticSolution << '\n';
    std::cout << timer.ellapsedTimeInSeconds() << std::endl;
  }

  std::cout << "\nSolution error" << '\n';
  std::cout << (solution.head(_structure.cellCount()) - analyticSolution)
                   .cwiseAbs()
                   .maxCoeff()
            << '\n';
  // Measure error
  return 0;
}
