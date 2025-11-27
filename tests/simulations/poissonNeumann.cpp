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
  return cos(target.x()) * cos(target.y());
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
  furoo::Timer timer;

  // Variables for numeric solution of Poisson equation
  Eigen::VectorXd rhs(_structure.cellCount()); // Function values to be solved
  Eigen::VectorXd solution, analyticSolution(_structure.cellCount());
  std::vector<Eigen::Triplet<double>> triplets; // to assemble matrix
  Eigen::SparseMatrix<double> A(
      _structure.cellCount(),
      _structure.cellCount()); // Number of cells + boundary faces
  //    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
  //                             Eigen::Lower |
  //                             Eigen::Upper,Eigen::IncompleteLUT<double>>
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                  Eigen::IncompleteLUT<double>>
      solver; // Biconjugate Gradient solver
  // =====
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
    if (cellIndex(cellId) == 0) {
      triplets.emplace_back(
          Eigen::Triplet<double>(cellIndex(cellId), cellIndex(cellId), 1));
      rhs[cellIndex(cellId)] =
          min_9_1_5_analytic(_structure.cellCenterPosition(cellId));
      analyticSolution[cellIndex(cellId)] =
          min_9_1_5_analytic(_structure.cellCenterPosition(cellId));
      return;
    }
    furoo::DifferentialRBF<2> rbf;
    rbf.setKernel(new furoo::CubicKernel<furoo::Point2d>());
    // rbf.setKernel(new furoo::GaussianKernel<furoo::Point2d>(rbfEps));
    rbf.setBasePolynomial(furoo::Definitions::PolynomialType::QUADRATIC);
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
              cellCenter.x() - cellRegionSize.x(), cellCenter.y()));
          break;
        case furoo::Definitions::Side::RIGHT:
          points.emplace_back(furoo::Point2d(
              cellCenter.x() + cellRegionSize.x(), cellCenter.y()));
          break;
        case furoo::Definitions::Side::TOP:
          points.emplace_back(furoo::Point2d(
              cellCenter.x(), cellCenter.y() + cellRegionSize.y()));
          break;
        case furoo::Definitions::Side::BOTTOM:
          points.emplace_back(furoo::Point2d(
              cellCenter.x(), cellCenter.y() - cellRegionSize.y()));
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
    double polynomialValue = 0.;
    // std::cout << "CellId: " <<cellIndex(cellId)<< '\n';
    if (false&&distance(points[0], points[1]) == distance(points[0], points[2]) &&
        distance(points[0], points[1]) == distance(points[0], points[3]) &&
        distance(points[0], points[1]) == distance(points[0], points[4])) {
      //        if (points.size() == 5) {
      char name[20];
      sprintf(name, "%zu-stencil", points.size());
      weights = rbf.laplacianWeights(distance(points[0], points[1]));
    } else {
      weights = rbf.laplacianWeights(points);
      int i = points.size();
      auto position = _structure.cellCenterPosition(cellId);
      //      polynomialValue += weights[i++] * 1;
      //      polynomialValue += weights[i++] * position.x();
      //      polynomialValue += weights[i++] * position.y();
      // polynomialValue += weights[i++] * position.x() * position.y();
      // polynomialValue += weights[i++] * position.x() * position.x();
      // polynomialValue += weights[i++] * position.y() * position.y();
    }

    // Assemble triplets for sparse matrix to solve Poisson system
    double centerWeight = weights[0];
    size_t i = 1;
    // for (size_t i = 0; i < neighborCells.size(); i++) {
    for (auto n : neighborCells) {
      if (n.id < 0) {
        centerWeight += weights[i];
      } else if (n.id >= 0) {
        triplets.emplace_back(Eigen::Triplet<double>(
            cellIndex(cellId), cellIndex(n.id), weights[i]));
      }
      i++;
    }

    triplets.emplace_back(Eigen::Triplet<double>(
        cellIndex(cellId), cellIndex(cellId), centerWeight));
    rhs[cellIndex(cellId)] = min_9_1_5(_structure.cellCenterPosition(cellId));
    rhs[cellIndex(cellId)] -= polynomialValue;
    analyticSolution[cellIndex(cellId)] =
        min_9_1_5_analytic(_structure.cellCenterPosition(cellId));
  });
  std::cout << "Took " << timer.ellapsedTimeInSeconds() << "seconds."
            << std::endl;

  // Solve
  timer.reset();
  A.setFromTriplets(triplets.begin(), triplets.end());
  std::cout << "Solving Poisson system" << '\n';
  Eigen::setNbThreads(8);
  solver.compute(A);
  //  std::cout << A << '\n';
  //  std::cout << rhs << '\n';
  //  std::cout << "\nAnalytic Solution\n" << analyticSolution << '\n';
  solution = solver.solve(rhs);
  //  std::cout << "\nNumerical Solution\n" << solution << '\n';
  std::cout << "Took " << timer.ellapsedTimeInSeconds() << "seconds."
            << std::endl;

  std::cout << "\nSolution error" << '\n';
  std::cout << (solution - analyticSolution).cwiseAbs().maxCoeff() << '\n';
  // Measure error
  return 0;
}
