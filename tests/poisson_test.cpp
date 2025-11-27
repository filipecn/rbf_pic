#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(Poisson, CellGraph) {
  { // ALL DIRICHLET
    CellGraph2 graph(BBox2d::squareBox(), [](const CellGraph2::Node &node) {
      return node.level() < 5;
    });
    std::cerr << graph.cellCount() << std::endl;

    size_t dirichletCount = 32 * 4;
    DifferentialRBF<2> rbf;
    rbf.setKernel(new GaussianKernel<Point2d>(0.25));
    Eigen::VectorXd rhs =
        Eigen::VectorXd::Zero(graph.cellCount() + dirichletCount);
    std::vector<Eigen::Triplet<double>> triplets;
    size_t dirBdId = graph.cellCount();

    graph.iterateCells([&](size_t cellId) {
      std::cerr << "processing " << cellId << std::endl;
      auto cellRegion = graph.cellRegion(cellId);
      auto center = graph.cellCenterPosition(cellId);
      std::cerr << center << std::endl;
      auto neighbors = graph.cellNeighbors(cellId);
      rhs[cellId - 1] =
          -1.25 * SQR(PI) * sin(PI * center.x()) * cos(0.5 * PI * center.y());
      std::vector<Point2d> points(neighbors.size() + 1);
      points[0] = center;
      std::queue<size_t> dirichletIndices;

      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0)
          points[i + 1] = graph.cellCenterPosition(neighbors[i].id);
        else {
          switch (neighbors[i].side) {
          case Definitions::Side::LEFT:
            points[i + 1] = Point2d(cellRegion.lower().x(), center.y());
            break;
          case Definitions::Side::RIGHT:
            points[i + 1] = Point2d(cellRegion.upper().x(), center.y());
            break;
          case Definitions::Side::BOTTOM:
            points[i + 1] = Point2d(center.x(), cellRegion.lower().y());
            rhs[dirBdId] = sin(PI * points[i + 1].x());
            break;
          case Definitions::Side::TOP:
            points[i + 1] = Point2d(center.x(), cellRegion.upper().y());
            break;
          default:
            THROW(false, "PicSimulation2: invalid neighbor side");
          }
          dirichletIndices.push(dirBdId++);
        }
      }

      std::vector<double> weights = rbf.laplacianWeights(points);
      std::cerr << "weights\n";
      for (auto w : weights)
        std::cerr << w << " ";
      std::cerr << std::endl;

      double centerWeight = weights[0];
      //      for (size_t i = 0; i < 3; i++) {
      //        size_t id = graph.cellCount() + dirichletCount + i;
      //        triplets.emplace_back(Eigen::Triplet<double>(cellId - 1, id,
      //                                                     weights[weights.size()
      //                                                         - i]));
      //      }

      for (size_t i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].id >= 0)
          triplets.emplace_back(Eigen::Triplet<double>(
              cellId - 1, neighbors[i].id - 1, weights[i + 1]));
        else { // if (neighbors[i].side == Definitions::Side::BOTTOM) {
          auto ind = dirichletIndices.front();
          dirichletIndices.pop();
          triplets.emplace_back(
              Eigen::Triplet<double>(cellId - 1, ind, weights[i + 1]));
          triplets.emplace_back(Eigen::Triplet<double>(ind, ind, 1.));
        }
      }

      ASSERT_FATAL(dirichletIndices.empty());
      triplets.emplace_back(
          Eigen::Triplet<double>(cellId - 1, cellId - 1, centerWeight));
    });
    //    for (size_t i = 0; i < 3; i++) {
    //      size_t id = graph.cellCount() + dirichletCount + i;
    //      triplets.emplace_back(Eigen::Triplet<double>(id, id, 1.));
    //    }
    std::cerr << dirBdId << std::endl;
    Eigen::SparseMatrix<double> A(graph.cellCount() + dirichletCount,
                                  graph.cellCount() + dirichletCount);
    A.setFromTriplets(triplets.begin(), triplets.end());
    std::cerr << "A\n" << A << std::endl;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>
        solver;

    solver.compute(A);
    std::cerr << "RHS\n" << rhs << std::endl << std::endl;
    Eigen::VectorXd x = solver.solve(rhs);
    std::cerr << x;
    double error = 0.;
    double error2 = 0.;

    graph.iterateCells([&](size_t cellId) {
      auto cp = graph.cellCenterPosition(cellId);
      double sol = sin(PI * cp.x()) * cos(0.5 * PI * cp.y());
      error = std::max(error, std::fabs(sol - x[cellId - 1]));
      error2 += SQR(sol - x[cellId - 1]);
    });

    std::cerr << "maxError " << error << std::endl;
    std::cerr << "Error2 " << std::sqrt(error2) << std::endl;
  }
}

TEST(Poisson, SimQuadTree) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  LinearVector s(qt.cellCount());
  LinearVector b(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearVector x(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearMatrix A(qt.cellCount() + bh.dirichletBoundaryCount());
  for (size_t i = 0; i < qt.cellCount(); i++) {
    double neumanNeighborsSum = 0.;
    static size_t dirBdId = qt.cellCount(); // Dirichlet boundary iterator
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0) {
            if (o == SimulationDomain2::Orientation::LEFT) {
              b[dirBdId] = cos(PI * p.y());
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            } else if (o == SimulationDomain2::Orientation::RIGHT) {
              b[dirBdId] = -cos(PI * p.y());
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            } else
              neumanNeighborsSum += w;
          } else
            A(i, id) = w;
        });
    // A(i, i) += neumanNeighborsSum;
    Point2d cp = qt.cellCenterPosition(i);
    b[i] = -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
    s[i] = cos(PI * cp.x()) * cos(PI * cp.y());
  }
  BiconjugateGradientLinearSolver bicg;
  bicg.solve(&A, &b, &x);
  LinearVector error = x - s;
  EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
  std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
}

TEST(Poisson, RegularGrid) {
  {
    SimRegularDomain2 qt(16);
    DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN,
                              DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                              DomainBoundaryHandler2::BoundaryType::NEUMANN);
    qt.setBoundaryHandler(&bh);
    bh.markFaces(&qt);
    LinearVector s(qt.cellCount());
    LinearVector b(qt.cellCount() + bh.dirichletBoundaryCount());
    LinearVector x(qt.cellCount() + bh.dirichletBoundaryCount());
    LinearMatrix A(qt.cellCount() + bh.dirichletBoundaryCount());
    for (size_t i = 0; i < qt.cellCount(); i++) {
      static size_t dirBdId = qt.cellCount(); // Dirichlet boundary iterator
      qt.iterateStencilAtCellCenter(i, [&](int id, double w, Point2d p,
                                           SimulationDomain2::Orientation o) {
        if (id < 0) {
          if (o == SimulationDomain2::Orientation::LEFT) {
            b[dirBdId] = cos(PI * p.y());
            A(dirBdId, dirBdId) = 1;
            A(i, dirBdId++) = w;
          } else if (o == SimulationDomain2::Orientation::RIGHT) {
            b[dirBdId] = -cos(PI * p.y());
            A(dirBdId, dirBdId) = 1;
            A(i, dirBdId++) = w;
          }
        } else
          A(i, id) = w;
      });
      // A(i, i) += neumanNeighborsSum;
      Point2d cp = qt.cellCenterPosition(i);
      b[i] = -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
      s[i] = cos(PI * cp.x()) * cos(PI * cp.y());
    }
    BiconjugateGradientLinearSolver bicg;
    bicg.solve(&A, &b, &x);
    LinearVector error = x - s;
    EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
    std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
  }
}

TEST(Poisson, SimQuadTreeFAceCount) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  LinearVector s(qt.cellCount());
  LinearVector b(qt.cellCount());
  LinearMatrix A(b.size());
  for (size_t i = 0; i < qt.cellCount(); i++) {
    double neumanNeighborsSum = 0.;
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0) {
            if (o == SimulationDomain2::Orientation::LEFT)
              b[i] -= w * cos(PI * p.y());
            else if (o == SimulationDomain2::Orientation::RIGHT)
              b[i] += w * cos(PI * p.y());
            else
              neumanNeighborsSum += w;
          } else
            A(i, id) = w;
        });
    // A(i, i) += neumanNeighborsSum;
    Point2d cp = qt.cellCenterPosition(i);
    b[i] += -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
    s[i] = cos(PI * cp.x()) * cos(PI * cp.y());
  }
  BiconjugateGradientLinearSolver bicg;
  LinearVector x(qt.cellCount());
  bicg.solve(&A, &b, &x);
  LinearVector error = x - s;
  EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
  std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
}

TEST(Poisson, RegularGridB) {
  SimRegularDomain2 qt(16);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::NEUMANN);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  LinearVector s(qt.cellCount());
  LinearVector b(qt.cellCount());
  LinearMatrix A(b.size());
  for (size_t i = 0; i < qt.cellCount(); i++) {
    double neumanNeighborsSum = 0.;
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0) {
            if (o == SimulationDomain2::Orientation::LEFT)
              b[i] -= w * cos(PI * p.y());
            else if (o == SimulationDomain2::Orientation::RIGHT)
              b[i] += w * cos(PI * p.y());
            else
              neumanNeighborsSum += w;
          } else
            A(i, id) = w;
        });
    // A(i, i) += neumanNeighborsSum;
    Point2d cp = qt.cellCenterPosition(i);
    b[i] += -2. * SQR(PI) * cos(PI * cp.x()) * cos(PI * cp.y());
    s[i] = cos(PI * cp.x()) * cos(PI * cp.y());
  }
  BiconjugateGradientLinearSolver bicg;
  LinearVector x(qt.cellCount());
  bicg.solve(&A, &b, &x);
  LinearVector error = x - s;
  EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
  std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
}

TEST(Poisson, SimQuadTreeDirichlet) {
  SimQuadTree qt(BBox2D::make_unit_bbox(), 4);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  LinearVector s(qt.cellCount());
  LinearVector b(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearVector x(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearMatrix A(qt.cellCount() + bh.dirichletBoundaryCount());
  for (size_t i = 0; i < qt.cellCount(); i++) {
    static size_t dirBdId = qt.cellCount(); // Dirichlet boundary iterator
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0) {
            if (o == SimulationDomain2::Orientation::BOTTOM) {
              b[dirBdId] = sin(PI * p.x());
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            } else {
              b[dirBdId] = 0;
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            }
          } else
            A(i, id) = w;
        });
    // A(i, i) += neumanNeighborsSum;
    Point2d cp = qt.cellCenterPosition(i);
    b[i] = -1.25 * SQR(PI) * sin(PI * cp.x()) * cos(0.5 * PI * cp.y());
    s[i] = sin(PI * cp.x()) * cos(0.5 * PI * cp.y());
  }
  BiconjugateGradientLinearSolver bicg;
  bicg.solve(&A, &b, &x);
  LinearVector error = x - s;
  EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
  std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
}

TEST(Poisson, RegularGridDirichlet) {
  SimRegularDomain2 qt(32);
  DomainBoundaryHandler2 bh(DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET,
                            DomainBoundaryHandler2::BoundaryType::DIRICHLET);
  qt.setBoundaryHandler(&bh);
  bh.markFaces(&qt);
  LinearVector s(qt.cellCount());
  LinearVector b(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearVector x(qt.cellCount() + bh.dirichletBoundaryCount());
  LinearMatrix A(qt.cellCount() + bh.dirichletBoundaryCount());
  for (size_t i = 0; i < qt.cellCount(); i++) {
    static size_t dirBdId = qt.cellCount(); // Dirichlet boundary iterator
    qt.iterateStencilAtCellCenter(
        i, [&](int id, double w, Point2d p, SimulationDomain2::Orientation o) {
          if (id < 0) {
            if (o == SimulationDomain2::Orientation::BOTTOM) {
              b[dirBdId] = sin(PI * p.x());
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            } else {
              b[dirBdId] = 0;
              A(dirBdId, dirBdId) = 1;
              A(i, dirBdId++) = w;
            }
          } else
            A(i, id) = w;
        });
    // A(i, i) += neumanNeighborsSum;
    Point2d cp = qt.cellCenterPosition(i);
    b[i] = -1.25 * SQR(PI) * sin(PI * cp.x()) * cos(0.5 * PI * cp.y());
    s[i] = sin(PI * cp.x()) * cos(0.5 * PI * cp.y());
  }
  // std::cout << A << '\n';
  // std::cout << b << '\n';
  // std::cout << s << '\n';
  BiconjugateGradientLinearSolver bicg;
  bicg.solve(&A, &b, &x);
  LinearVector error = x - s;
  EXPECT_EQ(Blas::norm(&error) / qt.cellCount() < 1e-3, true);
  std::cout << "ERROR " << Blas::norm(&error) / qt.cellCount() << std::endl;
}
