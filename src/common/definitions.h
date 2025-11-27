#ifndef FUROO_DEFINITIONS_H
#define FUROO_DEFINITIONS_H

#include <common/debug.h>
#include <functional>
#include <geometry/point.h>
#include <geometry/vector.h>
#include <iostream>
#include <vector>

namespace furoo {

/// Utility class containing most of the definitions of the library. Should not
/// be instatiated.
class Definitions {
public:
  enum class MeshLocation {
    CELL_CENTER,
    VERTEX_CENTER,
    FACE_CENTER,
    HORIZONTAL_FACE_CENTER,
    VERTICAL_FACE_CENTER,
    DEPTH_FACE_CENTER,
    CUSTOM
  };
  enum class PolynomialType {
    CONSTANT,
    LINEAR,
    QUADRATIC,
    CUBIC,
    ZERO,
    CUSTOM
  };
  enum class Boundary { DIRICHLET = 0, NEUMANN = 1, NONE = 2, CUSTOM = 3 };
  enum class BoundaryVelocity {
    NONE = 0,
    NO_SLIP_HORIZONTAL = 1,
    NO_SLIP_VERTICAL = 2,
    NO_SLIP_DEPTH = 3,
    FREE_SLIP_HORIZONTAL = 4,
    FREE_SLIP_VERTICAL = 5,
    FREE_SLIP_DEPTH = 6,
    NO_SLIP = 7,
    FREE_SLIP = 8,
    CUSTOM = 9
  };
  // TODO: fix this order!
  enum class Side {
    LEFT = 0,
    RIGHT = 1,
    TOP = 2,
    BOTTOM = 3,
    FRONT = 4,
    BACK = 5,
    CUSTOM = 6
  };
  enum class Orientation {
    HORIZONTAL = 0,
    VERTICAL = 1,
    DEPTH = 2,
    CUSTOM = 3,
    ANY = 4,
    XY_PLANE = 5,
    XZ_PLANE = 6,
    YZ_PLANE = 7
  };
  enum class Material { FLUID = 0, AIR = 1, SOLID = 2, CUSTOM = 3 };

  struct Neighbor {
    int id;
    Definitions::Side side;
  };

  template <int D> struct Intersection {
    Point<double, D> point;
    Vector<double, D> normal;
    bool isValid;
  };

  /// Checks if two orientations are equivalent. A vector orientation
  /// is equivalent to the orientation of the plane it
  /// \param a first orientation
  /// \param b second orientation
  /// \return true if equivalent
  static bool isEquivalent(Orientation a, Orientation b);

  /// \param id mesh location integer index
  /// \return mesh location associated to index
  static MeshLocation meshLocationFromId(size_t id);

  /// \param id boundary integer index
  /// \return boundary associated to index
  static Boundary boundaryFromId(size_t id);

  /// Given a numerical value of the material, returns it proper definitions
  /// \param id numerical value of the material
  /// \return Proper definition from Definition class
  static Material materialFromId(size_t id);

  /// Given a numerical value for the side, returns its proper definition
  /// \param id numerical value for side
  /// \return Proper definition for side
  static Side sideFromId(size_t id);

  /// Given two orientations, return the third orientations that is orhogonal to
  /// them
  /// \param axisA first orientation
  /// \param axisA second orientation
  /// \return third orientation orthogonal to the first two parameters given
  static Orientation orthogonalOrientation(Orientation axisA,
                                           Orientation axisB);

  /// Given an orientation, returns the corresponding sides that belong to that
  /// orientation \param orientation orientation which sides are been looked for
  /// \return a vector containing the two sides that belong to orientation
  /// (VERTICAL gives BOTTOM, TOP)
  static std::vector<Side> orientationSides(Orientation orientation);

  /// Given a side, returns the opposite side on the same orientation
  /// \param s side which opposite side is expected
  /// \return opposite side that belong to same orientation as side (ex: LEFT
  /// gives HORIZONTAL)
  static Side oppositeSide(Side s);

  /// Given a side, returns its corresponding orientation
  /// \param s side which orientation are expected
  /// \return Orientation corresponding to side
  static Orientation sideOrientation(Side s);

  friend std::ostream &operator<<(std::ostream &os, Side side) {
    switch (side) {
    case Side::LEFT:
      os << "Side::LEFT";
      break;
    case Side::RIGHT:
      os << "Side::RIGHT";
      break;
    case Side::TOP:
      os << "Side::TOP";
      break;
    case Side::BOTTOM:
      os << "Side::BOTTOM";
      break;
    case Side::FRONT:
      os << "Side::FRONT";
      break;
    case Side::BACK:
      os << "Side::BACK";
      break;
    default:
      os << "Side::CUSTOM";
    }
    return os;
  }

  friend std::ostream &operator<<(std::ostream &os, Orientation orientation) {
    switch (orientation) {
    case Orientation::HORIZONTAL:
      os << "Orientation::HORIZONTAL";
      break;
    case Orientation::VERTICAL:
      os << "Orientation::VERTICAL";
      break;
    case Orientation::DEPTH:
      os << "Orientation::DEPTH";
      break;
    case Orientation::ANY:
      os << "Orientation::ANY";
      break;
    default:
      os << "Orientation::CUSTOM";
    }
    return os;
  }

  /// This function should be used within rbf methods for 2D points. Given a
  /// PolynomialType, this method returns its proper polynomial base, following:
  /// CONSTANT: 1
  /// LINEAR: 1 x y
  /// QUADRATIC: 1 x y xy x^2 y^2
  /// CUBIC: 1 x y xy x^2 y^2 x^2y xy^2 x^3 y^3
  static std::function<std::vector<double>(const double *, size_t)>
  polynomial(PolynomialType baseType);

  /// This function should be used within rbf methods for 2D points. Given a
  /// PolynomialType, this method returns its proper firs derivative polynomial
  /// base, following:
  /// CONSTANT: 0
  /// LINEAR for dx: 0 1 0
  /// LINEAR for dy: 0 0 1
  /// QUADRATIC for dx: 0 1 0 y 2x 0
  /// QUADRATIC for dy: 0 0 2 x 0 2y
  /// CUBIC for dx: 0 1 0 y 2x 0 2xy y^2 3x^2 0
  /// CUBIC for dy: 0 0 1 x 0 2y x^2 x2y 0 3y^2
  static std::function<std::vector<double>(int, const double *, size_t)>
  polynomialDerivative(PolynomialType baseType);

  /// This function should be used within rbf methods for 2D points. Given a
  /// PolynomialType, this method returns its proper second derivative
  /// polynomial base, following:
  /// CONSTANT: 0
  /// LINEAR: 0 0 0
  /// QUADRATIC for dx: 0 0 0 0 2 0
  /// QUADRATIC for dy: 0 0 0 0 0 2
  /// CUBIC for dx: 0 0 0 0 2 0 2y 0 6x 0
  /// CUBIC for dy: 0 0 0 0 0 2 0 2x 0 6y
  static std::function<std::vector<double>(int, const double *, size_t)>
  polynomialSecondDerivative(PolynomialType baseType);
};

} // namespace furoo

#endif // FUROO_DEFINITIONS_H
