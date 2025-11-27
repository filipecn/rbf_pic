#include "definitions.h"
#include <common/debug.h>
#include <common/definitions.h>

namespace furoo {

Definitions::Material Definitions::materialFromId(size_t id) {
  switch (id) {
  case 0:
    return Material::FLUID;
  case 1:
    return Material::AIR;
  case 2:
    return Material::SOLID;
  default:
    return Material::CUSTOM;
  }
  return Material::CUSTOM;
}

Definitions::Side Definitions::sideFromId(size_t id) {
  switch (id) {
  case 0:
    return Definitions::Side::LEFT;
  case 1:
    return Definitions::Side::RIGHT;
  case 2:
    return Definitions::Side::TOP;
  case 3:
    return Definitions::Side::BOTTOM;
  case 4:
    return Definitions::Side::FRONT;
  case 5:
    return Definitions::Side::BACK;
  default:
    return Definitions::Side::CUSTOM;
  }
  return Definitions::Side::CUSTOM;
}

Definitions::Orientation
Definitions::orthogonalOrientation(Definitions::Orientation axisA,
                                   Definitions::Orientation axisB) {
  if (axisA == axisB)
    return Definitions::Orientation::CUSTOM;
  if ((axisA == Definitions::Orientation::HORIZONTAL &&
       axisB == Definitions::Orientation::VERTICAL) ||
      (axisB == Definitions::Orientation::HORIZONTAL &&
       axisA == Definitions::Orientation::VERTICAL))
    return Definitions::Orientation::DEPTH;
  if ((axisA == Definitions::Orientation::HORIZONTAL &&
       axisB == Definitions::Orientation::DEPTH) ||
      (axisB == Definitions::Orientation::HORIZONTAL &&
       axisA == Definitions::Orientation::DEPTH))
    return Definitions::Orientation::VERTICAL;
  return Definitions::Orientation::HORIZONTAL;
}

std::vector<Definitions::Side>
Definitions::orientationSides(Definitions::Orientation orientation) {
  switch (orientation) {
  case Definitions::Orientation::VERTICAL:
    return {Side::BOTTOM, Side::TOP};
  case Definitions::Orientation::HORIZONTAL:
    return {Side::LEFT, Side::RIGHT};
  case Definitions::Orientation::DEPTH:
    return {Side::BACK, Side::FRONT};
  case Definitions::Orientation::CUSTOM:
    return {};
  }
  return {};
}

Definitions::Side Definitions::oppositeSide(Definitions::Side s) {
  switch (s) {
  case Side::LEFT:
    return Side::RIGHT;
  case Side::RIGHT:
    return Side::LEFT;
  case Side::TOP:
    return Side::BOTTOM;
  case Side::BOTTOM:
    return Side::TOP;
  case Side::FRONT:
    return Side::BACK;
  case Side::BACK:
    return Side::FRONT;
  default:
    return Side::CUSTOM;
  }
  return Side::CUSTOM;
}

Definitions::Orientation Definitions::sideOrientation(Definitions::Side s) {
  switch (s) {
  case Definitions::Side::LEFT:
    return Definitions::Orientation::HORIZONTAL;
  case Definitions::Side::RIGHT:
    return Definitions::Orientation::HORIZONTAL;
  case Definitions::Side::TOP:
    return Definitions::Orientation::VERTICAL;
  case Definitions::Side::BOTTOM:
    return Definitions::Orientation::VERTICAL;
  case Definitions::Side::FRONT:
    return Definitions::Orientation::DEPTH;
  case Definitions::Side::BACK:
    return Definitions::Orientation::DEPTH;
  default:
    return Definitions::Orientation::CUSTOM;
  }
  return Definitions::Orientation::CUSTOM;
}

std::function<std::vector<double>(const double *, size_t)>
Definitions::polynomial(Definitions::PolynomialType baseType) {
  switch (baseType) {
  case Definitions::PolynomialType::ZERO:
    return [](const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(D);
      UNUSED_VARIABLE(point);
      return {};
    };
  case Definitions::PolynomialType::CONSTANT:
    return [](const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(D);
      UNUSED_VARIABLE(point);
      return {1.};
    };
    break;
  case Definitions::PolynomialType::LINEAR:
    return [](const double *point, size_t D) -> std::vector<double> {
      std::vector<double> base(1 + D, 1);
      for (size_t i = 0; i < D; i++)
        base[i + 1] = point[i];
      return base;
    };
    break;
  case Definitions::PolynomialType::QUADRATIC:
    return [](const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(D);
      double x = point[0];
      double y = point[1];
      if (D == 2) {
        std::vector<double> base = {1, x, y, x * y, x * x, y * y};
        return base;
      }
      double z = point[2];
      std::vector<double> base = {1,     x,     y,     z,     x * y,
                                  x * z, y * z, x * x, y * y, z * z};
      return base;
    };
    break;
  default:
    THROW(false, "Definitions: polynomial type not valid!");
  }
  return nullptr;
}

std::function<std::vector<double>(int, const double *, size_t)>
Definitions::polynomialDerivative(Definitions::PolynomialType baseType) {
  switch (baseType) {
  case Definitions::PolynomialType::CONSTANT:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(point);
      UNUSED_VARIABLE(d);
      UNUSED_VARIABLE(D);
      return {0.};
    };
    break;
  case Definitions::PolynomialType::LINEAR:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(point);
      std::vector<double> base(1 + D, 0);
      base[d + 1] = 1.;
      return base;
    };
    break;
  case Definitions::PolynomialType::QUADRATIC:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(D);
      double x = point[0];
      double y = point[1];
      if (D == 2) {
        if (d == 0)
          return {0, 1, 0, y, 2 * x, 0};
        THROW(d == 1, "Definitions::polynomialDerivative QUDRATIC type invalid "
                      "derivative dimension");
        return {0, 0, 1, x, 0, 2 * y};
      } else {
        double z = point[2];
        if (d == 0)
          return {0, 1, 0, 0, y, z, 0, 2 * x, 0, 0};
        if (d == 1)
          //     {1, x, y, z, x * y, x * z, y * z, x * x, y * y, z * z};
          return {0, 0, 1, 0, x, 0, z, 0, 2 * y, 0};
        return {0, 0, 0, 1, 0, x, y, 0, 0, 2 * z};
      }
    };
    break;
  default:
    THROW(false, "Definitions: polynomial type not valid!");
  }
  return nullptr;
}

std::function<std::vector<double>(int, const double *, size_t)>
Definitions::polynomialSecondDerivative(PolynomialType baseType) {
  switch (baseType) {
  case Definitions::PolynomialType::CONSTANT:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(point);
      UNUSED_VARIABLE(d);
      UNUSED_VARIABLE(D);
      return {0.};
    };
    break;
  case Definitions::PolynomialType::LINEAR:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(point);
      UNUSED_VARIABLE(d);
      if (D == 2)
        return {0., 0., 0.};
      return {0, 0, 0, 0};
    };
    break;
  case Definitions::PolynomialType::QUADRATIC:
    return [](int d, const double *point, size_t D) -> std::vector<double> {
      UNUSED_VARIABLE(point);
      if (D == 2) {
        if (d == 0)
          return {0, 0, 0, 0, 2, 0};
        THROW(d == 1, "Definitions::polynomialSecondDerivative invalid "
                      "derivative dimension");
        return {0, 0, 0, 0, 0, 2};
      } else {
        if (d == 0)
          return {0, 0, 0, 0, 0, 0, 0, 2, 0, 0};
        if (d == 1)
          //{1, x, y, z, x * y, x * z, y * z, x * x, y * y, z * z};
          return {0, 0, 0, 0, 0, 0, 0, 0, 2, 0};
        return {0, 0, 0, 0, 0, 0, 0, 0, 0, 2};
      }
    };
    break;
  default:
    THROW(false, "Definitions: polynomial type not valid!");
  }
  return nullptr;
}

Definitions::MeshLocation Definitions::meshLocationFromId(size_t id) {
  switch (id) {
  case 0:
    return MeshLocation::CELL_CENTER;
  case 1:
    return MeshLocation::VERTEX_CENTER;
  case 2:
    return MeshLocation::FACE_CENTER;
  case 3:
    return MeshLocation::HORIZONTAL_FACE_CENTER;
  case 4:
    return MeshLocation::VERTICAL_FACE_CENTER;
  case 5:
    return MeshLocation::DEPTH_FACE_CENTER;
  case 6:
    return MeshLocation::CUSTOM;
  default:
    return MeshLocation::CUSTOM;
  }
  return MeshLocation::CUSTOM;
}

Definitions::Boundary Definitions::boundaryFromId(size_t id) {
  switch (id) {
  case 0:
    return Definitions::Boundary::DIRICHLET;
  case 1:
    return Definitions::Boundary::NEUMANN;
  case 2:
    return Definitions::Boundary::NONE;
  case 3:
    return Definitions::Boundary::CUSTOM;
  default:
    return Definitions::Boundary::NONE;
  }
  return Definitions::Boundary::NONE;
}

} // namespace furoo
