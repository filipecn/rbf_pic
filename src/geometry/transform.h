#ifndef FUROO_GEOMETRY_TRANSFORM_H
#define FUROO_GEOMETRY_TRANSFORM_H

#include <common/debug.h>
#include <geometry/bbox.h>
#include <geometry/matrix.h>
#include <geometry/point.h>
#include <geometry/vector.h>

namespace furoo {

template<int D>
class Transform {
 public:
  Transform();
  /// \param m matrix of transformation
  /// \param m_inv inverse of matrix of transformation
  Transform(Matrix<double, D + 1> m, Matrix<double, D + 1> m_inv);
  /// \return transformation matrix
  Matrix<double, D + 1> matrix() const;
  /// \return inverse transformation matrix
  Matrix<double, D + 1> inverseMatrix() const;
  /// \return inverse transformation
  Transform<D> inverse() const;
  /// Apply transform to point
  /// \param p point
  /// \return T(p)
  Point<double, D> operator()(const Point<double, D> &p) const;
  /// Apply transform to vector
  /// \param v vector
  /// \return T(v)
  Vector<double, D> operator()(const Vector<double, D> &v) const;
  /// \param b bounding box
  /// \return T(b)
  BBox<double, D> operator()(const BBox<double, D> &b) const;
  /// \param v vector to be reflected
  /// \param n normal
  /// \return reflection of v over normal n
  static Vector<double, D> applyReflection(const Vector<double, D> &v,
                                           const Vector<double, D> &n);
  /// \param s scale vector
  /// \return scale transform based on s
  static Transform<D> scale(Vector<double, D> s);

  /// Translation transform
  /// \param t translation vector
  /// \return translation transform based on t
  static Transform<D> translate(Vector<double, D> t);

  /// Rotation over cartesian axis transform
  /// \param axis (0 = x, 1 = y and 2 = z)
  /// \param angle in degrees
  /// \return rotation transform
  static Transform<D> rotate(size_t axis,  double angle);

  /// Computes the transformation from a origin box to a unit box [0,1]^D
  /// \param b origin box
  /// \return transformation
  static Transform<D> toUnitBox(BBox<double, D> b);

 private:
  Matrix<double, D + 1> _m, _m_inv;
};

#include "geometry/transform.inl"

typedef Transform<2> Transform2;
typedef Transform<3> Transform3;

class Transform2D {
 public:
  Transform2D();
  Transform2D(const Transform2D &other);
  Transform2D(const Matrix3d &mat, const Matrix3d inv_mat);
  /** Builds a transform that takes bbox to the unit bbox
   * \param bbox
   */
  Transform2D(const BBox2D &bbox);
  void reset();
  Transform2D applyTranslate(const Vector2d &d);
  Transform2D applyScale(double x, double y);
  Transform2D applyRotate(double angle);
  Transform2D inverse() const;
  void operator()(const Point2d &p, Point2d *r) const;
  void operator()(const Vector2d &v, Vector2d *r) const;
  Vector2d operator()(const Vector2d &v) const;
  Point2d operator()(const Point2d &p) const;
  Point2d operator()(const Point2i &p) const;
  BBox2D operator()(const BBox2D &b) const;
  Transform2D &operator=(const Transform2D &t);
  Transform2D &operator*=(const Transform2D &t);
  Transform2D operator*(const Transform2D &t) const;
  Matrix3d getMatrix() const;

  static Transform2D scale(double x, double y);
  static Transform2D rotate(double angle);
  static Transform2D translate(const Vector2d &v);
  /** \brief  reflects v off the curve with normal n
   * \param v
   * \param n normalized normal
   * \returns reflection
   */
  static Vector2d applyReflection(const Vector2d &v, const Vector2d &n);

 private:
  Matrix3d _m, _m_inv;
};

class Transform3D {
 public:
  Transform3D();
  Transform3D(const Transform3D &other);
  Transform3D(const Matrix4d &mat, const Matrix4d inv_mat);
  /** Builds a transform that takes bbox to the unit bbox
   * \param bbox
   */
  Transform3D(const BBox3D &bbox);
  void reset();
  Transform3D applyTranslate(const Vector3d &d);
  Transform3D applyScale(double x, double y, double z);
  Transform3D applyRotate(double angle, Vector3d axis);
  Transform3D inverse() const;
  void operator()(const Point3d &p, Point3d *r) const;
  void operator()(const Vector3d &v, Vector3d *r) const;
  Vector3d operator()(const Vector3d &v) const;
  Point3d operator()(const Point3d &p) const;
  Point3d operator()(const Point3i &p) const;
  BBox3D operator()(const BBox3D &b) const;
  Transform3D &operator=(const Transform3D &t);
  Transform3D &operator*=(const Transform3D &t);
  Transform3D operator*(const Transform3D &t) const;
  Matrix4d getMatrix() const;

  static Transform3D scale(double x, double y, double z);
  static Transform3D rotate(double angle, Vector3d axis);
  static Transform3D translate(const Vector3d &v);
  /** \brief  reflects v off the curve with normal n
   * \param v
   * \param n normalized normal
   * \returns reflection
   */
  static Vector3d applyReflection(const Vector3d &v, const Vector3d &n);

 private:
  Matrix4d _m, _m_inv;
};

} // furoo namespace

#endif
