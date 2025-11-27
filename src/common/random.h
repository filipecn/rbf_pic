#ifndef FUROO_COMMON_RANDOM_H
#define FUROO_COMMON_RANDOM_H

#include <geometry/bbox.h>
#include <geometry/numeric.h>

namespace furoo {

/** \brief Random Number Generator
 * Implements the "Halton Sequence".
 */
class HaltonSequence {
 public:
  /** \brief Default constructor.
   */
  HaltonSequence() : _base(2), _ind(1) {}
  /** \brief Constructor.
   * \param b base ( > 1)
   */
  HaltonSequence(unsigned int b) : _base(b), _ind(1) {}
  /** \brief
   * \param b base ( > 1)
   */
  void setBase(unsigned int b) {
    _base = b;
    _ind = 1;
  }
  double randomDouble() {
    double result = 0.;
    double f = 1.;
    unsigned int i = _ind++;
    while (i > 0) {
      f /= _base;
      result += f * (i % _base);
      i /= _base;
    }
    return result;
  }

  double randomDouble(double a, double b) { return lerp(randomDouble(), a, b); }

 private:
  unsigned int _base, _ind;
};

class BoxSampler {
 public:
  BoxSampler(size_t bx = 2, size_t by = 3, size_t bz = 5) {
    _hsx = HaltonSequence(bx);
    _hsy = HaltonSequence(by);
    _hsz = HaltonSequence(bz);
  }

  Point2d sample(const BBox2D &bbox) {
    return Point2d(_hsx.randomDouble(bbox.lower().x(), bbox.upper().x()),
                   _hsy.randomDouble(bbox.lower().y(), bbox.upper().y()));
  }

  Point2d sample(const BBox2d &bbox) {
    return Point2d(_hsx.randomDouble(bbox.lower().x(), bbox.upper().x()),
                   _hsy.randomDouble(bbox.lower().y(), bbox.upper().y()));
  }

  Point3d sample(const BBox3d &bbox) {
    return Point3d(_hsx.randomDouble(bbox.lower().x(), bbox.upper().x()),
                   _hsy.randomDouble(bbox.lower().y(), bbox.upper().y()),
                   _hsz.randomDouble(bbox.lower().z(), bbox.upper().z()));
  }

  Point2d sample() { return Point2d(_hsx.randomDouble(), _hsy.randomDouble()); }

 private:
  HaltonSequence _hsx, _hsy, _hsz;
};

} // furoo namespace

#endif // FUROO_COMMON_RANDOM_H
