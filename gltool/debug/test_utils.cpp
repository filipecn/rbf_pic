#include "test_utils.h"
#include "external/lodepng.h"

#include <aergia/aergia.h>
using namespace furoo;

void Injector::setSpacing(double s) {
  _spacing = s;
}

void Injector::setupBoxShape(Point2d pMin, Point2d pMax) {
  for (double y = pMin.y(); y <= pMax.y(); y += _spacing)
    for (double x = pMin.x(); x <= pMax.x(); x += _spacing)
      _particleSystem->addParticle(Point2d(x, y));
}

void Injector::setupCircleShape(Point2d center, double radius) {
  Point2d pMin = center - Vector2d(radius, radius);
  Point2d pMax = center + Vector2d(radius, radius);
  for (double y = pMin.y(); y <= pMax.y(); y += _spacing)
    for (double x = pMin.x(); x <= pMax.x(); x += _spacing)
      if (distance(center, Point2d(x, y)) <= radius)
        _particleSystem->addParticle(Point2d(x, y));
}

void Injector::loadFromFile(const char *filename) {
//  IO::loadFromPV(filename, _particleSystem);
  UNUSED_VARIABLE(filename);
}

void saveFramebuffer(const char *filename,
                     unsigned int WIDTH,
                     unsigned int HEIGHT) {
  auto *imageData = (unsigned char *) malloc(WIDTH * HEIGHT * 4);
  glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, imageData);
  for (uint x = 0; x < WIDTH; x++)
    for (uint y = 0; y < HEIGHT / 2; y++)
      for (uint k = 0; k < 4; k++) {
        unsigned char tmp = imageData[4 * WIDTH * (HEIGHT - 1 - y) + 4 * x + k];
        imageData[4 * WIDTH * (HEIGHT - 1 - y) + 4 * x + k] =
            imageData[4 * WIDTH * y + 4 * x + k];
        imageData[4 * WIDTH * y + 4 * x + k] = tmp;
      }
  //Encode the image
  std::cout << filename << std::endl;
  unsigned error = lodepng::encode(filename, imageData, WIDTH, HEIGHT);
  //if there's an error, display it
  if (error)
    std::cout << "encoder error " << error << ": " << lodepng_error_text(error)
              << std::endl;
  free(imageData);
}

Point2d advectExact(Point2d p,
                    double dt,
                    std::function<Vector2d(Point2d)> velocity,
                    unsigned int subStepsCount) {
  unsigned int numSubSteps = subStepsCount;
  double subdt = dt / numSubSteps;
  Point2d pt = p;
  for (unsigned int t = 0; t < numSubSteps; t++) {
    Vector2d v0 = velocity(pt);
    // Mid-point rule
    Point2d midPt = pt + 0.5 * subdt * v0;
    Vector2d midVel = velocity(midPt);
    pt += subdt * midVel;
  }
  return pt;
}