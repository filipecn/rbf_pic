#ifndef FUROO_ADAPTIVE_CENTERED_PIC2_H
#define FUROO_ADAPTIVE_CENTERED_PIC2_H

#include <apps/adaptive_pic2.h>

namespace furoo {

class AdaptiveCenteredPic2 : public AdaptivePic2 {
 public:
  /// \param region simulation domain
  /// \param minLevel minimum tree level for fluid cells
  /// \param surfaceLevel surface level
  explicit AdaptiveCenteredPic2(BBox2d region,
                                size_t minLevel,
                                size_t surfaceLevel);
  void computeExternalForces(double dt) override;
  void defineBoundaryConditions(Definitions::Side side, Vector2d vel) override;
  void applyBoundaryCondition() override;
  void transferFromParticlesToGrid() override;
  void transferFromGridToParticles() override;
  void solvePressure(double dt) override;
  double cfl() const override;
  void loadFromFile(const char *filename);
 private:
  size_t _cellGradientFieldIds[2];
  size_t _cellVelocityFieldIds[2];
};

}

#endif //FUROO_ADAPTIVE_CENTERED_PIC2_H
