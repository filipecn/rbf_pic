#ifndef FUROO_ADAPTIVE_PIC2_H
#define FUROO_ADAPTIVE_PIC2_H

#include <apps/pic_simulation2.h>
#include <set>

namespace furoo {

class AdaptivePic2 : public PicSimulation2 {
public:
  AdaptivePic2() = default;
  /// \param region simulation domain
  /// \param minLevel minimum tree level for fluid cells
  /// \param surfaceLevel surface level
  explicit AdaptivePic2(BBox2d region, size_t minLevel, size_t surfaceLevel);
  void resetTree(size_t level);
  std::set<size_t> trackSurface();
  void refineSurface();
  void refineFluid(size_t fluidLevel);

  /// \param doCoarse returns true if cell must be coarsed
  void coarseTree(std::function<bool(Definitions::Material, size_t)> doCoarse);

  /// Must be called before simulation
  void buildTree();

  /// Must be called every simulation step
  void updateTree();

  /// Rearrange particles inside fluid cells for better performance and
  /// interpolation
  virtual void reseedParticles();

  void solvePressure(double dt) override;
  double cfl() const override;

protected:
  size_t _minLevel;
  size_t _surfaceLevel;
};

} // namespace furoo

#endif // FUROO_ADAPTIVE_PIC2_H
