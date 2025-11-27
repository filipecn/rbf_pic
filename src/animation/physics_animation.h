#ifndef FUROO_ANIMATION_PHYSICS_ANIMATION_H
#define FUROO_ANIMATION_PHYSICS_ANIMATION_H

#include <animation/animation.h>

namespace furoo {

class PhysicsAnimation : public Animation {
public:
  PhysicsAnimation();
  virtual ~PhysicsAnimation();
  void advanceFrame();
  Frame currentFrame() const;
  void setCurrentFrame(const Frame &);
  double currentTime() const;

protected:
  virtual void onAdvanceTimeStep(double) = 0;
  virtual unsigned int numberOfSubTimeSteps(double) const = 0;
  virtual void onInitialize();

private:
  Frame _curr;
  double _currentTime = 0.;

  void onUpdate(const Frame &);
  void advanceTimeStep(double);
  void initialize();
};

} // furoo namespace

#endif // FUROO_ANIMATION_PHYSICS_ANIMATION_H
