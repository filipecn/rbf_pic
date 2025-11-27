#include "animation/physics_animation.h"
#include <iostream>

namespace furoo {

PhysicsAnimation::PhysicsAnimation() {}

PhysicsAnimation::~PhysicsAnimation() {}

void PhysicsAnimation::advanceFrame() {
  Frame f = _curr;
  update(++f);
}

Frame PhysicsAnimation::currentFrame() const { return _curr; }

void PhysicsAnimation::setCurrentFrame(const Frame &f) { _curr = f; }

double PhysicsAnimation::currentTime() const { return _currentTime; }

void PhysicsAnimation::onInitialize() {}

void PhysicsAnimation::onUpdate(const Frame &frame) {
  if (frame.index() > _curr.index()) {
    if (_curr.index() == 0) {
      initialize();
    }

    unsigned int numberOfFrames = frame.index() - _curr.index();

    for (unsigned int i = 0; i < numberOfFrames; ++i)
      advanceTimeStep(frame.timeInterval());

    _curr = frame;
  }
}

void PhysicsAnimation::advanceTimeStep(double dt) {
  _currentTime = _curr.ellapsedTime();
  // Perform adaptive time-stepping
  double remainingTime = dt;
  // TODO check this 1e-8 later
  while (remainingTime > 1e-8) {
    unsigned int numSteps = numberOfSubTimeSteps(remainingTime);
    double actualTimeInterval = remainingTime / static_cast<double>(numSteps);
    onAdvanceTimeStep(actualTimeInterval);
    remainingTime -= actualTimeInterval;
    _currentTime += actualTimeInterval;
  }
}

void PhysicsAnimation::initialize() { onInitialize(); }

} // furoo namespace
