#include "animation/animation.h"

namespace furoo {

Animation::Animation()
    : _animationFilesPath(std::string("animation/")),
      _fps(60u),
      _startFrame(0u),
      _endFrame(120u) {}

Animation::~Animation() = default;

Animation::Animation(double startTime, double duration, size_t fps)
    : Animation() {
  _fps = fps;
  _startFrame = startTime * _fps;
  _endFrame = _startFrame + duration * _fps;
}

Animation::Animation(size_t startFrame, size_t frameCount) : Animation() {
  _startFrame = startFrame;
  _endFrame = _startFrame + frameCount;
}

void Animation::setFPS(size_t fps) {
  _startFrame = _startFrame * _fps / fps;
  _endFrame = _endFrame * _fps / fps;
  _fps = fps;
}

void Animation::setAnimationFilesPath(const char *path) {
  _animationFilesPath = path;
}

void Animation::setAnimationTime(double startTime, double duration) {
  _startFrame = startTime * _fps;
  _endFrame = _startFrame + duration * _fps;
}

void Animation::setAnimationTime(size_t startFrame, size_t frameCount) {
  _startFrame = startFrame;
  _endFrame = _startFrame + frameCount;
}

void Animation::run() {
  double dt = 1. / _fps;
  for (size_t i = 0; i < _endFrame; i++) {
    update(dt);
    if (i >= _startFrame && i <= _endFrame)
      write(_animationFilesPath, i);
  }
}

void Animation::update(const Frame &frame) { onUpdate(frame); }
} // furoo namespace
