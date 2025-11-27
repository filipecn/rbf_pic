#ifndef FUROO_ANIMATION_ANIMATION_H
#define FUROO_ANIMATION_ANIMATION_H

#include <animation/frame.h>
#include <common/debug.h>
#include <string>

namespace furoo {

class Animation {
 public:
  Animation();
  /// \param startTime starting time (seconds)
  /// \param duration duration time (seconds)
  /// \param fps frames per second
  Animation(double startTime, double duration, size_t fps = 60);
  /// \param startFrame starting frame
  /// \param frameCount total number of frames
  Animation(size_t startFrame, size_t frameCount);
  virtual ~Animation();
  /// It also updates startFrame and endFrame values
  /// \param fps frames per second
  void setFPS(size_t fps);
  /// \param path directory where animation files will be stored
  void setAnimationFilesPath(const char *path);
  /// \param startTime starting time (seconds)
  /// \param duration duration time (seconds)
  void setAnimationTime(double startTime, double duration);
  /// \param startFrame starting frame
  /// \param frameCount total number of frames
  void setAnimationTime(size_t startFrame, size_t frameCount);
  /// Run entire animation (0s, startTime + duration), saving files from
  /// specified frames.
  void run();
  /// \param dirPath directory path for storing files
  /// \param frame frame number
  virtual void write(std::string dirPath, size_t frame) const {
    // TODO: this method has to become abstract on the furoote!
    UNUSED_VARIABLE(dirPath);
    UNUSED_VARIABLE(frame);
  };
  /// \param dt time duration of frame
  virtual void update(double dt) {
    // TODO: this method has to become abstract on the furoote!
    UNUSED_VARIABLE(dt);
  }

  void update(const Frame &); // TODO: todie!

 protected:
  std::string _animationFilesPath;
  size_t _fps;
  size_t _startFrame;
  size_t _endFrame;

  virtual void onUpdate(const Frame &) = 0; // TODO: todie!
};

} // furoo namespace

#endif // FUROO_ANIMATION_ANIMATION_H
