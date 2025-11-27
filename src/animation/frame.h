#ifndef FUROO_ANIMATION_FRAME_H
#define FUROO_ANIMATION_FRAME_H

namespace furoo {

class Frame {
public:
  Frame();
  Frame(unsigned int, double);
  virtual ~Frame();

  unsigned int index() const;
  double timeInterval() const;
  double ellapsedTime() const;
  void advance();
  void advance(unsigned int);
  Frame &operator++();
  Frame operator++(int);

private:
  unsigned int _curr = 0;
  double _timeInterval = 1. / 60.;
};

} // furoo namespace

#endif // FUROO_ANIMATION_FRAME_H
