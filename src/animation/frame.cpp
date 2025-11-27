#include "animation/frame.h"
#include <common/debug.h>

namespace furoo {

Frame::Frame() {}

Frame::Frame(unsigned int newIndex, double newTimeInterval)
    : _curr(newIndex), _timeInterval(newTimeInterval) {}

Frame::~Frame() {}

unsigned int Frame::index() const { return _curr; }

double Frame::ellapsedTime() const { return _curr * _timeInterval; }

double Frame::timeInterval() const { return _timeInterval; }

void Frame::advance() { ++_curr; }

void Frame::advance(unsigned int delta) { _curr += delta; }

Frame &Frame::operator++() {
  advance();
  return *this;
}

Frame Frame::operator++(int i) {
  UNUSED_VARIABLE(i);
  Frame result = *this;
  advance();
  return result;
}

} // furoo namespace
