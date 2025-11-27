#ifndef FUROO_COMMON_TIMER_H
#define FUROO_COMMON_TIMER_H

#include <chrono>
#include <map>
#include <queue>
#include <iostream>

namespace furoo {

class Timer {
 public:
  Timer();
  double ellapsedTimeInSeconds() const;
  void reset();

 private:
  std::chrono::high_resolution_clock _clock;
  std::chrono::high_resolution_clock::time_point _startingPoint;
};

class TimeLogger {
 public:
  TimeLogger();
  void reset();
  void log(const char *name);
  void report();
  void createSnapshot();
  double operator[](const char *name) const;
 private:
  Timer _timer;
  std::map<std::string, double> _log;
};

} // furoo namespace

#endif // FUROO_COMMON_TIMER_H
