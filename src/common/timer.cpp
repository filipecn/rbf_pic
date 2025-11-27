#include "common/timer.h"
#include <iomanip>
#include <sstream>
#include <common/code_profiler.h>

namespace furoo {

Timer::Timer() { _startingPoint = _clock.now(); }

double Timer::ellapsedTimeInSeconds() const {
  auto end = std::chrono::high_resolution_clock::now();
  auto count = std::chrono::duration_cast<std::chrono::milliseconds>(
      end - _startingPoint)
      .count();
  return count / 1000.0;
}

void Timer::reset() { _startingPoint = _clock.now(); }

TimeLogger::TimeLogger() {
  _timer.reset();
}

void TimeLogger::reset() {
  for (auto &l : _log) l.second = 0.;
  _timer.reset();
}

void TimeLogger::log(const char *name) {
  std::string n(name);
  double ellapsedTime = _timer.ellapsedTimeInSeconds();
  auto it = _log.find(n);
  if (it == _log.end())
    _log[name] = ellapsedTime;
  else it->second += ellapsedTime;
  _timer.reset();
}

void TimeLogger::report() {
  std::priority_queue<std::pair<double, std::string>> q;
  double totalTime = 0.;
  for (auto l : _log) {
    std::pair<double, std::string> p(l.second, l.first);
    q.push(p);
    totalTime += l.second;
  }
  std::cerr << "Total time: " << totalTime << " seconds\n";
  std::ostringstream percentLine, timeLine;
  while (!q.empty()) {
    percentLine << q.top().second << " - " << std::setprecision(3)
                << 100 * (q.top().first / totalTime) << std::fixed
                << "%\t";
    timeLine << q.top().second << " - " << std::setprecision(3)
             << q.top().first << std::fixed
             << "s\t";
    q.pop();
  }
  std::cerr << timeLine.str() << std::endl << percentLine.str() << std::endl;
}

void TimeLogger::createSnapshot() {
  for (auto l : _log)
    CodeProfiler::instance()[l.first.c_str()] = l.second;
}

double TimeLogger::operator[](const char *name) const {
  auto it = _log.find(std::string(name));
  if (it == _log.end())
    return 0;
  return it->second;
}

} // furoo namespace
