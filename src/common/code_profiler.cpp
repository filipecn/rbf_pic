#include "code_profiler.h"

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>

namespace furoo {

CodeProfiler CodeProfiler::_codeProfiler;

CodeProfiler &CodeProfiler::instance() { return _codeProfiler; }

void CodeProfiler::profile(std::string name, double value) {
  instance()[name] = value;
}

void CodeProfiler::reset() {
  for (auto &s : instance()._curSnapshot.statistics)
    s.second = 0.;
  for (auto &s : instance()._curSnapshot.series)
    s.second.clear();
}

void CodeProfiler::count(std::string name, size_t n) { instance()[name] += n; }

double &CodeProfiler::operator[](std::string name) {
  if (_curSnapshot.statistics.find(name) == _curSnapshot.statistics.end())
    _curSnapshot.statistics[name] = 0.;
  return _curSnapshot.statistics[name];
}

void CodeProfiler::addToSeries(std::string name, double value) {
  if (_curSnapshot.series.find(name) == _curSnapshot.series.end())
    _curSnapshot.series[name] = std::vector<double>();
  _curSnapshot.series[name].emplace_back(value);
}

void CodeProfiler::createSnapshot(std::string name, std::string filename) {
  instance()._curSnapshot.name = name;
  instance()._history.push_back(instance()._curSnapshot);
  if (filename.empty())
    return;
  std::ofstream fp(filename, std::ofstream::out | std::ofstream::app);
  if (!fp.good())
    return;
  instance().dumpSnapshot(fp, instance()._curSnapshot);
  fp.close();
}

void CodeProfiler::saveHistory(std::string filename, bool csv) {
  if (csv)
    return saveHistoryToCsv(filename);

  std::ofstream fp(filename, std::ofstream::out);
  if (!fp.good())
    return;
  for (auto s : instance()._history)
    instance().dumpSnapshot(fp, s);
  fp.close();
}

void CodeProfiler::saveHistoryToCsv(std::string filename)
{
  std::ofstream fp(filename, std::ofstream::out);
  if (!fp.good())
    return;
  if (instance()._history.size() > 2)
    instance().writeCsvHeader(fp, instance()._history[1]);
  for (auto s : instance()._history)
    instance().writeCsvRow(fp, s);
  fp.close();
}

void CodeProfiler::writeCsvHeader(std::ostream &fp, Snapshot &snapshot)
{
  fp << "frame";
  for (auto &s : snapshot.statistics)
    fp << "," << s.first.c_str();
  fp << "\n";
}

void CodeProfiler::writeCsvRow(std::ostream &fp, Snapshot &snapshot)
{
  fp << snapshot.name;
  for (auto &s : snapshot.statistics)
    fp << "," << s.second;
  fp << "\n";
}

void CodeProfiler::dumpSnapshot(std::ostream &fp, Snapshot &snapshot) {
  fp << snapshot.name << "\t" << snapshot.statistics.size() << "\n";
  for (auto &s : snapshot.statistics)
    fp << s.first.c_str() << "\t" << s.second << "\n";
  fp << snapshot.series.size() << " ";
  for (auto &s : snapshot.series) {
    fp << s.first << "\t" << s.second.size() << "\n";
    for (auto v : s.second)
      fp << v << "\n";
  }
  fp << std::endl;
}

CodeProfiler::CodeProfiler() = default;

} // namespace furoo
