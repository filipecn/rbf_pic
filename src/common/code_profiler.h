#ifndef FUROO_CODE_PROFILER_H
#define FUROO_CODE_PROFILER_H

#include <map>
#include <string>
#include <vector>

namespace furoo {

/// Singleton class for easy code profiling. Statistics such as counting, timing
/// and other measures can be added to a log system and be dumped into files.
/// Useful for debugging purposes.
class CodeProfiler {
public:
  /// \return class instance
  static CodeProfiler &instance();

  static void profile(std::string name, double value);

  static void count(std::string name, size_t n = 1);

  /// Deleted public constructor
  CodeProfiler(CodeProfiler &) = delete;

  /// Deleted public assign operator
  void operator=(const CodeProfiler &) = delete;

  /// Sets field values to zero
  static void reset();

  /// Read or write a statistics
  /// \param name statistics name
  /// \return reference to statistics value
  double &operator[](std::string name);

  /// Adds a new value to a specific series of statistics
  /// \param name series name
  /// \param value value to be added
  void addToSeries(std::string name, double value);

  /// Saves the current snapshot to history
  /// \param name if passed, appends the snapshot into file
  static void createSnapshot(std::string name, std::string filename = "");

  /// Saves history into file
  /// \param filename path to file
  static void saveHistory(std::string filename, bool csv = false);

  static void saveHistoryToCsv(std::string filename);

private:
  CodeProfiler();
  struct Snapshot {
    std::string name;
    std::map<std::string, double> statistics;
    std::map<std::string, std::vector<double>> series;
  };
  void dumpSnapshot(std::ostream &fp, Snapshot &snapshot);
  void writeCsvHeader(std::ostream &fp, Snapshot &snapshot);
  void writeCsvRow(std::ostream &fp, Snapshot &snapshot);

  static CodeProfiler _codeProfiler;
  Snapshot _curSnapshot;
  std::vector<Snapshot> _history;
};

} // namespace furoo

#endif // FUROO_CODE_PROFILER_H
