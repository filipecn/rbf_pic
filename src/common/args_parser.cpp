#include "args_parser.h"
#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <utility>

namespace furoo {

void ArgsParser::parse(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if (_boolArguments.find(argv[i]) != _boolArguments.end()) {
      bool v = _boolArguments[argv[i]];
      std::istringstream(argv[i + 1]) >> v;
      _boolArguments[argv[i]] = v;
      i++;
    } else if (_doubleArguments.find(argv[i]) != _doubleArguments.end()) {
      double v = _doubleArguments[argv[i]];
      std::istringstream(argv[i + 1]) >> v;
      _doubleArguments[argv[i]] = v;
      i++;
    } else if (_intArguments.find(argv[i]) != _intArguments.end()) {
      int v = _intArguments[argv[i]];
      std::istringstream(argv[i + 1]) >> v;
      _intArguments[argv[i]] = v;
      i++;
    } else if (_stringArguments.find(argv[i]) != _stringArguments.end()) {
      std::string v = _stringArguments[argv[i]];
      std::istringstream(argv[i + 1]) >> v;
      _stringArguments[argv[i]] = v;
      i++;
    }
  }
}

void ArgsParser::parse(std::string filename) {
  std::ifstream fp(filename, std::ifstream::in);
  std::cerr << "================= PARSING " << filename
            << " =================\n";
  auto split = [](const std::string &s,
                  char delim) -> std::vector<std::string> {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (getline(ss, item, delim))
      if (item.size())
        elems.push_back(item);
    return elems;
  };
  while (fp.good()) {
    std::string name;
    fp >> name;
    if (name[0] == '#') {
      std::string line;
      std::getline(fp, line);
      std::cerr << "COMMENT " << line << std::endl;

    } else if (_stringListArguments.find(name) != _stringListArguments.end()) {
      std::string line;
      std::getline(fp, line);
      _stringListArguments[name] = split(line, ' ');
      std::cerr << "vector<string> " << name << " [";
      for (auto s : _stringListArguments[name])
        std::cerr << s << ',';
      std::cerr << "]\n";
    } else if (_doubleListArguments.find(name) != _doubleListArguments.end()) {
      std::string line;
      std::getline(fp, line);
      _doubleListArguments[name] = std::vector<double>();
      auto elems = split(line, ' ');
      for (auto e : elems) {
        std::stringstream s;
        s << e;
        double value = 0.;
        s >> value;
        _doubleListArguments[name].emplace_back(value);
      }
      std::cerr << "vector<double> " << name << " [";
      for (auto s : _doubleListArguments[name])
        std::cerr << s << ',';
      std::cerr << "]\n";
    } else if (_intListArguments.find(name) != _intListArguments.end()) {
      std::string line;
      std::getline(fp, line);
      _intListArguments[name] = std::vector<int>();
      auto elems = split(line, ' ');
      for (auto e : elems) {
        std::stringstream s;
        s << e;
        int value = 0;
        s >> value;
        _intListArguments[name].emplace_back(value);
      }
      std::cerr << "vector<int> " << name << " [";
      for (auto s : _intListArguments[name])
        std::cerr << s << ',';
      std::cerr << "]\n";
    } else if (_boolArguments.find(name) != _boolArguments.end()) {
      bool v = _boolArguments[name];
      fp >> v;
      _boolArguments[name] = v;
      std::cerr << "bool " << name << " " << v << std::endl;
    } else if (_doubleArguments.find(name) != _doubleArguments.end()) {
      double v = _doubleArguments[name];
      fp >> v;
      _doubleArguments[name] = v;
      std::cerr << "double " << name << " " << v << std::endl;
    } else if (_intArguments.find(name) != _intArguments.end()) {
      int v = _intArguments[name];
      fp >> v;
      _intArguments[name] = v;
      std::cerr << "int " << name << " " << v << std::endl;
    } else if (_stringArguments.find(name) != _stringArguments.end()) {
      std::string v = _stringArguments[name];
      fp >> v;
      _stringArguments[name] = v;
      std::cerr << "string " << name << " " << v << std::endl;
    } else {
      std::string line;
      std::getline(fp, line);
      std::cerr << "UNDEFINED " << line << std::endl;
    }
  }
  std::cerr << "================= END PARSING =================\n";
  fp.close();
}

void ArgsParser::addStringListArgument(std::string arg) {
  _stringListArguments[arg] = std::vector<std::string>();
}

void ArgsParser::addDoubleListArgument(std::string arg) {
  _doubleListArguments[arg] = std::vector<double>();
}

void ArgsParser::addIntListArgument(std::string arg) {
  _intListArguments[arg] = std::vector<int>();
}

void ArgsParser::addBoolArgument(std::string arg, bool value) {
  _boolArguments[arg] = value;
}

void ArgsParser::addIntArgument(std::string arg, int value) {
  _intArguments[arg] = value;
}

void ArgsParser::addDoubleArgument(std::string arg, double value) {
  _doubleArguments[arg] = value;
}

void ArgsParser::addStringArgument(std::string arg, std::string value) {
  _stringArguments[arg] = std::move(value);
}

bool ArgsParser::getBool(std::string arg) const {
  auto it = _boolArguments.find(arg);
  if (it != _boolArguments.end())
    return it->second;
  return false;
}

int ArgsParser::getInt(std::string arg) const {
  auto it = _intArguments.find(arg);
  if (it != _intArguments.end())
    return it->second;
  return 0;
}

double ArgsParser::getDouble(std::string arg) const {
  auto it = _doubleArguments.find(arg);
  if (it != _doubleArguments.end())
    return it->second;
  return 0.;
}

std::string ArgsParser::getString(std::string arg) const {
  auto it = _stringArguments.find(arg);
  if (it != _stringArguments.end())
    return it->second;
  return std::string();
}

std::vector<std::string> ArgsParser::getStringList(std::string arg) const {
  auto it = _stringListArguments.find(arg);
  if (it != _stringListArguments.end())
    return it->second;
  return std::vector<std::string>();
}

std::vector<double> ArgsParser::getDoubleList(std::string arg) const {
  auto it = _doubleListArguments.find(arg);
  if (it != _doubleListArguments.end())
    return it->second;
  return std::vector<double>();
}

std::vector<int> ArgsParser::getIntList(std::string arg) const {
  auto it = _intListArguments.find(arg);
  if (it != _intListArguments.end())
    return it->second;
  return std::vector<int>();
}

} // namespace furoo