#ifndef FUROO_ARGS_PARSER_H
#define FUROO_ARGS_PARSER_H

#include <map>
#include <vector>
#include <string>

namespace furoo {

/// Helper class to easily program arguments and configuration files
class ArgsParser {
public:
  ArgsParser() = default;

  /// Parses a file content that where lines can be:
  /// variableName<1 space>value
  /// vectorName<1 space>value1<1 space>...<1 space>valuen
  /// \param filename path/to/file
  void parse(std::string filename);

  /// Parses the argument list
  /// \param argc number of words
  /// \param argv word list
  void parse(int argc, char **argv);

  /// Registers a vector of strings variable
  /// \param arg variable name
  void addStringListArgument(std::string arg);

  /// Registers a vector of doubles variable
  /// \param arg variable name
  void addDoubleListArgument(std::string arg);

  /// Registers a vector of integers variable
  /// \param arg variable name
  void addIntListArgument(std::string arg);

  /// Registers a boolean variable
  /// \param arg variable name
  /// \param value default value (in case variable is not parsed)
  void addBoolArgument(std::string arg, bool value);

  /// Registers an integer variable
  /// \param arg variable name
  /// \param value default value (in case variable is not parsed)
  void addIntArgument(std::string arg, int value);

  /// Registers a double variable
  /// \param arg variable name
  /// \param value default value (in case variable is not parsed)
  void addDoubleArgument(std::string arg, double value);

  /// Registers a string variable
  /// \param arg variable name
  /// \param value default value (in case variable is not parsed)
  void addStringArgument(std::string arg, std::string value);

  /// \param arg boolean argument name
  /// \return argument bool value
  bool getBool(std::string arg) const;

  /// \param arg integer argument name
  /// \return argument int value
  int getInt(std::string arg) const;

  /// \param arg double argument name
  /// \return argument double value
  double getDouble(std::string arg) const;

  /// \param arg string argument name
  /// \return argument std::string value
  std::string getString(std::string arg) const;

  /// \param arg string list argument name
  /// \return argument std::vector<std::string> value
  std::vector<std::string> getStringList(std::string arg) const;

  /// \param arg double list argument name
  /// \return argument std::vector<std::string> value
  std::vector<double> getDoubleList(std::string arg) const;

  /// \param arg int list argument name
  /// \return argument std::vector<std::string> value
  std::vector<int> getIntList(std::string arg) const;

private:
  std::map<std::string, std::vector<std::string>> _stringListArguments;
  std::map<std::string, std::vector<double>> _doubleListArguments;
  std::map<std::string, std::vector<int>> _intListArguments;
  std::map<std::string, bool> _boolArguments;
  std::map<std::string, int> _intArguments;
  std::map<std::string, double> _doubleArguments;
  std::map<std::string, std::string> _stringArguments;
};

} // namespace furoo

#endif // FUROO_ARGS_PARSER_H
