#ifndef AOCOMMON_UNITS_FLUX_DENSITY_H_
#define AOCOMMON_UNITS_FLUX_DENSITY_H_

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace aocommon {
namespace units {
class FluxDensity {
 public:
  enum Unit { kKiloJansky, kJansky, kMilliJansky, kMicroJansky, kNanoJansky };
  /**
   * Parse the string to a flux density, possibly with unit specification, and
   * return in Jansky.
   * @return The flux density (Jy)
   */
  static double Parse(const std::string& s,
                      const std::string& value_description,
                      Unit default_unit = Unit::kJansky);

  static std::string ToNiceString(double value_jansky);

 private:
  static size_t FindNumberEnd(const std::string& s);
  static bool IsDigit(const char c) { return c >= '0' && c <= '9'; }
  static bool IsWhitespace(const char c) { return c == ' ' || c == '\t'; }
};

inline std::string FluxDensity::ToNiceString(double value_jansky) {
  std::ostringstream str;
  if (value_jansky == 0.0)
    return "0 Jy";
  else {
    if (value_jansky < 0.0) {
      str << "-";
      value_jansky = -value_jansky;
    }
    if (value_jansky >= 1000.0)
      str << std::round(value_jansky * 0.1) / 100.0 << " KJy";
    else if (value_jansky >= 1.0)
      str << std::round(value_jansky * 1e2) / 100.0 << " Jy";
    else if (value_jansky >= 1e-3)
      str << std::round(value_jansky * 1e5) / 100.0 << " mJy";
    else if (value_jansky >= 1e-6)
      str << std::round(value_jansky * 1e8) / 100.0 << " µJy";
    else if (value_jansky >= 1e-9)
      str << std::round(value_jansky * 1e11) / 100.0 << " nJy";
    else
      str << value_jansky << " Jy";
    return str.str();
  }
}

inline double FluxDensity::Parse(const std::string& s,
                                 const std::string& value_description,
                                 Unit default_unit) {
  size_t end = FindNumberEnd(s);
  if (end == 0) {
    throw std::runtime_error("Error parsing " + value_description);
  }
  std::string number = s.substr(0, end);
  const double val = std::atof(number.c_str());
  // Skip whitespace after number
  const char* c = s.c_str();
  while (IsWhitespace(c[end])) {
    ++end;
  }
  const std::string unit_string = std::string(&c[end]);
  std::string unit_string_lower = unit_string;
  std::for_each(unit_string_lower.begin(), unit_string_lower.end(),
                [](char& c) { c = std::tolower(c); });

  // Unit string empty? Than use default unit.
  if (unit_string.empty()) {
    switch (default_unit) {
      case Unit::kKiloJansky:
        return val * 1000.0;
      case Unit::kJansky:
        return val;
      case Unit::kMilliJansky:
        return val / 1e3;
      case Unit::kMicroJansky:
        return val / 1e6;
      case Unit::kNanoJansky:
        return val / 1e9;
    }
  } else if (unit_string_lower == "jy" || unit_string_lower == "jansky")
    return val;
  else if (unit_string == "mjy" || unit_string == "mJy" ||
           unit_string_lower == "millijansky")
    return val * 1e-3;
  else if (unit_string == "KJy" || unit_string == "kjy" ||
           unit_string_lower == "kilojansky")
    return val * 1e3;
  else if (unit_string == "µJy" || unit_string == "µjy" ||
           unit_string_lower == "microjansky")
    return val * 1e-6;
  else if (unit_string == "njy" || unit_string == "nJy" ||
           unit_string_lower == "nanojansky")
    return val * 1e-9;

  throw std::runtime_error(
      "Invalid unit specification in flux density given for " +
      value_description);
}

inline size_t FluxDensity::FindNumberEnd(const std::string& s) {
  const char* c = s.c_str();
  size_t pos = 0;
  while (IsWhitespace(c[pos])) ++pos;
  if (c[pos] == '-') ++pos;
  while (IsDigit(c[pos])) ++pos;
  if (c[pos] == '.') ++pos;
  while (IsDigit(c[pos])) ++pos;
  if (c[pos] == 'e' || c[pos] == 'E') {
    ++pos;
    if (c[pos] == '-' || c[pos] == '+') ++pos;
    while (IsDigit(c[pos])) ++pos;
  }
  return pos;
}
}  // namespace units
}  // namespace aocommon

#endif
