#ifndef AOCOMMON_UNITS_ANGLE_H_
#define AOCOMMON_UNITS_ANGLE_H_

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace aocommon {
namespace units {
class Angle {
 public:
  enum Unit { kRadians, kDegrees, kArcminutes, kArcseconds, kMilliarcseconds };
  /**
   * Parse the string as an angle, possibly with unit specification, and return
   * in radians.
   * @return The angle (rad)
   */
  static double Parse(const std::string& s,
                      const std::string& value_description,
                      Unit default_unit = Unit::kRadians);

  static std::string ToNiceString(double angle_rad);

 private:
  static size_t FindNumberEnd(const std::string& s);
  static bool IsDigit(const char c) { return c >= '0' && c <= '9'; }
  static bool IsWhitespace(const char c) { return c == ' ' || c == '\t'; }
};

inline std::string Angle::ToNiceString(double angle_rad) {
  std::ostringstream str;
  double angle_deg = angle_rad * 180.0 / M_PI;
  if (angle_deg < 0.0) {
    str << "-";
    angle_deg = -angle_deg;
  }
  if (angle_deg >= 2.0) {
    str << std::round(angle_deg * 100.0) / 100.0 << " deg";
  } else {
    const double angle_arcmin = angle_rad * 180.0 * 60.0 / M_PI;
    if (angle_arcmin >= 2.0) {
      str << std::round(angle_arcmin * 100.0) / 100.0 << "'";
    } else {
      const double angle_arcsec = angle_rad * 180.0 * 60.0 * 60.0 / M_PI;
      if (angle_arcsec >= 1.0) {
        str << std::round(angle_arcsec * 100.0) / 100.0 << "''";
      } else {
        str << std::round(angle_arcsec * 100.0 * 1000.0) / 100.0 << " masec";
      }
    }
  }
  return str.str();
}

inline double Angle::Parse(const std::string& s,
                           const std::string& value_description,
                           Unit default_unit) {
  size_t end = FindNumberEnd(s);
  if (end == 0) throw std::runtime_error("Error parsing " + value_description);
  std::string number = s.substr(0, end);
  const double val = std::atof(number.c_str());
  // Skip whitespace after number
  const char* c = s.c_str();
  while (IsWhitespace(c[end])) ++end;
  std::string unit_string = std::string(&c[end]);
  std::for_each(unit_string.begin(), unit_string.end(),
                [](char& c) { c = std::tolower(c); });

  // Unit string empty? Than use default unit.
  if (unit_string.empty()) {
    switch (default_unit) {
      case Unit::kRadians:
        return val;
      case Unit::kDegrees:
        return val * M_PI / 180.0;
      case Unit::kArcminutes:
        return val * M_PI / (180.0 * 60.0);
      case Unit::kArcseconds:
        return val * M_PI / (180.0 * 60.0 * 60.0);
      case Unit::kMilliarcseconds:
        return val * M_PI / (180.0 * 60.0 * 60.0 * 1000.0);
    }
  }

  // In degrees?
  else if (unit_string == "deg" || unit_string == "degrees")
    return val * M_PI / 180.0;

  // In arcmin?
  else if (unit_string == "amin" || unit_string == "arcmin" ||
           unit_string == "\'")
    return val * M_PI / (180.0 * 60.0);

  // In arcsec?
  else if (unit_string == "asec" || unit_string == "arcsec" ||
           unit_string == "\'\'")
    return val * M_PI / (180.0 * 60.0 * 60.0);

  // In marcsec?
  else if (unit_string == "mas" || unit_string == "masec" ||
           unit_string == "marcsec")
    return val * M_PI / (180.0 * 60.0 * 60.0 * 1000.0);

  // In radians
  else if (unit_string == "rad" || unit_string == "radians")
    return val;

  throw std::runtime_error("Invalid unit specification in angle given for " +
                           value_description);
}

inline size_t Angle::FindNumberEnd(const std::string& s) {
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
