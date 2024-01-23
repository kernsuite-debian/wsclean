#ifndef AOCOMMON_COORDINATE_SYSTEM_H_
#define AOCOMMON_COORDINATE_SYSTEM_H_

#include <cstring>

namespace aocommon {

/**
 * The fields in a CoordinateSystem define the image coordinate system
 * as used in standard radio interferometric images.
 */
struct CoordinateSystem {
  /// Width in pixels
  std::size_t width;
  /// Height in pixels
  std::size_t height;
  /// J2000 right ascension, in radians
  double ra;
  /// J2000 declination, in radians
  double dec;
  /// Pixel size in horizontal direction, in radians
  double dl;
  /// Pixel size in vertical direction, in radians
  double dm;
  /// Tangential shift in the horizontal direction, in radians
  double l_shift;
  /// Tangential shift in the vertical direction, in radians
  double m_shift;
};

}  // namespace aocommon

#endif
