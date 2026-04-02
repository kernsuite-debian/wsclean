#ifndef AOCOMMON_IMAGE_COORDINATES_H_
#define AOCOMMON_IMAGE_COORDINATES_H_

#include <algorithm>
#include <cmath>
#include <span>
#include <type_traits>
#include <vector>

namespace aocommon {

/**
 * This class collects all the LM coordinate transform as defined in
 * Perley (1999)'s "imaging with non-coplaner arrays".
 */
class ImageCoordinates {
 public:
  template <typename T>
  static void RaDecToLM(T ra, T dec, T phaseCentreRa, T phaseCentreDec,
                        T& destL, T& destM) {
    const T deltaAlpha = ra - phaseCentreRa;
    const T sinDeltaAlpha = std::sin(deltaAlpha);
    const T cosDeltaAlpha = std::cos(deltaAlpha);
    const T sinDec = std::sin(dec);
    const T cosDec = std::cos(dec);
    const T sinDec0 = std::sin(phaseCentreDec);
    const T cosDec0 = std::cos(phaseCentreDec);

    destL = cosDec * sinDeltaAlpha;
    destM = sinDec * cosDec0 - cosDec * sinDec0 * cosDeltaAlpha;
  }

  template <typename T>
  static T RaDecToN(T ra, T dec, T phaseCentreRa, T phaseCentreDec) {
    const T cosDeltaAlpha = std::cos(ra - phaseCentreRa);
    const T sinDec = std::sin(dec);
    const T cosDec = std::cos(dec);
    const T sinDec0 = std::sin(phaseCentreDec);
    const T cosDec0 = std::cos(phaseCentreDec);

    return sinDec * sinDec0 + cosDec * cosDec0 * cosDeltaAlpha;
  }

  template <typename T>
  static void LMToRaDec(T l, T m, T phaseCentreRa, T phaseCentreDec, T& destRa,
                        T& destDec) {
    const T cosDec0 = std::cos(phaseCentreDec);
    const T sinDec0 = std::sin(phaseCentreDec);
    const T lmTerm = std::sqrt(T(1.0) - l * l - m * m);
    const T deltaAlpha = std::atan2(l, lmTerm * cosDec0 - m * sinDec0);

    destRa = deltaAlpha + phaseCentreRa;
    destDec = std::asin(m * cosDec0 + lmTerm * sinDec0);
  }

  template <typename T>
  static void XYToLM(size_t x, size_t y, T pixelSizeX, T pixelSizeY,
                     size_t width, size_t height, T& l, T& m) {
    const T midX = T(width) / 2.0;
    const T midY = T(height) / 2.0;
    l = (midX - (T)x) * pixelSizeX;
    m = ((T)y - midY) * pixelSizeY;
  }

  template <typename T>
  static void LMToXY(T l, T m, T pixelSizeX, T pixelSizeY, size_t width,
                     size_t height, int& x, int& y) {
    const T midX = T(width) / 2.0;
    const T midY = T(height) / 2.0;
    x = std::round(-l / pixelSizeX) + midX;
    y = std::round(m / pixelSizeY) + midY;
  }

  template <typename T>
  static void LMToXYfloat(T l, T m, T pixelSizeX, T pixelSizeY, size_t width,
                          size_t height, T& x, T& y) {
    const T midX = T(width) / 2.0;
    const T midY = T(height) / 2.0;
    x = -l / pixelSizeX + midX;
    y = m / pixelSizeY + midY;
  }

  template <typename T>
  static T AngularDistance(T ra1, T dec1, T ra2, T dec2) {
    const T sinDec1 = std::sin(dec1);
    const T cosDec1 = std::cos(dec1);
    const T sinDec2 = std::sin(dec2);
    const T cosDec2 = std::cos(dec2);
    const T cosVal =
        sinDec1 * sinDec2 + cosDec1 * cosDec2 * std::cos(ra1 - ra2);
    // Rounding errors sometimes cause cosVal to be slightly larger than 1,
    // which would cause a NaN return value.
    return cosVal <= 1.0 ? std::acos(cosVal) : 0.0;
  }

  /**
   * This function calculates the mean ra, dec position from a list of ra, dec
   * positions. It does this by transforming the equatorial coordinates to
   * 3-dimensional Cartesian positions on the unit sphere, then takes the
   * average and then converts back to equatorial coordinates. In the unlikely
   * case that the mean position falls in the origin of the sphere (e.g. with
   * two positions 180 degree apart), a value of ra 0, dec 0 is returned.
   * @tparam CoordinateType A std::pair or std::array of size 2 with numeric
   * values.
   */
  template <typename CoordinateType>
  static CoordinateType MeanPosition(
      std::span<CoordinateType> equatorial_coordinates) {
    using NumType =
        decltype(std::get<0>(std::remove_const_t<CoordinateType>()));
    // The x-axis is in the direction ra=0, dec=0.
    // The y-axis is in the direction ra=90, dec=0.
    // The z-axis is in the direction dec=90.
    if (equatorial_coordinates.empty())
      return CoordinateType{NumType(0.0), NumType(0.0)};
    NumType sum_x(0.0);
    NumType sum_y(0.0);
    NumType sum_z(0.0);
    for (const CoordinateType& coordinate : equatorial_coordinates) {
      sum_z += std::sin(std::get<1>(coordinate));
      const NumType cos_dec = std::cos(std::get<1>(coordinate));
      sum_x += cos_dec * std::cos(std::get<0>(coordinate));
      sum_y += cos_dec * std::sin(std::get<0>(coordinate));
    }
    const NumType ra = std::atan2(sum_y, sum_x);
    const NumType r = std::sqrt(sum_x * sum_x + sum_y * sum_y + sum_z * sum_z);
    if (r == 0.0) return CoordinateType{NumType(0.0), NumType(0.0)};
    const NumType dec = 0.5 * M_PI - std::acos(sum_z / r);
    return CoordinateType{ra, dec};
  }

  /**
   * Same as the overload with a std::span parameter, but for a vector. This is
   * for convenience: A std::vector is normally implicitly convertable to a
   * std::span, but because of the templated item type, this must be done
   * explicitly.
   */
  template <typename CoordinateType>
  static CoordinateType MeanPosition(
      const std::vector<CoordinateType>& equatorial_coordinates) {
    return MeanPosition(std::span(equatorial_coordinates));
  }

  template <typename T>
  static T MeanRA(const std::vector<T>& raValues) {
    std::vector<T> sorted(raValues);
    for (size_t i = 0; i != sorted.size(); ++i) {
      while (sorted[i] >= 2 * M_PI) sorted[i] -= 2.0 * M_PI;
      while (sorted[i] < 0.0) sorted[i] += 2.0 * M_PI;
    }
    std::sort(sorted.begin(), sorted.end());
    T gapSize = 0.0;
    T gapCentre = 0.0;
    for (size_t i = 0; i != sorted.size(); ++i) {
      double dist;
      if (i == sorted.size() - 1)
        dist = 2.0 * M_PI + sorted.front() - sorted.back();
      else
        dist = sorted[i + 1] - sorted[i];
      if (dist > gapSize) {
        gapSize = dist;
        gapCentre = sorted[i] + gapSize * 0.5;
      }
    }
    if (gapCentre >= 2.0 * M_PI) gapCentre -= 2.0 * M_PI;
    T sum = 0.0;
    for (size_t i = 0; i != sorted.size(); ++i) {
      if (sorted[i] < gapCentre)
        sum += sorted[i];
      else
        sum += sorted[i] - 2.0 * M_PI;
    }
    sum /= sorted.size();
    if (sum < 0.0)
      return sum + 2.0 * M_PI;
    else
      return sum;
  }

 private:
  ImageCoordinates() = delete;
};

}  // namespace aocommon

#endif
