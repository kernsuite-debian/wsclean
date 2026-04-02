#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>

#include "drawgaussian.h"

using schaapcommon::math::Ellipse;

BOOST_AUTO_TEST_SUITE(draw_gaussian)

BOOST_AUTO_TEST_CASE(to_lm) {
  constexpr size_t kHeight = 200;
  constexpr size_t kWidth = 500;
  aocommon::Image image(kWidth, kHeight, 0.0f);
  constexpr long double kRa = 20.0 * (M_PI / 180.0);   // 20 degrees
  constexpr long double kDec = 60.0 * (M_PI / 180.0);  // 60 degrees
  constexpr long double kPixelScaleL = 1.0 * (M_PI / 180.0 / 60.0);  // 1 amin
  constexpr long double kPixelScaleM = 1.0 * (M_PI / 180.0 / 60.0);  // 1 amin
  constexpr long double kLShift = 0.0;
  constexpr long double kMShift = 0.0;
  constexpr long double kSourceRa = kRa;
  constexpr long double kSourceDec = kDec;
  Ellipse ellipse;
  ellipse.major = kPixelScaleL * 20.0;
  ellipse.minor = kPixelScaleM * 5.0;
  // Rotate Gaussian ten degrees to the east
  ellipse.position_angle = 10.0 * (M_PI / 180.0);
  constexpr long double kFlux = 150.0;

  DrawGaussianToLm(image.Data(), kWidth, kHeight, kRa, kDec, kPixelScaleL,
                   kPixelScaleM, kLShift, kMShift, kSourceRa, kSourceDec,
                   ellipse, kFlux);

  BOOST_CHECK_CLOSE_FRACTION(image.Sum(), kFlux, 1e-5);

  constexpr size_t kCentralPixel = (kHeight / 2) * kWidth + kWidth / 2;
  BOOST_CHECK_CLOSE_FRACTION(image[kCentralPixel], 1.32381356, 1e-5);

  constexpr size_t kLeftOfCenterPixel = 90 * kWidth + 247;
  BOOST_CHECK_CLOSE_FRACTION(image[kLeftOfCenterPixel], 0.0631099567, 1e-5);

  constexpr size_t kRightOfCenterPixel = 90 * kWidth + 253;
  BOOST_CHECK_CLOSE_FRACTION(image[kRightOfCenterPixel], 0.532994568, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
