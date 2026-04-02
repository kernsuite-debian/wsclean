#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>

#include "paddedconvolution.h"

using aocommon::Image;

namespace schaapcommon::math {
namespace {

void Check(size_t width, size_t height) {
  // Use a padding factor of 3, so that both padded even sizes remain even
  // and padded odd sizes remain odd.
  constexpr size_t kPadding = 3;

  Image image(width, height, 0.0f);
  const Image zero(width, height, 0.0f);

  PaddedConvolution(image, zero, width * kPadding, height * kPadding);
  BOOST_CHECK_LT(image.RMS(), 1e-5);

  Image delta(width, height, 0.0f);
  delta.Value(width / 2, height / 2) = 1.0f;

  PaddedConvolution(image, delta, width * kPadding, height * kPadding);
  BOOST_CHECK_LT(image.RMS(), 1e-5);

  for (size_t i = 0; i != image.Size(); ++i) image[i] = i + 3;
  PaddedConvolution(image, delta, width * kPadding, height * kPadding);
  for (size_t i = 0; i != image.Size(); ++i)
    BOOST_CHECK_CLOSE_FRACTION(image[i], i + 3, 1e-6);

  delta *= 2.0f;
  PaddedConvolution(image, delta, width * kPadding, height * kPadding);
  for (size_t i = 0; i != image.Size(); ++i)
    BOOST_CHECK_CLOSE_FRACTION(image[i], (i + 3) * 2, 1e-5);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(padded_convolution)

BOOST_AUTO_TEST_CASE(simple_even) { Check(10, 8); }

BOOST_AUTO_TEST_CASE(simple_odd) { Check(9, 11); }

BOOST_AUTO_TEST_CASE(full_example) {
  Image image(8, 6, 0.0f);
  image.Value(1, 3) = 3.0;
  image.Value(4, 5) = 1.5;
  image.Value(6, 5) = -1.0;

  Image psf(8, 6, 0.0f);
  psf.Value(4, 3) = 0.5;
  psf.Value(3, 3) = 3.0;
  psf.Value(5, 3) = 7.0;
  psf.Value(4, 2) = 10.5;
  psf.Value(4, 4) = 37.0;

  PaddedConvolution(image, psf, 20, 30);
  const Image reference(
      8, 6,
      {
          0.0f, 0.0f,   0.0f,  0.0f, 0.0f,   0.0f, 0.0f,   0.0f,  // row 0
          0.0f, 0.0f,   0.0f,  0.0f, 0.0f,   0.0f, 0.0f,   0.0f,  // row 1
          0.0f, 31.5f,  0.0f,  0.0f, 0.0f,   0.0f, 0.0f,   0.0f,  // row 2
          9.0f, 1.5f,   21.0f, 0.0f, 0.0f,   0.0f, 0.0f,   0.0f,  // row 3
          0.0f, 111.0f, 0.0f,  0.0f, 15.75f, 0.0f, -10.5f, 0.0f,  // row 4
          0.0f, 0.0f,   0.0f,  4.5f, 0.75f,  7.5f, -0.5f,  -7.0f  // row 5
      });
  for (size_t i = 0; i != reference.Size(); ++i) {
    if (reference[i] == 0.0f)
      BOOST_CHECK_LT(std::fabs(image[i]), 1e-5);
    else
      BOOST_CHECK_CLOSE_FRACTION(reference[i], image[i], 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace schaapcommon::math
