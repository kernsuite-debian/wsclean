// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "facetimage.h"
#include "tfacet.h"

#include <aocommon/uvector.h>

#include <iostream>

using schaapcommon::facets::BoundingBox;
using schaapcommon::facets::Coord;
using schaapcommon::facets::Facet;
using schaapcommon::facets::FacetImage;
using schaapcommon::facets::PixelPosition;

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_SUITE(facetimage)

namespace {

const double kScale = 0.001;
const size_t kImageSize = 100;
const size_t kNTerms = 2;
const double kPadding = 1.25;
const size_t kAlign = 4;

/// Creates a facet with pixels (0, 0), (size, 0), (size, size), (0, size).
Facet CreateSquareFacet(size_t size) {
  Facet::InitializationData data(kScale, kImageSize);
  data.padding = kPadding;
  data.align = kAlign;
  data.make_square = true;

  const double radecZero = kScale * kImageSize / 2;
  const double radecSize = radecZero - size * kScale;
  const std::vector<Coord> coordinates{{radecZero, -radecZero},
                                       {radecSize, -radecZero},
                                       {radecSize, -radecSize},
                                       {radecZero, -radecSize}};

  const Facet facet(data, std::move(coordinates));
  const std::vector<PixelPosition>& pixels = facet.GetPixels();
  BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
  BOOST_CHECK_EQUAL(pixels[0], PixelPosition(0, 0));
  BOOST_CHECK_EQUAL(pixels[1], PixelPosition(size, 0));
  BOOST_CHECK_EQUAL(pixels[2], PixelPosition(size, size));
  BOOST_CHECK_EQUAL(pixels[3], PixelPosition(0, size));

  return facet;
}

/**
 * Fill the data for a facet with a value.
 * Use a different value in the padding area.
 */
void FillFacetData(float* data, const Facet& facet, float value) {
  const float pad_value = -42.0f;

  const BoundingBox unpadded_box(facet.GetPixels());
  const BoundingBox& padded_box = facet.GetUntrimmedBoundingBox();

  // Check that the padded box is larger on all sides.
  BOOST_REQUIRE_LT(padded_box.Min().x, unpadded_box.Min().x);
  BOOST_REQUIRE_LT(padded_box.Min().y, unpadded_box.Min().y);
  BOOST_REQUIRE_GT(padded_box.Max().x, unpadded_box.Max().x);
  BOOST_REQUIRE_GT(padded_box.Max().y, unpadded_box.Max().y);

  const int pad_top = unpadded_box.Min().y - padded_box.Min().y;
  for (int y = 0; y < pad_top; ++y) {
    std::fill_n(data, padded_box.Width(), pad_value);
    data += padded_box.Width();
  }

  const int pad_left = unpadded_box.Min().x - padded_box.Min().x;
  const int pad_right = padded_box.Max().x - unpadded_box.Max().x;
  for (int y = unpadded_box.Min().y; y < unpadded_box.Max().y; ++y) {
    std::fill_n(data, pad_left, pad_value);
    data += pad_left;
    std::fill_n(data, unpadded_box.Width(), value);
    data += unpadded_box.Width();
    std::fill_n(data, pad_right, pad_value);
    data += pad_right;
  }

  const int pad_bottom = padded_box.Max().y - unpadded_box.Max().y;
  for (int y = 0; y < pad_bottom; ++y) {
    std::fill_n(data, padded_box.Width(), pad_value);
    data += padded_box.Width();
  }
}

std::string MaskToString(const aocommon::Image& mask) {
  std::ostringstream str;
  const float* iterator = mask.Data();
  for (size_t y = 0; y != mask.Height(); ++y) {
    for (size_t x = 0; x != mask.Width(); ++x) {
      BOOST_CHECK(*iterator == 0.0 || *iterator == 1.0);
      str << (*iterator == 0.0 ? ' ' : 'X');
      ++iterator;
    }
    str << '\n';
  }
  return str.str();
}

FacetImage MakeTriangleFacetImage() {
  Facet::InitializationData data(0.004, 0.004, 50, 50);
  data.padding = kPadding;
  data.align = 1;
  data.make_square = false;

  std::vector<Coord> coordinates{{-0.08, 0.08}, {-0.05, 0.07}, {-0.06, 0.04}};
  const Facet facet(data, std::move(coordinates));
  FacetImage facet_image(data.image_width, data.image_height, 1);
  facet_image.SetFacet(facet, false);
  return facet_image;
}

}  // namespace

BOOST_AUTO_TEST_CASE(constructor) {
  FacetImage fi(kImageSize, kImageSize, kNTerms);
  BOOST_CHECK_EQUAL(fi.Width(), 0);
  BOOST_CHECK_EQUAL(fi.Height(), 0);
  BOOST_CHECK_EQUAL(fi.OffsetX(), 0);
  BOOST_CHECK_EQUAL(fi.OffsetY(), 0);
}

BOOST_AUTO_TEST_CASE(set_facet) {
  const Facet facet = CreateSquareFacet(40);
  BOOST_CHECK_EQUAL(facet.GetTrimmedBoundingBox().Min(), PixelPosition(0, 0));
  BOOST_CHECK_EQUAL(facet.GetTrimmedBoundingBox().Max(), PixelPosition(40, 40));

  // Padding increases the facet to -5...45.
  // Alignment rounds -5 down to -6 and 45 up to 46, yielding (46 -- 6) % 4 = 0
  BOOST_CHECK_EQUAL(facet.GetUntrimmedBoundingBox().Min(),
                    PixelPosition(-6, -6));
  BOOST_CHECK_EQUAL(facet.GetUntrimmedBoundingBox().Max(),
                    PixelPosition(46, 46));

  FacetImage fi(kImageSize, kImageSize, kNTerms);
  // Facet is uninitialized upon construction
  BOOST_CHECK(!fi.GetFacet());
  // Trimmed box
  fi.SetFacet(facet, true);

  // Now the facet is initialized
  BOOST_CHECK(fi.GetFacet());

  // Repeat checks above on facet returned by FacetImage::GetFacet()
  BOOST_CHECK_EQUAL(fi.GetFacet()->GetTrimmedBoundingBox().Min(),
                    PixelPosition(0, 0));
  BOOST_CHECK_EQUAL(fi.GetFacet()->GetTrimmedBoundingBox().Max(),
                    PixelPosition(40, 40));
  BOOST_CHECK_EQUAL(fi.GetFacet()->GetUntrimmedBoundingBox().Min(),
                    PixelPosition(-6, -6));
  BOOST_CHECK_EQUAL(fi.GetFacet()->GetUntrimmedBoundingBox().Max(),
                    PixelPosition(46, 46));

  // Check image dimensions derived from facet
  BOOST_CHECK_EQUAL(fi.Width(), 40);
  BOOST_CHECK_EQUAL(fi.Height(), 40);
  BOOST_CHECK_EQUAL(fi.OffsetX(), 0);
  BOOST_CHECK_EQUAL(fi.OffsetY(), 0);

  // Untrimmed, padded box
  fi.SetFacet(facet, false);
  BOOST_CHECK_EQUAL(fi.Width(), 52);
  BOOST_CHECK_EQUAL(fi.Height(), 52);
  BOOST_CHECK_EQUAL(fi.OffsetX(), -6);
  BOOST_CHECK_EQUAL(fi.OffsetY(), -6);
}

BOOST_AUTO_TEST_CASE(copy_to_facet) {
  const size_t kFacetSize = 30;

  const Facet facet = CreateSquareFacet(kFacetSize);

  std::vector<aocommon::UVector<float>> input_images;
  for (size_t term = 0; term < kNTerms; ++term) {
    aocommon::UVector<float> image(kImageSize * kImageSize);
    for (size_t y = 0; y < kImageSize; ++y) {
      for (size_t x = 0; x < kImageSize; ++x) {
        image[y * kImageSize + x] = x + y + term + 0.5;
      }
    }
    input_images.push_back(std::move(image));
  }

  FacetImage facet_image(kImageSize, kImageSize, kNTerms);
  facet_image.SetFacet(facet, false);
  facet_image.CopyToFacet(input_images);

  for (size_t term = 0; term < kNTerms; ++term) {
    const float* data = facet_image.Data(term);
    for (int y = facet.GetUntrimmedBoundingBox().Min().y;
         y < facet.GetUntrimmedBoundingBox().Max().y; ++y) {
      for (int x = facet.GetUntrimmedBoundingBox().Min().x;
           x < facet.GetUntrimmedBoundingBox().Max().x; ++x) {
        // Values in the padding and alignment areas should be zero.
        float ref_value = 0.0;
        if (x >= 0 && x < static_cast<int>(kFacetSize) && y >= 0 &&
            y < static_cast<int>(kFacetSize)) {
          ref_value = x + y + term + 0.5;
        }
        BOOST_CHECK_CLOSE(*data, ref_value, 1e-8);
        ++data;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(multiply) {
  const size_t kFacetSize = 30;

  const Facet facet = CreateSquareFacet(kFacetSize);

  std::vector<aocommon::UVector<float>> input_images(kNTerms);
  for (size_t term = 0; term < kNTerms; ++term) {
    input_images[term].assign(kImageSize * kImageSize, term + 1.0f);
  }

  FacetImage facet_image(kImageSize, kImageSize, kNTerms);
  facet_image.SetFacet(facet, false);
  facet_image.CopyToFacet(input_images);

  BOOST_CHECK_NO_THROW(facet_image *= 2.0);

  for (size_t term = 0; term < kNTerms; ++term) {
    const float* data = facet_image.Data(term);
    for (int y = facet.GetUntrimmedBoundingBox().Min().y;
         y < facet.GetUntrimmedBoundingBox().Max().y; ++y) {
      for (int x = facet.GetUntrimmedBoundingBox().Min().x;
           x < facet.GetUntrimmedBoundingBox().Max().x; ++x) {
        // Values in the padding and alignment areas should be zero.
        float ref_value = 0.0;
        if (x >= 0 && x < static_cast<int>(kFacetSize) && y >= 0 &&
            y < static_cast<int>(kFacetSize)) {
          ref_value = 2.0f * (term + 1.0f);
        }
        BOOST_CHECK_CLOSE(*data, ref_value, 1e-8);
        ++data;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(add_to_image) {
  // Properties of the main image:
  Facet::InitializationData data(0.002, 0.002, kImageSize * 2, kImageSize);
  data.padding = kPadding;
  data.align = kAlign;
  data.make_square = false;

  // Make a facet in the lower left corner and upper right corner
  // (please note coord system conventions!)
  std::vector<Coord> coordinates1{
      {0.1, -0.05}, {0.05, -0.05}, {0.05, -0.03}, {0.1, -0.03}};
  std::vector<Coord> coordinates2;

  for (const Coord& coordinate : coordinates1) {
    // Mirror in origin
    coordinates2.emplace_back(-1 * coordinate.ra, -1 * coordinate.dec);
  }

  const Facet facet1(data, std::move(coordinates1));
  data.make_square = true;
  const Facet facet2(data, std::move(coordinates2));

  std::vector<std::vector<float>> output_images;
  std::vector<float*> output_ptrs;
  for (size_t term = 0; term < kNTerms; ++term) {
    output_images.emplace_back(data.image_width * data.image_height, 0.0f);
    output_ptrs.push_back(output_images.back().data());
  }

  FacetImage facet_image(data.image_width, data.image_height, kNTerms);

  facet_image.SetFacet(facet1, false);
  // Initialize facet1 with odd numbers.
  for (size_t term = 0; term < kNTerms; ++term) {
    FillFacetData(facet_image.Data(term), facet1, 1 + term * 2);
  }
  facet_image.AddToImage(output_ptrs);

  facet_image.SetFacet(facet2, false);
  // Initialize facet2 with even numbers.
  for (size_t term = 0; term < kNTerms; ++term) {
    FillFacetData(facet_image.Data(term), facet2, 2 + term * 2);
  }
  facet_image.AddToImage(output_ptrs);

  // Create unpadded bounding boxes for the facets: AddToImage should not use
  // data in the padding area.
  const BoundingBox box1(facet1.GetPixels());
  const BoundingBox box2(facet2.GetPixels());

  for (size_t term = 0; term < kNTerms; ++term) {
    for (int y = 0; y < static_cast<int>(data.image_height); ++y) {
      for (int x = 0; x < static_cast<int>(data.image_width); ++x) {
        const float result_pixel =
            output_images[term][y * data.image_width + x];
        if (y >= box1.Min().y && y < box1.Max().y && x >= box1.Min().x &&
            x < box1.Max().x) {
          // Check if "facet1" value was assigned.
          BOOST_CHECK_CLOSE(result_pixel, 1 + term * 2, 1e-8);
        } else if (y >= box2.Min().y && y < box2.Max().y && x >= box2.Min().x &&
                   x < box2.Max().x) {
          // Check if "facet2" value was assigned.
          BOOST_CHECK_CLOSE(result_pixel, 2 + term * 2, 1e-8);
        } else {
          // Check if the value is still zero.
          BOOST_CHECK_CLOSE(result_pixel, 0.0f, 1e-8);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(make_mask) {
  FacetImage facet_image = MakeTriangleFacetImage();
  FillFacetData(facet_image.Data(0), *facet_image.GetFacet(), 1);
  const aocommon::Image image = facet_image.MakeMask();
  BOOST_CHECK_EQUAL(image.Width(), 10);
  BOOST_CHECK_EQUAL(image.Height(), 12);
  const std::string reference =
      "          \n"
      "          \n"
      "   X      \n"
      "   XX     \n"
      "  XXX     \n"
      "  XXXX    \n"
      " XXXXX    \n"
      " XXXXXX   \n"
      " XXXXXX   \n"
      "   XXXXX  \n"
      "      XX  \n"
      "          \n";
  BOOST_CHECK_EQUAL(MaskToString(image), reference);
}

BOOST_AUTO_TEST_CASE(add_with_mask) {
  FacetImage facet_image = MakeTriangleFacetImage();
  std::fill_n(facet_image.Data(0), facet_image.Width() * facet_image.Height(),
              4.0f);

  const size_t full_width = facet_image.FullImageWidth();
  const size_t full_height = facet_image.FullImageHeight();
  aocommon::Image full_image(full_width, full_height, 3.0f);
  aocommon::Image weight_image(full_width, full_height, 10.0f);
  aocommon::Image mask = facet_image.MakeMask();
  mask *= 0.5;
  mask[0] = 1.0;  // Make one weight outside bounding box different to test this
                  // works too
  facet_image.AddWithMask(full_image.Data(), weight_image.Data(), mask, 0);
  const size_t offset =
      facet_image.OffsetX() + facet_image.OffsetY() * full_width;
  size_t image_modifications = 0;
  for (size_t y = 0; y != facet_image.Height(); ++y) {
    for (size_t x = 0; x != facet_image.Width(); ++x) {
      if (x == 0 && y == 0) {
        BOOST_CHECK_CLOSE_FRACTION(full_image[offset], 7.0f, 1e-6);
        BOOST_CHECK_CLOSE_FRACTION(weight_image[offset], 11.0f, 1e-6);
      } else {
        bool image_changed =
            std::fabs(full_image[offset + x + y * full_width] - 3.0f) > 1e-6;
        bool weight_changed =
            std::fabs(weight_image[offset + x + y * full_width] - 10.0f) > 1e-6;
        BOOST_CHECK(image_changed == weight_changed);
        if (image_changed) {
          ++image_modifications;
          BOOST_CHECK_CLOSE_FRACTION(full_image[offset + x + y * full_width],
                                     5.0f, 1e-6);
          BOOST_CHECK_CLOSE_FRACTION(weight_image[offset + x + y * full_width],
                                     10.5f, 1e-6);
        }
      }
    }
  }
  // Facet contains 34 pixels (see nr of Xes in the make_mask test)
  BOOST_CHECK_EQUAL(image_modifications, 34);
}

BOOST_AUTO_TEST_CASE(multiply_image_inside_facet) {
  // Properties of the main image:
  Facet::InitializationData data(0.002, 0.002, kImageSize * 2, kImageSize);
  data.padding = kPadding;
  data.align = kAlign;
  data.make_square = false;

  // Make a facet in the lower left corner and upper right corner
  // (please note coord system conventions!)
  std::vector<Coord> coordinates1{
      {0.1, -0.05}, {0.05, -0.05}, {0.05, -0.03}, {0.1, -0.03}};
  std::vector<Coord> coordinates2;

  for (const Coord& coordinate : coordinates1) {
    // Mirror in origin
    coordinates2.emplace_back(-1 * coordinate.ra, -1 * coordinate.dec);
  }

  const Facet facet1(data, std::move(coordinates1));
  data.make_square = true;
  const Facet facet2(data, std::move(coordinates2));

  std::vector<std::vector<float>> output_images;
  std::vector<float*> output_ptrs;
  for (size_t term = 0; term < kNTerms; ++term) {
    output_images.emplace_back(data.image_width * data.image_height, 1.0f);
    output_ptrs.push_back(output_images.back().data());
  }

  FacetImage facet_image(data.image_width, data.image_height, kNTerms);

  facet_image.SetFacet(facet1, false);
  BOOST_CHECK_NO_THROW(facet_image.MultiplyImageInsideFacet(output_ptrs, 2.0f));

  facet_image.SetFacet(facet2, false);
  BOOST_CHECK_NO_THROW(facet_image.MultiplyImageInsideFacet(output_ptrs, 4.0f));

  // Create unpadded bounding boxes for the facets: AddToImage should not use
  // data in the padding area.
  const BoundingBox box1(facet1.GetPixels());
  const BoundingBox box2(facet2.GetPixels());

  for (size_t term = 0; term < kNTerms; ++term) {
    for (int y = 0; y < static_cast<int>(data.image_height); ++y) {
      for (int x = 0; x < static_cast<int>(data.image_width); ++x) {
        const float result_pixel =
            output_images[term][y * data.image_width + x];
        if (y >= box1.Min().y && y < box1.Max().y && x >= box1.Min().x &&
            x < box1.Max().x) {
          // Check if "facet1" value was assigned.
          BOOST_CHECK_CLOSE(result_pixel, 2.0f, 1e-8);
        } else if (y >= box2.Min().y && y < box2.Max().y && x >= box2.Min().x &&
                   x < box2.Max().x) {
          // Check if "facet2" value was assigned.
          BOOST_CHECK_CLOSE(result_pixel, 4.0f, 1e-8);
        } else {
          // Check if the value is still one.
          BOOST_CHECK_CLOSE(result_pixel, 1.0f, 1e-8);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
