// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "tfacet.h"
#include "facet.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::BoundingBox;
using schaapcommon::facets::Coord;
using schaapcommon::facets::Facet;
using schaapcommon::facets::PixelPosition;

namespace {
const double kScale = 0.001;
const size_t kImageSize = 100;

struct SquareFixture {
  SquareFixture() : data(kScale, kImageSize) {
    data.padding = 1.4;
    data.align = 4;
    data.make_square = true;
  }

  Facet::InitializationData data;
};

void CheckEncapsulates(const BoundingBox& untrimmed_box,
                       const BoundingBox& trimmed_box) {
  BOOST_CHECK_LE(untrimmed_box.Min().x, trimmed_box.Min().x);
  BOOST_CHECK_GE(untrimmed_box.Max().x, trimmed_box.Max().x);
  BOOST_CHECK_LE(untrimmed_box.Min().y, trimmed_box.Min().y);
  BOOST_CHECK_GE(untrimmed_box.Max().y, trimmed_box.Max().y);
}

void CheckEncapsulates(const BoundingBox& trimmed_box,
                       const PixelPosition& min_coord,
                       const PixelPosition& max_coord) {
  // Checks are slightly counter-intuitive, but due to the
  // coordinate system conventions. See documentation for CalculatePixels
  const int image_size =
      static_cast<int>(kImageSize);  // Use signed arithmatic.
  BOOST_CHECK_LE(trimmed_box.Min().x, -min_coord.x + image_size / 2);
  BOOST_CHECK_LE(trimmed_box.Min().y, min_coord.y + image_size / 2);
  BOOST_CHECK_GE(trimmed_box.Max().x, -max_coord.x + image_size / 2);
  BOOST_CHECK_GE(trimmed_box.Max().y, max_coord.y + image_size / 2);
}

/**
 * Creates a facet with a diamond shape.
 */
Facet CreateDiamondFacet(int offset_x = 0, int offset_y = 0,
                         bool set_custom_direction = false) {
  Facet::InitializationData data(kScale, kImageSize);
  data.phase_centre.ra = offset_x * kScale;
  data.phase_centre.dec = -offset_y * kScale;

  std::vector<Coord> coords{
      {0.02, 0.0}, {0.0, 0.01}, {-0.02, 0.0}, {0.0, -0.01}};

  std::optional<Coord> direction;
  if (set_custom_direction) {
    direction = Coord(-0.02, 0.01);
  }

  return Facet(data, std::move(coords), std::move(direction));
}

/**
 * Creates a facet with a rectangular shape.
 */
Facet CreateRectangularFacet(double padding, size_t align, bool make_square,
                             size_t feather_size) {
  Facet::InitializationData data(kScale, kImageSize);
  data.padding = padding;
  data.align = align;
  data.make_square = make_square;
  data.feather_size = feather_size;

  std::vector<Coord> coords{
      {0.01, -0.02}, {0.01, 0.02}, {-0.01, 0.02}, {-0.01, -0.02}};
  const Facet facet(data, std::move(coords));
  const std::vector<PixelPosition>& pixels = facet.GetPixels();
  BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
  BOOST_CHECK_EQUAL(pixels[0], PixelPosition(40, 30));
  BOOST_CHECK_EQUAL(pixels[1], PixelPosition(40, 70));
  BOOST_CHECK_EQUAL(pixels[2], PixelPosition(60, 70));
  BOOST_CHECK_EQUAL(pixels[3], PixelPosition(60, 30));

  return facet;
}

/**
 * Creates a concave facet at given x/y offset, with shape
 *    o   o
 *   / \ / \
 *  o   o   o
 *   \     /
 *    \   /
 *      o
 */
Facet CreateConcaveFacet(int offset_x = 0, int offset_y = 0) {
  Facet::InitializationData data(kScale, kImageSize);
  data.phase_centre.ra = offset_x * kScale;
  data.phase_centre.dec = -offset_y * kScale;

  std::vector<Coord> coords{{0.02, 0.0},   {0.0, -0.02}, {-0.02, 0.0},
                            {-0.01, 0.01}, {0.0, 0.0},   {0.01, 0.01}};
  return Facet(data, std::move(coords));
}

/**
 * Creates a pair of connected facets, to check whether the
 * calculated facet-line intersections seamlessly match across facet boundaries
 *  o--- o---------o
 *  |     \        |
 *  |   1  \   2   |
 *  |       \      |
 *  o -------o-----o
 *
 */
std::pair<Facet, Facet> CreateConnectedFacets() {
  Facet::InitializationData data(kScale, kImageSize);

  const std::vector<Coord> coords1{
      {0.05, -0.05}, {0.05, -0.043}, {0.043, -0.041}, {0.041, -0.05}};
  const std::vector<Coord> coords2{
      {0.043, -0.041}, {0.034, -0.043}, {0.034, -0.05}, {0.041, -0.05}};

  return std::make_pair(Facet(data, coords1), Facet(data, coords2));
}

void CheckBoundingBoxes(const Facet& facet, const PixelPosition& trimmed_min,
                        const PixelPosition& trimmed_max,
                        const PixelPosition& untrimmed_min,
                        const PixelPosition& untrimmed_max,
                        const PixelPosition& convolution_min,
                        const PixelPosition& convolution_max) {
  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Min(), trimmed_min);
  BOOST_CHECK_EQUAL(trimmed_box.Max(), trimmed_max);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  BOOST_CHECK_EQUAL(untrimmed_box.Min(), untrimmed_min);
  BOOST_CHECK_EQUAL(untrimmed_box.Max(), untrimmed_max);

  const BoundingBox& convolution_box = facet.GetConvolutionBox();
  BOOST_CHECK_EQUAL(convolution_box.Min(), convolution_min);
  BOOST_CHECK_EQUAL(convolution_box.Max(), convolution_max);
}

void CheckBoundingBoxes(const Facet& facet, const PixelPosition& trimmed_min,
                        const PixelPosition& trimmed_max,
                        const PixelPosition& untrimmed_min,
                        const PixelPosition& untrimmed_max) {
  CheckBoundingBoxes(facet, trimmed_min, trimmed_max, untrimmed_min,
                     untrimmed_max, trimmed_min, trimmed_max);
}

void CheckIntersections(const Facet& facet, int y,
                        const std::vector<std::pair<int, int>>& ref) {
  std::vector<std::pair<int, int>> isects = facet.HorizontalIntersections(y);
  BOOST_CHECK_EQUAL(isects.size(), ref.size());
  for (size_t i = 0; i != isects.size(); ++i) {
    BOOST_CHECK_EQUAL(isects[i].first, ref[i].first);
    BOOST_CHECK_EQUAL(isects[i].second, ref[i].second);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(facet)

BOOST_AUTO_TEST_CASE(polygon_empty_intersection) {
  const std::vector<PixelPosition> poly1{{0, 0}, {0, 10}, {10, 10}, {10, 0}};
  const std::vector<PixelPosition> poly2{
      {11, 11}, {11, 20}, {20, 20}, {20, 11}};
  BOOST_CHECK_EQUAL(Facet::PolygonIntersection(poly1, poly2).size(), 0);
}

BOOST_AUTO_TEST_CASE(polygon_one_intersection) {
  // Ordering for poly1 and poly2 is on purpose clockwise and anti-clockwise,
  // respectively
  const std::vector<PixelPosition> poly1{{0, 0}, {0, 10}, {10, 10}, {10, 0}};
  const std::vector<PixelPosition> poly2{{5, 5}, {15, 5}, {15, 15}, {5, 15}};
  std::vector<PixelPosition> poly3 = Facet::PolygonIntersection(poly1, poly2);
  BOOST_CHECK_EQUAL(poly3[0], PixelPosition(5, 10));
  BOOST_CHECK_EQUAL(poly3[1], PixelPosition(10, 10));
  BOOST_CHECK_EQUAL(poly3[2], PixelPosition(10, 5));
  BOOST_CHECK_EQUAL(poly3[3], PixelPosition(5, 5));
}

BOOST_AUTO_TEST_CASE(polygon_two_intersections) {
  const std::vector<PixelPosition> poly1{
      {0, 0}, {0, 10}, {10, 10}, {2, 5}, {10, 0}};
  const std::vector<PixelPosition> poly2{{5, 0}, {5, 10}, {15, 10}, {15, 0}};
  // Intersection would result in two polygons, which is not allowed
  BOOST_CHECK_THROW(Facet::PolygonIntersection(poly1, poly2),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(initialization_data_constructor_2_arguments) {
  const Facet::InitializationData data(kScale, kImageSize);
  BOOST_CHECK_EQUAL(data.phase_centre.ra, 0.0);
  BOOST_CHECK_EQUAL(data.phase_centre.dec, 0.0);
  BOOST_CHECK_EQUAL(data.pixel_scale_x, kScale);
  BOOST_CHECK_EQUAL(data.pixel_scale_y, kScale);
  BOOST_CHECK_EQUAL(data.image_width, kImageSize);
  BOOST_CHECK_EQUAL(data.image_height, kImageSize);
  BOOST_CHECK_EQUAL(data.l_shift, 0.0);
  BOOST_CHECK_EQUAL(data.m_shift, 0.0);
  BOOST_CHECK_EQUAL(data.padding, 1.0);
  BOOST_CHECK_EQUAL(data.align, 1);
  BOOST_CHECK_EQUAL(data.make_square, false);
}

BOOST_AUTO_TEST_CASE(initialization_data_constructor_4_arguments) {
  const Facet::InitializationData data(kScale, kScale + 1, kImageSize,
                                       kImageSize + 42);
  BOOST_CHECK_EQUAL(data.phase_centre.ra, 0.0);
  BOOST_CHECK_EQUAL(data.phase_centre.dec, 0.0);
  BOOST_CHECK_EQUAL(data.pixel_scale_x, kScale);
  BOOST_CHECK_EQUAL(data.pixel_scale_y, kScale + 1);
  BOOST_CHECK_EQUAL(data.image_width, kImageSize);
  BOOST_CHECK_EQUAL(data.image_height, kImageSize + 42);
  BOOST_CHECK_EQUAL(data.l_shift, 0.0);
  BOOST_CHECK_EQUAL(data.m_shift, 0.0);
  BOOST_CHECK_EQUAL(data.padding, 1.0);
  BOOST_CHECK_EQUAL(data.align, 1);
  BOOST_CHECK_EQUAL(data.make_square, false);
}

BOOST_AUTO_TEST_CASE(constructor_coordinates) {
  const Facet::InitializationData data(kScale, kImageSize);

  const std::vector<Coord> coordinates{
      {0.0, 0.0}, {0.0, 20 * kScale}, {10 * kScale, 0.0}};
  const Coord direction{42.0, 43.0};
  const Facet facet(data, coordinates, direction);

  BOOST_REQUIRE_EQUAL(facet.GetCoords().size(), coordinates.size());
  for (size_t i = 0; i < coordinates.size(); ++i) {
    BOOST_CHECK_EQUAL(facet.GetCoords()[i].ra, coordinates[i].ra);
    BOOST_CHECK_EQUAL(facet.GetCoords()[i].dec, coordinates[i].dec);
  }
  BOOST_CHECK_EQUAL(facet.GetPixels().size(), coordinates.size());
  BOOST_CHECK_EQUAL(facet.RA(), direction.ra);
  BOOST_CHECK_EQUAL(facet.Dec(), direction.dec);
  BOOST_CHECK(facet.DirectionLabel().empty());
}

BOOST_AUTO_TEST_CASE(constructor_boundingbox_full_image_square) {
  const size_t kImageWidth = 40;
  const size_t kImageHeight = 60;
  Facet::InitializationData data(kScale, kScale, kImageWidth, kImageHeight);
  data.phase_centre = Coord(0.42, 0.84);
  data.make_square = true;

  const BoundingBox box{{{0, 0}, {kImageWidth, kImageHeight}}};
  const Facet facet(data, box);
  CheckBoundingBoxes(facet, {-10, 0}, {50, 60}, {-10, 0}, {50, 60});

  const std::vector<PixelPosition> kExpectedPixels{
      {0, 0}, {0, kImageHeight}, {kImageWidth, kImageHeight}, {kImageWidth, 0}};

  // Calculate exact ra+dec coordinates for the pixels.
  std::vector<Coord> expected_coordinates;
  for (const PixelPosition& pixel : kExpectedPixels) {
    double l;
    double m;
    double ra;
    double dec;
    aocommon::ImageCoordinates::XYToLM(pixel.x, pixel.y, kScale, kScale,
                                       kImageWidth, kImageHeight, l, m);
    aocommon::ImageCoordinates::LMToRaDec(l, m, data.phase_centre.ra,
                                          data.phase_centre.dec, ra, dec);
    expected_coordinates.emplace_back(ra, dec);
  }

  BOOST_REQUIRE_EQUAL(facet.GetPixels().size(), kExpectedPixels.size());
  BOOST_REQUIRE_EQUAL(facet.GetCoords().size(), kExpectedPixels.size());
  for (size_t i = 0; i < kExpectedPixels.size(); ++i) {
    BOOST_CHECK_EQUAL(facet.GetPixels()[i], kExpectedPixels[i]);
    BOOST_CHECK_CLOSE(facet.GetCoords()[i].ra, expected_coordinates[i].ra,
                      1.0e-9);
    BOOST_CHECK_CLOSE(facet.GetCoords()[i].dec, expected_coordinates[i].dec,
                      1.0e-9);
  }
  BOOST_CHECK_CLOSE(facet.RA(), data.phase_centre.ra, 1.0e-9);
  BOOST_CHECK_CLOSE(facet.Dec(), data.phase_centre.dec, 1.0e-9);

  BOOST_CHECK(facet.DirectionLabel().empty());
}

BOOST_AUTO_TEST_CASE(constructor_boundingbox_partial_image) {
  Facet::InitializationData data(kScale, kImageSize);
  const BoundingBox box{{{20, 70}, {50, 50}}};
  const Facet facet(data, box);

  const std::vector<PixelPosition> kExpectedPixels{
      {20, 50}, {20, 70}, {50, 70}, {50, 50}};
  BOOST_REQUIRE_EQUAL(facet.GetPixels().size(), kExpectedPixels.size());
  for (size_t i = 0; i < kExpectedPixels.size(); ++i) {
    BOOST_CHECK_EQUAL(facet.GetPixels()[i], kExpectedPixels[i]);
  }

  // This test uses approximate values for the the left ra and bottom dec
  // values. The full image test above already checks that the Facet uses exact
  // coordinates. The right ra and top dec values are at the center of the
  // image, so those values should be exactly equal to 0.0.
  const double kLeft = kScale * (50 - 20);
  const double kTop = kScale * (70 - 50);
  BOOST_REQUIRE_EQUAL(facet.GetCoords().size(), kExpectedPixels.size());
  BOOST_CHECK_CLOSE(facet.GetCoords()[0].ra, kLeft, 0.1);
  BOOST_CHECK_CLOSE(facet.GetCoords()[1].ra, kLeft, 0.1);
  BOOST_CHECK_CLOSE(facet.GetCoords()[2].ra, 0.0, 1.0e-9);
  BOOST_CHECK_CLOSE(facet.GetCoords()[3].ra, 0.0, 1.0e-9);
  BOOST_CHECK_CLOSE(facet.GetCoords()[0].dec, 0.0, 1.0e-9);
  BOOST_CHECK_CLOSE(facet.GetCoords()[1].dec, kTop, 0.1);
  BOOST_CHECK_CLOSE(facet.GetCoords()[2].dec, kTop, 0.1);
  BOOST_CHECK_CLOSE(facet.GetCoords()[3].dec, 0.0, 1.0e-9);

  // Use exact coordinates for the center pixel in this test.
  double ra;
  double dec;
  double l;
  double m;
  aocommon::ImageCoordinates::XYToLM(35, 60, kScale, kScale, kImageSize,
                                     kImageSize, l, m);
  aocommon::ImageCoordinates::LMToRaDec(l, m, 0.0, 0.0, ra, dec);
  BOOST_CHECK_CLOSE(facet.RA(), ra, 1.0e-9);
  BOOST_CHECK_CLOSE(facet.Dec(), dec, 1.0e-9);

  BOOST_CHECK(facet.DirectionLabel().empty());
}

BOOST_AUTO_TEST_CASE(constructor_boundingbox_shift_align_pad) {
  // Use different x and y pixel scales in this test.
  const double kScaleX = kScale;
  const double kScaleY = kScale * 1.2;
  Facet::InitializationData data(kScaleX, kScaleY, kImageSize, kImageSize);
  data.l_shift = -0.01;
  data.m_shift = -0.02;
  data.align = 4;
  data.padding = 1.2;

  const BoundingBox box{{{10, 10}, {32, 42}}};
  const Facet facet(data, box);
  CheckBoundingBoxes(facet, {9, 10}, {33, 42}, {7, 6}, {35, 46});

  const std::vector<PixelPosition> kExpectedPixels{
      {10, 10}, {10, 42}, {32, 42}, {32, 10}};

  // Calculate exact ra+dec coordinates for the pixels.
  std::vector<Coord> expected_coordinates;
  for (const PixelPosition& pixel : kExpectedPixels) {
    double l;
    double m;
    double ra;
    double dec;
    aocommon::ImageCoordinates::XYToLM(pixel.x, pixel.y, kScaleX, kScaleY,
                                       kImageSize, kImageSize, l, m);
    l += data.l_shift;
    m += data.m_shift;
    aocommon::ImageCoordinates::LMToRaDec(l, m, 0.0, 0.0, ra, dec);
    expected_coordinates.emplace_back(ra, dec);
  }

  BOOST_REQUIRE_EQUAL(facet.GetPixels().size(), kExpectedPixels.size());
  BOOST_REQUIRE_EQUAL(facet.GetCoords().size(), kExpectedPixels.size());
  for (size_t i = 0; i < kExpectedPixels.size(); ++i) {
    BOOST_CHECK_EQUAL(facet.GetPixels()[i], kExpectedPixels[i]);
    BOOST_CHECK_CLOSE(facet.GetCoords()[i].ra, expected_coordinates[i].ra,
                      1.0e-9);
    BOOST_CHECK_CLOSE(facet.GetCoords()[i].dec, expected_coordinates[i].dec,
                      1.0e-9);
  }

  {
    double ra;
    double dec;
    double l;
    double m;
    aocommon::ImageCoordinates::XYToLM(21, 26, kScaleX, kScaleY, kImageSize,
                                       kImageSize, l, m);
    l += data.l_shift;
    m += data.m_shift;
    aocommon::ImageCoordinates::LMToRaDec(l, m, 0.0, 0.0, ra, dec);
    BOOST_CHECK_CLOSE(facet.RA(), ra, 1.0e-9);
    BOOST_CHECK_CLOSE(facet.Dec(), dec, 1.0e-9);
  }

  BOOST_CHECK(facet.DirectionLabel().empty());
}

BOOST_AUTO_TEST_CASE(constructor_boundingbox_outside_image) {
  const Facet::InitializationData data(kScale, 100);
  const BoundingBox kAbove{{{0, -1}, {50, 50}}};
  const BoundingBox kBelow{{{50, 50}, {100, 101}}};
  const BoundingBox kLeft{{{-1, 50}, {50, 100}}};
  const BoundingBox kRight{{{50, 0}, {101, 50}}};
  BOOST_CHECK_THROW(Facet(data, kAbove), std::runtime_error);
  BOOST_CHECK_THROW(Facet(data, kBelow), std::runtime_error);
  BOOST_CHECK_THROW(Facet(data, kLeft), std::runtime_error);
  BOOST_CHECK_THROW(Facet(data, kRight), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(pixel_positions) {
  Facet::InitializationData data(kScale, kImageSize);

  const std::vector<Coord> coords{{0, 0}, {0, 0.02}, {0.02, 0.02}, {0.02, 0}};

  {
    const Facet facet(data, coords);
    const std::vector<PixelPosition>& pixels = facet.GetPixels();
    BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
    BOOST_CHECK_EQUAL(pixels[0], PixelPosition(50, 50));
    BOOST_CHECK_EQUAL(pixels[1], PixelPosition(50, 70));
    BOOST_CHECK_EQUAL(pixels[2], PixelPosition(30, 70));
    BOOST_CHECK_EQUAL(pixels[3], PixelPosition(30, 50));
  }

  {
    data.l_shift = 0.01;
    data.m_shift = 0.02;
    const Facet facet(data, coords);
    const std::vector<PixelPosition>& pixels = facet.GetPixels();
    BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
    BOOST_CHECK_EQUAL(pixels[0], PixelPosition(60, 30));
    BOOST_CHECK_EQUAL(pixels[1], PixelPosition(60, 50));
    BOOST_CHECK_EQUAL(pixels[2], PixelPosition(40, 50));
    BOOST_CHECK_EQUAL(pixels[3], PixelPosition(40, 30));
  }

  {
    data.phase_centre = Coord(0.02, 0.01);
    data.l_shift = data.m_shift = 0.0;
    const Facet facet(data, coords);
    const std::vector<PixelPosition>& pixels = facet.GetPixels();
    BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
    BOOST_CHECK_EQUAL(pixels[0], PixelPosition(70, 40));
    BOOST_CHECK_EQUAL(pixels[1], PixelPosition(70, 60));
    BOOST_CHECK_EQUAL(pixels[2], PixelPosition(50, 60));
    BOOST_CHECK_EQUAL(pixels[3], PixelPosition(50, 40));
  }
}

// Create a diamond facet in the center, so no clipping occurs.
BOOST_AUTO_TEST_CASE(no_clipping) {
  Facet facet = CreateDiamondFacet();

  const std::vector<PixelPosition>& pixels = facet.GetPixels();
  BOOST_REQUIRE_EQUAL(pixels.size(), 4u);
  BOOST_CHECK_EQUAL(pixels[0], PixelPosition(30, 50));
  BOOST_CHECK_EQUAL(pixels[1], PixelPosition(50, 60));
  BOOST_CHECK_EQUAL(pixels[2], PixelPosition(70, 50));
  BOOST_CHECK_EQUAL(pixels[3], PixelPosition(50, 40));

  // Check facet centroid
  const PixelPosition centroid = facet.Centroid();
  BOOST_CHECK_EQUAL(centroid, PixelPosition(50, 50));
  BOOST_CHECK_EQUAL(facet.RA(), 0.0);
  BOOST_CHECK_EQUAL(facet.Dec(), 0.0);
}

// Create a diamond in the top right corner of the image.
BOOST_AUTO_TEST_CASE(clip_top_right) {
  const int offset_x = 50;
  const int offset_y = offset_x;
  Facet facet = CreateDiamondFacet(offset_x, offset_y);

  const std::vector<PixelPosition>& pixels = facet.GetPixels();
  // Facet clipped to triangle
  BOOST_REQUIRE_EQUAL(pixels.size(), 3u);
  BOOST_CHECK_EQUAL(pixels[0], PixelPosition(100, 90));
  BOOST_CHECK_EQUAL(pixels[1], PixelPosition(80, 100));
  BOOST_CHECK_EQUAL(pixels[2], PixelPosition(100, 100));
  PixelPosition centroid_ref((pixels[0].x + pixels[1].x + pixels[2].x) / 3,
                             (pixels[0].y + pixels[1].y + pixels[2].y) / 3);
  BOOST_CHECK_EQUAL(facet.Centroid(), centroid_ref);

  // Manually convert centroid to ra,dec coords, and check
  const double phase_centre_ra = offset_x * kScale;
  const double phase_centre_dec = -offset_y * kScale;
  double l;
  double m;
  double ra;
  double dec;
  aocommon::ImageCoordinates::XYToLM(facet.Centroid().x, facet.Centroid().y,
                                     kScale, kScale, kImageSize, kImageSize, l,
                                     m);

  aocommon::ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec,
                                        ra, dec);
  BOOST_CHECK_CLOSE(facet.RA(), ra, 1e-6);
  BOOST_CHECK_CLOSE(facet.Dec(), dec, 1e-6);
}

// Create a diamond in the bottom left corner of the image.
BOOST_AUTO_TEST_CASE(clip_bottom_left) {
  const bool set_custom_direction = true;
  Facet facet = CreateDiamondFacet(-50, -50, set_custom_direction);

  const std::vector<PixelPosition>& pixels = facet.GetPixels();
  // Facet clipped to triangle
  BOOST_REQUIRE_EQUAL(pixels.size(), 3u);
  BOOST_CHECK_EQUAL(pixels[0], PixelPosition(0, 10));
  BOOST_CHECK_EQUAL(pixels[1], PixelPosition(20, 0));
  BOOST_CHECK_EQUAL(pixels[2], PixelPosition(0, 0));
  PixelPosition centroid_ref((pixels[0].x + pixels[1].x + pixels[2].x) / 3,
                             (pixels[0].y + pixels[1].y + pixels[2].y) / 3);
  BOOST_CHECK_EQUAL(facet.Centroid(), centroid_ref);
  // Custom (ra, dec) direction reproduced?
  BOOST_CHECK_EQUAL(facet.RA(), -0.02);
  BOOST_CHECK_EQUAL(facet.Dec(), 0.01);
}

BOOST_AUTO_TEST_CASE(point_in_polygon) {
  const Facet facet = CreateDiamondFacet();

  BOOST_CHECK(facet.Contains(PixelPosition(50, 50)));
  // Pixel on edge that is owned by the facet
  BOOST_CHECK(facet.Contains(PixelPosition(30, 50)));
  // Pixel on edge that is not owned by the facet
  BOOST_CHECK(!facet.Contains(PixelPosition(50, 60)));
  BOOST_CHECK(!facet.Contains(PixelPosition(70, 50)));
  // Pixel outside facet
  BOOST_CHECK(!facet.Contains(PixelPosition(29, 50)));
  BOOST_CHECK(!facet.Contains(PixelPosition(50, 61)));
}

BOOST_AUTO_TEST_CASE(facet_bounding_boxes) {
  // Invalid padding
  BOOST_CHECK_THROW(CreateRectangularFacet(0.99, 1, false, 0),
                    std::invalid_argument);

  // Invalid alignment
  BOOST_CHECK_THROW(CreateRectangularFacet(1.0, 42, false, 0),
                    std::invalid_argument);

  // No padding, no alignment, no squaring.
  CheckBoundingBoxes(CreateRectangularFacet(1.0, 1, false, 0), {40, 30},
                     {60, 70}, {40, 30}, {60, 70});

  // Only enable padding.
  CheckBoundingBoxes(CreateRectangularFacet(1.5, 1, false, 0), {40, 30},
                     {60, 70}, {35, 20}, {65, 80});

  // Enable alignment, facet should not change
  CheckBoundingBoxes(CreateRectangularFacet(1.0, 4, false, 0), {40, 30},
                     {60, 70}, {40, 30}, {60, 70});

  // Only enable squaring.
  CheckBoundingBoxes(CreateRectangularFacet(1.0, 1, true, 0), {30, 30},
                     {70, 70}, {30, 30}, {70, 70});

  // Enable everything and use a non-power-of-two alignment.
  CheckBoundingBoxes(CreateRectangularFacet(1.5, 25, true, 0), {25, 25},
                     {75, 75}, {13, 13}, {88, 88});
}

BOOST_AUTO_TEST_CASE(feather_space) {
  // Add 5 extra pixels (everything fits inside the main image)
  CheckBoundingBoxes(CreateRectangularFacet(1.0, 1, false, 5), {35, 25},
                     {65, 75}, {35, 25}, {65, 75});

  // Add 35 extra pixels (for the trimmed bounding box: make sure space
  // is not added beyond the full image)
  CheckBoundingBoxes(CreateRectangularFacet(1.0, 1, false, 35), {5, 0},
                     {95, 100}, {5, 0}, {95, 100}, {5, -5}, {95, 105});

  // Add 5 extra pixels with factor 1.5 padding. Padding should only change
  // the untrimmed bounding box.
  CheckBoundingBoxes(CreateRectangularFacet(1.5, 1, false, 5), {35, 25},
                     {65, 75}, {28, 13}, {72, 87}, {35, 25}, {65, 75});

  // Add 5 extra pixels and request square bounding box.
  CheckBoundingBoxes(CreateRectangularFacet(1.5, 1, true, 5), {25, 25},
                     {75, 75}, {13, 13}, {87, 87}, {25, 25}, {75, 75});
}

BOOST_AUTO_TEST_CASE(horizontal_intersections_rectangle) {
  const Facet facet = CreateRectangularFacet(1.0, 1, false, 0);

  // Specified y is below the facet.
  CheckIntersections(facet, 20, {});

  // Specified y is at the bottom of the facet or intersects the facet.
  for (int y = 30; y < 70; ++y) {
    CheckIntersections(facet, y, {std::make_pair(40, 60)});
  }

  // Specified y is at the top of the facet.
  // Since the facet ranges *until* the top, the result should be empty.
  CheckIntersections(facet, 70, {});

  // Specified y is above the facet.
  CheckIntersections(facet, 90, {});
}

BOOST_AUTO_TEST_CASE(horizontal_intersections_diamond) {
  // See the clip_large_box test for the pixel coordinates of the diamond.
  const Facet facet = CreateDiamondFacet();

  // Specified y is below the facet.
  CheckIntersections(facet, 20, {});

  // Specified y is at the bottom of the facet.
  // Result is empty (0, 0), since only valid intersections are half-open
  // intervals.
  CheckIntersections(facet, 40, {});

  // Specified y intersects at 1/4 of the facet height.
  CheckIntersections(facet, 45, {std::make_pair(40, 60)});

  // Specified y intersects middle of the facet.
  CheckIntersections(facet, 50, {std::make_pair(30, 70)});

  // Specified y intersects at 3/4 of the facet height.
  CheckIntersections(facet, 55, {std::make_pair(40, 60)});

  // Specified y is at the top of the facet.
  // Since the facet ranges *until* the top, the result should be empty.
  CheckIntersections(facet, 60, {});

  // Specified y is above the facet.
  CheckIntersections(facet, 90, {});
}

BOOST_AUTO_TEST_CASE(horizontal_intersections_clipped_facet) {
  // Create diamond facet in lower left corner, with pixels extending
  // beyond image boundaries, resulting in a clipped, triangular facet
  Facet facet = CreateDiamondFacet(-50, -50);

  // Specified y is at the top of the facet.
  // Since the facet ranges *until* the top, the result should be empty.
  CheckIntersections(facet, 10, {});

  // Specified y is at the bottom of the facet.
  CheckIntersections(facet, 0, {std::make_pair(0, 20)});

  // Specified y intersects at half the facet height
  CheckIntersections(facet, 5, {std::make_pair(0, 10)});
}

BOOST_AUTO_TEST_CASE(horizontal_intersections_concave_facet) {
  Facet facet = CreateConcaveFacet();
#ifdef HAVE_BOOST_LT_166
  // Intersections for concave facet not supported
  BOOST_CHECK_THROW(CheckIntersections(facet, 30, {}), std::runtime_error);
#else
  // Specified y is at bottom of facet
  CheckIntersections(facet, 30, {});

  // Specified y cuts halfway through "convex" part
  CheckIntersections(facet, 40, {std::make_pair(40, 60)});

  // Just before "concavity" starts
  CheckIntersections(facet, 50, {std::make_pair(30, 70)});

  // Halfway "concave facets"
  CheckIntersections(facet, 55,
                     {std::make_pair(35, 45), std::make_pair(55, 65)});

  // Just below max-y vertices
  CheckIntersections(facet, 59,
                     {std::make_pair(39, 41), std::make_pair(59, 61)});

  // Empty at top
  CheckIntersections(facet, 60, {});
  BOOST_CHECK_EQUAL(facet.Centroid(), PixelPosition(50, 46));
#endif
}

BOOST_AUTO_TEST_CASE(facet_continuity) {
  std::pair<Facet, Facet> facets = CreateConnectedFacets();

  for (int y = 0; y != 8; ++y) {
    std::vector<std::pair<int, int>> isects1 =
        facets.first.HorizontalIntersections(y);
    std::vector<std::pair<int, int>> isects2 =
        facets.second.HorizontalIntersections(y);

    BOOST_CHECK_EQUAL(isects1[0].first, 0u);
    BOOST_CHECK_EQUAL(isects2[0].second, 16u);
    // Checks continuity between facets
    BOOST_CHECK_EQUAL(isects1[0].second, isects2[0].first);
  }
}

BOOST_FIXTURE_TEST_CASE(square_bounding_box, SquareFixture) {
  /**
   * These and the next 4 auto test cases check the make_square=true option. In
   * particular they make sure that if a bounding box requires extra padding to
   * make it square-shaped whether a) the bounding boxes still really
   * encapsulate all coordinates; b) the bounding boxes are still properly
   * aligned; and c) the bounding boxes are square shaped.
   */
  const PixelPosition min_coord(0, 2);
  const PixelPosition mid_coord(5, 5);
  const PixelPosition max_coord(10, 8);
  std::vector<Coord> coords{{min_coord.x * kScale, min_coord.y * kScale},
                            {mid_coord.x * kScale, mid_coord.y * kScale},
                            {max_coord.x * kScale, max_coord.y * kScale}};
  const Facet facet(data, std::move(coords));

  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Width(), trimmed_box.Height());
  CheckEncapsulates(trimmed_box, min_coord, max_coord);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  CheckEncapsulates(untrimmed_box, trimmed_box);
  BOOST_CHECK_EQUAL(untrimmed_box.Width(), untrimmed_box.Height());
  BOOST_CHECK_EQUAL(untrimmed_box.Width() % 4, 0);
}

BOOST_FIXTURE_TEST_CASE(square_bounding_box_near_right_edge, SquareFixture) {
  const PixelPosition min_coord(-49, -49);
  const PixelPosition max_coord(-48, 2);
  std::vector<Coord> coords{{min_coord.x * kScale, min_coord.y * kScale},
                            {min_coord.x * kScale, max_coord.y * kScale},
                            {max_coord.x * kScale, max_coord.y * kScale}};
  const Facet facet(data, std::move(coords));

  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Width(), trimmed_box.Height());
  CheckEncapsulates(trimmed_box, min_coord, max_coord);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  CheckEncapsulates(untrimmed_box, trimmed_box);
  BOOST_CHECK_EQUAL(untrimmed_box.Width(), untrimmed_box.Height());
  BOOST_CHECK_EQUAL(untrimmed_box.Width() % 4, 0);
}

BOOST_FIXTURE_TEST_CASE(square_bounding_box_near_bottom_edge, SquareFixture) {
  const PixelPosition min_coord(2, -49);
  const PixelPosition max_coord(49, -48);
  std::vector<Coord> coords{{min_coord.x * kScale, min_coord.y * kScale},
                            {max_coord.x * kScale, min_coord.y * kScale},
                            {max_coord.x * kScale, max_coord.y * kScale}};
  const Facet facet(data, std::move(coords));

  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Width(), trimmed_box.Height());
  CheckEncapsulates(trimmed_box, min_coord, max_coord);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  CheckEncapsulates(untrimmed_box, trimmed_box);
  BOOST_CHECK_EQUAL(untrimmed_box.Width(), untrimmed_box.Height());
  BOOST_CHECK_EQUAL(untrimmed_box.Width() % 4, 0);
}

BOOST_FIXTURE_TEST_CASE(square_bounding_box_near_left_edge, SquareFixture) {
  const PixelPosition min_coord(48, 2);
  const PixelPosition max_coord(50, 49);
  std::vector<Coord> coords{{min_coord.x * kScale, min_coord.y * kScale},
                            {max_coord.x * kScale, min_coord.y * kScale},
                            {max_coord.x * kScale, max_coord.y * kScale}};
  const Facet facet(data, std::move(coords));

  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Width(), trimmed_box.Height());
  CheckEncapsulates(trimmed_box, min_coord, max_coord);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  CheckEncapsulates(untrimmed_box, trimmed_box);
  BOOST_CHECK_EQUAL(untrimmed_box.Width(), untrimmed_box.Height());
  BOOST_CHECK_EQUAL(untrimmed_box.Width() % 4, 0);
}

BOOST_FIXTURE_TEST_CASE(square_bounding_box_near_top_edge, SquareFixture) {
  const PixelPosition min_coord(2, 48);
  const PixelPosition max_coord(49, 50);
  std::vector<Coord> coords{{min_coord.x * kScale, min_coord.y * kScale},
                            {min_coord.x * kScale, max_coord.y * kScale},
                            {max_coord.x * kScale, max_coord.y * kScale}};
  const Facet facet(data, std::move(coords));

  const BoundingBox& trimmed_box = facet.GetTrimmedBoundingBox();
  BOOST_CHECK_EQUAL(trimmed_box.Width(), trimmed_box.Height());
  CheckEncapsulates(trimmed_box, min_coord, max_coord);

  const BoundingBox& untrimmed_box = facet.GetUntrimmedBoundingBox();
  CheckEncapsulates(untrimmed_box, trimmed_box);
  BOOST_CHECK_EQUAL(untrimmed_box.Width(), untrimmed_box.Height());
  BOOST_CHECK_EQUAL(untrimmed_box.Width() % 4, 0);
}

BOOST_AUTO_TEST_CASE(serialization) {
  // Since a facet contains Coord, Boundingbox and Pixel objects, this test
  // covers serialization of those objects, too.

  const double kRa = 42.0;
  const double kDec = 43.0;
  const std::string kDirection = "FirstLeftThenRight";

  // Create a facet with different bounding boxes. See facet_bounding_boxes
  // test.
  Facet input = CreateRectangularFacet(1.5, 25, true, 30);
  BOOST_CHECK_EQUAL(input.RA(), 0.0);
  BOOST_CHECK_EQUAL(input.Dec(), 0.0);
  input.SetRA(kRa);
  input.SetDec(kDec);
  input.SetDirectionLabel(kDirection);
  aocommon::SerialOStream ostr;
  input.Serialize(ostr);

  aocommon::SerialIStream istr(std::move(ostr));
  const Facet output(istr);

  BOOST_REQUIRE_EQUAL(output.GetCoords().size(), 4u);
  BOOST_REQUIRE_EQUAL(output.GetPixels().size(), 4u);
  for (size_t i = 0; i < input.GetCoords().size(); ++i) {
    BOOST_CHECK_EQUAL(input.GetCoords()[i].ra, output.GetCoords()[i].ra);
    BOOST_CHECK_EQUAL(input.GetCoords()[i].dec, output.GetCoords()[i].dec);
    BOOST_CHECK_EQUAL(input.GetPixels()[i], output.GetPixels()[i]);
  }
  BOOST_CHECK_EQUAL(output.RA(), kRa);
  BOOST_CHECK_EQUAL(output.Dec(), kDec);
  BOOST_CHECK_EQUAL(input.Centroid(), output.Centroid());
  BOOST_CHECK_EQUAL(output.DirectionLabel(), kDirection);
  BOOST_CHECK_EQUAL(input.GetTrimmedBoundingBox().Min(),
                    output.GetTrimmedBoundingBox().Min());
  BOOST_CHECK_EQUAL(input.GetTrimmedBoundingBox().Max(),
                    output.GetTrimmedBoundingBox().Max());
  BOOST_CHECK_EQUAL(input.GetUntrimmedBoundingBox().Min(),
                    output.GetUntrimmedBoundingBox().Min());
  BOOST_CHECK_EQUAL(input.GetUntrimmedBoundingBox().Max(),
                    output.GetUntrimmedBoundingBox().Max());
}

BOOST_AUTO_TEST_SUITE_END()
