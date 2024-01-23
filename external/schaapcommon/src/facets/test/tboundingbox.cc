// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "boundingbox.h"
#include "tfacet.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::BoundingBox;
using schaapcommon::facets::PixelPosition;

BOOST_AUTO_TEST_SUITE(bounding_box)

BOOST_AUTO_TEST_CASE(pixel_operators) {
  const PixelPosition p1(4, 2);
  const PixelPosition p2(2, 1);

  BOOST_CHECK(p1 != p2);
  const PixelPosition p3 = p2 + p2;
  BOOST_CHECK(p1 == p3);
  const PixelPosition p4 = p3 - p2;
  BOOST_CHECK(p4 == p2);
}

BOOST_AUTO_TEST_CASE(equality_and_inequality) {
  const BoundingBox a;
  const BoundingBox b({{1, 0}, {0, 1}});
  const BoundingBox c({{1, 0}, {0, 2}});
  BOOST_CHECK(a == a);
  BOOST_CHECK(!(a != a));
  BOOST_CHECK(b == b);
  BOOST_CHECK(!(b != b));
  BOOST_CHECK(!(a == b));
  BOOST_CHECK(a != b);
  BOOST_CHECK(!(b == c));
  BOOST_CHECK(b != c);
}

BOOST_AUTO_TEST_CASE(bounding_box_empty) {
  const BoundingBox box;
  BOOST_CHECK_EQUAL(box.Min(), PixelPosition(0, 0));
  BOOST_CHECK_EQUAL(box.Max(), PixelPosition(0, 0));
  BOOST_CHECK_EQUAL(box.Width(), 0);
  BOOST_CHECK_EQUAL(box.Height(), 0);
  BOOST_CHECK_EQUAL(box.Centre(), PixelPosition(0, 0));
}

BOOST_AUTO_TEST_CASE(bounding_box_no_alignment) {
  const std::vector<PixelPosition> pixels{{-1, -20}, {4, 2}, {0, 5}};
  const BoundingBox box(pixels);
  BOOST_CHECK_EQUAL(box.Min(), PixelPosition(-1, -20));
  BOOST_CHECK_EQUAL(box.Max(), PixelPosition(4, 5));
  BOOST_CHECK_EQUAL(box.Width(), 5);
  BOOST_CHECK_EQUAL(box.Height(), 25);
  BOOST_CHECK_EQUAL(box.Centre(),
                    PixelPosition(1, -7));  // Centre is rounded towards 0.
}

BOOST_AUTO_TEST_CASE(bounding_box_aligned) {
  const std::vector<PixelPosition> pixels{{-1, -20}, {4, 2}, {0, 5}};
  const size_t align = 4;
  const BoundingBox box(pixels, align);
  BOOST_CHECK_EQUAL(box.Min(), PixelPosition(-2, -21));
  BOOST_CHECK_EQUAL(box.Max(), PixelPosition(6, 7));
  BOOST_CHECK_EQUAL(box.Width() % align, 0);
  BOOST_CHECK_EQUAL(box.Height() % align, 0);
  BOOST_CHECK_EQUAL(box.Centre(), PixelPosition(2, -7));
}

BOOST_AUTO_TEST_CASE(point_in_bounding_box) {
  const std::vector<PixelPosition> pixels{{0, 0}, {0, 10}, {10, 10}, {10, 0}};
  const BoundingBox box(pixels);

  BOOST_CHECK(box.Contains(PixelPosition(5, 5)));
  BOOST_CHECK(!box.Contains(PixelPosition(10, 10)));
}

BOOST_AUTO_TEST_SUITE_END()
