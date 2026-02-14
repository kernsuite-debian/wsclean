// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "facet.h"

#include <aocommon/imagecoordinates.h>
#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

#include <cassert>
#include <cmath>
#include <stdexcept>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/algorithms/is_convex.hpp>

// Even though PixelPositions are int, let boost perform the intersection
// calculations based on floats to avoid rounding issues
BOOST_GEOMETRY_REGISTER_POINT_2D(schaapcommon::facets::PixelPosition, float,
                                 cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_RING(std::vector<schaapcommon::facets::PixelPosition>)

namespace schaapcommon {
namespace facets {

void Coord::Serialize(aocommon::SerialOStream& stream) const {
  stream.Double(ra).Double(dec);
}

void Coord::Unserialize(aocommon::SerialIStream& stream) {
  stream.Double(ra).Double(dec);
}

std::vector<std::pair<int, int>> Facet::HorizontalIntersections(
    const int y_intersect) const {
  std::vector<std::pair<int, int>> result;
  if (y_intersect < min_y_ || y_intersect >= max_y_) {
    return result;
  }

  std::vector<PixelPosition> poly = pixels_;
  boost::geometry::correct(poly);

#ifdef HAVE_BOOST_LT_166
  if (!boost::geometry::is_convex(poly))
    throw std::runtime_error(
        "Concave facets are not supported for Boost < 1.66! Make all facets "
        "convex, or use a newer version of Boost");

  bool found_x1 = false;
  int x1 = 0;
  int x2 = 0;
  size_t i;
  for (i = 0; i != pixels_.size(); ++i) {
    const PixelPosition& p1 = pixels_[i];
    const PixelPosition& p2 = pixels_[(i + 1) % pixels_.size()];
    if ((p1.y <= y_intersect && p2.y > y_intersect) ||
        (p2.y <= y_intersect && p1.y > y_intersect)) {
      int x;
      if (p1.y == y_intersect) {
        x = p1.x;
      } else if (p2.y == y_intersect) {
        x = p2.x;
      } else {
        const double beta = double(p2.x - p1.x) / double(p2.y - p1.y);
        const double xfl = p1.x + beta * (y_intersect - p1.y);
        x = round(xfl);
      }
      if (!found_x1) {
        x1 = x;
        found_x1 = true;
      } else {
        x2 = x;
        break;
      }
    }
  }

  // The loop should have found x1 and x2, and then stopped using 'break'.
  assert(i != pixels_.size());
  if (x1 != x2) result.push_back(std::minmax(x1, x2));
#else
  using Line =
      boost::geometry::model::linestring<schaapcommon::facets::PixelPosition>;

  Line l;
  l.push_back(
      schaapcommon::facets::PixelPosition(trimmed_box_.Min().x, y_intersect));
  l.push_back(
      schaapcommon::facets::PixelPosition(trimmed_box_.Max().x, y_intersect));
  std::vector<Line> intersections;
  boost::geometry::intersection(l, poly, intersections);

  for (auto intersection : intersections) {
    result.emplace_back(intersection[0].x, intersection[1].x);
  }
#endif
  return result;
}

PixelPosition Facet::Centroid() const {
  using point_f =
      boost::geometry::model::point<float, 2, boost::geometry::cs::cartesian>;
  boost::geometry::model::polygon<point_f> poly;

  for (const PixelPosition& pixel : pixels_) {
    boost::geometry::append(poly, point_f(pixel.x, pixel.y));
  }

  point_f x{};
  boost::geometry::centroid(poly, x);
  return PixelPosition(x.get<0>(), x.get<1>());
}

namespace {
void CheckData(const Facet::InitializationData& data) {
  if (data.padding < 1.0) {
    throw std::invalid_argument("Padding factor should be >= 1.0");
  }
  if ((data.align > 1) &&
      ((data.image_width % data.align) || (data.image_height % data.align))) {
    throw std::invalid_argument("Image is not aligned");
  }
}

BoundingBox CalculateUntrimmedBox(const Facet::InitializationData& data,
                                  const BoundingBox& trimmed_box) {
  auto width =
      static_cast<size_t>(std::ceil(trimmed_box.Width() * data.padding));
  auto height =
      static_cast<size_t>(std::ceil(trimmed_box.Height() * data.padding));

  // Calculate padding. Divide by two since the padding occurs on both sides.
  const PixelPosition pad((width - trimmed_box.Width()) / 2,
                          (height - trimmed_box.Height()) / 2);

  // Create the padded, squared and aligned bounding box for the facet.
  return BoundingBox({trimmed_box.Min() - pad, trimmed_box.Max() + pad},
                     data.align, data.make_square);
}

PixelPosition RaDecToXY(const Facet::InitializationData& data,
                        const Coord& coordinate) {
  PixelPosition pixel;
  double l;
  double m;
  aocommon::ImageCoordinates::RaDecToLM(coordinate.ra, coordinate.dec,
                                        data.phase_centre.ra,
                                        data.phase_centre.dec, l, m);
  l -= data.l_shift;
  m -= data.m_shift;
  aocommon::ImageCoordinates::LMToXY(l, m, data.pixel_scale_x,
                                     data.pixel_scale_y, data.image_width,
                                     data.image_height, pixel.x, pixel.y);
  return pixel;
}

Coord XYToRaDec(const Facet::InitializationData& data,
                const PixelPosition& pixel) {
  Coord coordinates;
  double l;
  double m;
  aocommon::ImageCoordinates::XYToLM(pixel.x, pixel.y, data.pixel_scale_x,
                                     data.pixel_scale_y, data.image_width,
                                     data.image_height, l, m);
  l += data.l_shift;
  m += data.m_shift;
  aocommon::ImageCoordinates::LMToRaDec(l, m, data.phase_centre.ra,
                                        data.phase_centre.dec, coordinates.ra,
                                        coordinates.dec);
  return coordinates;
}
}  // namespace

Facet::Facet(const InitializationData& data, std::vector<Coord> coordinates,
             std::optional<Coord> direction)
    : coords_(std::move(coordinates)),
      pixels_(),
      min_y_(0),
      max_y_(0),
      dir_(),
      trimmed_box_(),
      untrimmed_box_() {
  CheckData(data);
  if (coords_.size() < 3) {
    throw std::runtime_error("Number of coordinates < 3, facet incomplete!");
  }

  pixels_.reserve(coords_.size());
  bool need_clip = false;
  for (const Coord& coordinate : coords_) {
    pixels_.push_back(RaDecToXY(data, coordinate));
    const int x = pixels_.back().x;
    const int y = pixels_.back().y;
    if (!need_clip && (x < 0 || x > static_cast<int>(data.image_width) ||
                       y < 0 || y > static_cast<int>(data.image_height))) {
      need_clip = true;
    }
  }

  if (need_clip) {
    std::vector<PixelPosition> image_box{
        PixelPosition(0, 0), PixelPosition(0, data.image_height),
        PixelPosition(data.image_width, data.image_height),
        PixelPosition(data.image_width, 0), PixelPosition(0, 0)};
    pixels_ = PolygonIntersection(pixels_, image_box);
  }

  if (!pixels_.empty()) {
    // Calculate bounding box for the pixels only, and set min_y_ and max_y_.
    const BoundingBox pixel_box(pixels_);
    min_y_ = pixel_box.Min().y;
    max_y_ = pixel_box.Max().y;

    // Calculate the trimmed_box_.
    const std::vector<PixelPosition> boundary_positions{pixel_box.Min(),
                                                        pixel_box.Max()};
    trimmed_box_ = BoundingBox(boundary_positions, data.align, data.make_square,
                               data.feather_size, data.image_width,
                               data.image_height, false);

    convolution_box_ = BoundingBox(boundary_positions, data.align,
                                   data.make_square, data.feather_size,
                                   data.image_width, data.image_height, true);

    untrimmed_box_ = CalculateUntrimmedBox(data, trimmed_box_);

    dir_ = direction.value_or(XYToRaDec(data, Centroid()));
  }
}

Facet::Facet(const InitializationData& data, const BoundingBox& box)
    : coords_(),
      pixels_{box.Min(),
              {box.Min().x, box.Max().y},
              box.Max(),
              {box.Max().x, box.Min().y}},
      min_y_(box.Min().y),
      max_y_(box.Max().y),
      dir_(XYToRaDec(data, Centroid())),
      direction_label_(),
      trimmed_box_(box),
      untrimmed_box_(CalculateUntrimmedBox(data, trimmed_box_)) {
  CheckData(data);
  if (trimmed_box_.Min().x < 0 || trimmed_box_.Min().y < 0 ||
      trimmed_box_.Max().x > static_cast<int>(data.image_width) ||
      trimmed_box_.Max().y > static_cast<int>(data.image_height)) {
    throw std::runtime_error(
        "Facet bounding box extends beyond image boundaries");
  }
  if (data.make_square || data.align > 1) {
    trimmed_box_ =
        BoundingBox({box.Min(), box.Max()}, data.align, data.make_square);
  }
  convolution_box_ = trimmed_box_;
  coords_.reserve(pixels_.size());
  for (const PixelPosition& pixel : pixels_) {
    coords_.push_back(XYToRaDec(data, pixel));
  }
}

std::vector<PixelPosition> Facet::PolygonIntersection(
    std::vector<PixelPosition> poly1, std::vector<PixelPosition> poly2) {
  // Make polygons clockwise and append closing point when needed.
  // This is the reason why poly1 and poly2 are passed by value.
  boost::geometry::correct(poly1);
  boost::geometry::correct(poly2);

  std::vector<std::vector<PixelPosition>> poly_results;
  boost::geometry::intersection<std::vector<PixelPosition>,
                                std::vector<PixelPosition>>(poly1, poly2,
                                                            poly_results);
  if (poly_results.size() == 1) {
    // Return intersection points, except for closing point.
    poly_results.front().resize(poly_results.front().size() - 1);
    return std::move(poly_results.front());
  } else if (poly_results.empty()) {
    return {};
  } else {
    throw std::runtime_error(
        "Expected 0 or 1 intersecting polygons, but found " +
        std::to_string(poly_results.size()));
  }
}

bool Facet::Contains(const PixelPosition& pixel) const {
  std::vector<PixelPosition> polygon = pixels_;
  boost::geometry::correct(polygon);

  bool inClosedPolygon = boost::geometry::covered_by(pixel, polygon);

  if (!inClosedPolygon) {
    return false;
  } else {
    // Pixel is in the closed polygon, but we should exclude it from
    // zero-length corners and edges that are "owned" by another facet
    const std::vector<std::pair<int, int>> intersections =
        HorizontalIntersections(pixel.y);
    if (intersections.empty()) return false;

    for (const auto& isect : intersections) {
      if (isect.second == pixel.x) {
        return false;
      }
    }
    return true;
  }
}

void Facet::Serialize(aocommon::SerialOStream& stream) const {
  stream.ObjectVector(coords_)
      .ObjectVector(pixels_)
      .UInt32(min_y_)
      .UInt32(max_y_)
      .Object(dir_)
      .String(direction_label_)
      .Object(trimmed_box_)
      .Object(untrimmed_box_)
      .Object(convolution_box_);
}

void Facet::Unserialize(aocommon::SerialIStream& stream) {
  stream.ObjectVector(coords_)
      .ObjectVector(pixels_)
      .UInt32(min_y_)
      .UInt32(max_y_)
      .Object(dir_)
      .String(direction_label_)
      .Object(trimmed_box_)
      .Object(untrimmed_box_)
      .Object(convolution_box_);
}

}  // namespace facets
}  // namespace schaapcommon
