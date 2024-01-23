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

// Even though Pixels are int, let boost perform the intersection
// calculations based on floats to avoid rounding issues
BOOST_GEOMETRY_REGISTER_POINT_2D(schaapcommon::facets::Pixel, float,
                                 cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_RING(std::vector<schaapcommon::facets::Pixel>)

namespace schaapcommon {
namespace facets {

void Coord::Serialize(aocommon::SerialOStream& stream) const {
  stream.Double(ra).Double(dec);
}

void Coord::Unserialize(aocommon::SerialIStream& stream) {
  stream.Double(ra).Double(dec);
}

void Pixel::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(x).UInt32(y);
}

void Pixel::Unserialize(aocommon::SerialIStream& stream) {
  stream.UInt32(x).UInt32(y);
}

BoundingBox::BoundingBox(const std::vector<Pixel>& pixels, size_t align,
                         bool make_square) {
  if (pixels.empty()) {
    throw std::invalid_argument("Cannot create boundingbox for 0 pixels");
  }

  min_ = max_ = pixels.front();
  for (auto i = pixels.begin() + 1; i != pixels.end(); ++i) {
    min_.x = std::min(min_.x, i->x);
    max_.x = std::max(max_.x, i->x);
    min_.y = std::min(min_.y, i->y);
    max_.y = std::max(max_.y, i->y);
  }

  if (make_square) {
    const size_t width = max_.x - min_.x;
    const size_t height = max_.y - min_.y;
    if (width > height) {
      // Adapt height
      min_.y -= (width - height) / 2;
      max_.y = min_.y + width;
    } else {
      // Adapt width
      min_.x -= (height - width) / 2;
      max_.x = min_.x + height;
    }
  }

  if (align > 1) {
    const size_t width = max_.x - min_.x;
    const size_t height = max_.y - min_.y;
    const size_t align_x = width % align ? align - width % align : 0u;
    const size_t align_y = height % align ? align - height % align : 0u;
    min_.x -= align_x / 2;
    min_.y -= align_y / 2;
    max_.x += (align_x + 1) / 2;
    max_.y += (align_y + 1) / 2;
  }
}

void BoundingBox::Serialize(aocommon::SerialOStream& stream) const {
  stream.Object(min_).Object(max_);
}

void BoundingBox::Unserialize(aocommon::SerialIStream& stream) {
  stream.Object(min_).Object(max_);
}

std::vector<std::pair<int, int>> Facet::HorizontalIntersections(
    const int y_intersect) const {
  std::vector<std::pair<int, int>> result;
  if (y_intersect < min_y_ || y_intersect >= max_y_) {
    return result;
  }

  std::vector<Pixel> poly = pixels_;
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
    const Pixel& p1 = pixels_[i];
    const Pixel& p2 = pixels_[(i + 1) % pixels_.size()];
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
  using Line = boost::geometry::model::linestring<schaapcommon::facets::Pixel>;

  Line l;
  l.push_back(schaapcommon::facets::Pixel(trimmed_box_.Min().x, y_intersect));
  l.push_back(schaapcommon::facets::Pixel(trimmed_box_.Max().x, y_intersect));
  std::vector<Line> intersections;
  boost::geometry::intersection(l, poly, intersections);

  for (auto intersection : intersections) {
    result.emplace_back(intersection[0].x, intersection[1].x);
  }
#endif
  return result;
}

Pixel Facet::Centroid() const {
  using point_f =
      boost::geometry::model::point<float, 2, boost::geometry::cs::cartesian>;
  boost::geometry::model::polygon<point_f> poly;

  point_f x;
  for (const Pixel& pixel : pixels_) {
    boost::geometry::append(poly, point_f(pixel.x, pixel.y));
  }
  boost::geometry::centroid(poly, x);
  return Pixel(x.get<0>(), x.get<1>());
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
  const Pixel pad((width - trimmed_box.Width()) / 2,
                  (height - trimmed_box.Height()) / 2);

  // Create the padded, squared and aligned bounding box for the facet.
  return BoundingBox({trimmed_box.Min() - pad, trimmed_box.Max() + pad},
                     data.align, data.make_square);
}

Pixel RaDecToXY(const Facet::InitializationData& data,
                const Coord& coordinate) {
  Pixel pixel;
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

Coord XYToRaDec(const Facet::InitializationData& data, const Pixel& pixel) {
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
    std::vector<Pixel> image_box{Pixel(0, 0), Pixel(0, data.image_height),
                                 Pixel(data.image_width, data.image_height),
                                 Pixel(data.image_width, 0), Pixel(0, 0)};
    pixels_ = PolygonIntersection(pixels_, image_box);
  }

  if (!pixels_.empty()) {
    // Calculate bounding box for the pixels only, and set min_y_ and max_y_.
    BoundingBox pixel_box(pixels_);
    min_y_ = pixel_box.Min().y;
    max_y_ = pixel_box.Max().y;

    // Calculate the trimmed_box_.
    trimmed_box_ = BoundingBox({pixel_box.Min(), pixel_box.Max()}, data.align,
                               data.make_square);

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
  coords_.reserve(pixels_.size());
  for (const Pixel& pixel : pixels_) {
    coords_.push_back(XYToRaDec(data, pixel));
  }
}

std::vector<Pixel> Facet::PolygonIntersection(std::vector<Pixel> poly1,
                                              std::vector<Pixel> poly2) {
  // Make polygons clockwise and append closing point when needed.
  // This is the reason why poly1 and poly2 are passed by value.
  boost::geometry::correct(poly1);
  boost::geometry::correct(poly2);

  std::vector<std::vector<Pixel>> poly_results;
  boost::geometry::intersection<std::vector<Pixel>, std::vector<Pixel>>(
      poly1, poly2, poly_results);
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

bool Facet::Contains(const Pixel& pixel) const {
  std::vector<Pixel> polygon = pixels_;
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
      .Object(untrimmed_box_);
}

void Facet::Unserialize(aocommon::SerialIStream& stream) {
  stream.ObjectVector(coords_)
      .ObjectVector(pixels_)
      .UInt32(min_y_)
      .UInt32(max_y_)
      .Object(dir_)
      .String(direction_label_)
      .Object(trimmed_box_)
      .Object(untrimmed_box_);
}

}  // namespace facets
}  // namespace schaapcommon
