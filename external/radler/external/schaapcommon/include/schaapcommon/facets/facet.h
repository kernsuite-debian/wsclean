// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FACETS_FACET_H_
#define SCHAAPCOMMON_FACETS_FACET_H_

#include <aocommon/io/serialstreamfwd.h>

#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "boundingbox.h"

namespace schaapcommon {
namespace facets {

/// Structure for holding ra and dec coordinates.
struct Coord {
  constexpr Coord() = default;
  constexpr Coord(double _ra, double _dec) : ra(_ra), dec(_dec) {}

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

  double ra{0.0};
  double dec{0.0};
};

class Facet {
 public:
  /**
   * Data required for initializing a Facet.
   * It is only used during initialization and not stored inside a facet.
   * This struct avoids long argument lists in the Facet constructor.
   */
  struct InitializationData {
    /**
     * Constructor that accepts a single pixel scale and image size that
     * apply to both the x and the y dimensions.
     * @param scale Pixel resolution for the l and m directions (rad).
     * @param image_size Width and height of the main image (number of pixels).
     */
    explicit InitializationData(double scale, size_t image_size)
        : pixel_scale_x(scale),
          pixel_scale_y(scale),
          image_width(image_size),
          image_height(image_size) {}

    /**
     * Constructor that sets the pixel scale and image size, which are always
     * required.
     * @param scale_x Pixel resolution in l-direction (rad).
     * @param scale_y Pixel resolution in m-direction (rad).
     * @param width Width of main image (number of pixels).
     * @param height Height of main image (number of pixels).
     */
    explicit InitializationData(double scale_x, double scale_y, size_t width,
                                size_t height)
        : pixel_scale_x(scale_x),
          pixel_scale_y(scale_y),
          image_width(width),
          image_height(height) {}

    Coord phase_centre;    ///< Phase centre of main image.
    double pixel_scale_x;  ///< Pixel resolution in l-direction (rad).
    double pixel_scale_y;  ///< Pixel resolution in m-direction (rad).
    size_t image_width;    ///< Width of main image (number of pixels).
    size_t image_height;   ///< Height of main image (number of pixels).
    double l_shift = 0.0;  ///< Shift of phase centre in l-direction (rad)
    double m_shift = 0.0;  ///< Shift of phase centre in m-direction (rad)
    double padding =
        1.0;  ///< Padding factor for the bounding box. Should be >= 1.0.
    size_t feather_size =
        0;  ///< Number of pixels added to bounding box for feathering
    size_t align =
        1u;  ///< Bounding box alignment. Typically a small power of two.
    bool make_square = false;  ///< If true, create a square bounding box.
  };

  /**
   * Constructor that creates a facet using a coordinate list.
   *
   * It converts the ra+dec vertex coordinates into x+y pixel coordinates.
   *
   * Note that the following coordinate systems are adopted:
   * ra/dec:
   *
   *           ^ dec
   *           |
   *           |
   *   ra <----+
   *
   * x/y (image coords):
   *   y
   *   ^
   *   |
   *   |
   *   o --> x
   * where "o" is either the lower left corner of the main image, if
   * origin_at_centre = false, or the center pixel of the main image, if
   * origin_at_centre = true
   *
   * This function clips the facet to the image borders, if it falls outside
   * of the range of the image.
   *
   * Besides calculating the pixel coordinates, this function:
   *  - calculates a bounding box for the pixel coordinates using padding,
   * alignment and squaring parameters in the InitializationData argument.
   *  - computes the (ra, dec) position of the facet centroid in case the
   * direction argument is empty.
   *
   */
  explicit Facet(const InitializationData& data, std::vector<Coord> coordinates,
                 std::optional<Coord> direction = {});

  /**
   * Constructor that creates a facet using a boundingbox of an image.
   */
  explicit Facet(const InitializationData& data, const BoundingBox& box);

  /**
   * Constructor that deserializes facet data from an input stream.
   */
  explicit Facet(aocommon::SerialIStream& stream) { Unserialize(stream); }

  /**
   * @brief Computes the x-coordinates of the intersection points with
   * the polygonal shaped facet, given a specified y-coordinate.
   *
   * @param y_intersect y-coordinate for which to compute intersection points
   * @return A vector with pairs of x-coordinates for the first and second
   * intersection point, second > first. If no intersections were found, an
   * empty vector is returned.
   */
  std::vector<std::pair<int, int>> HorizontalIntersections(
      const int y_intersect) const;

  /**
   * Right ascension value that points inside this facet.
   * It is not necessarily the centroid, but rather a
   * point in the facet where e.g. a label can be placed to
   * identify the facet.
   */
  double RA() const { return dir_.ra; }

  /**
   * Declination value that points inside this facet.
   */
  double Dec() const { return dir_.dec; }

  bool Empty() const { return pixels_.empty(); }

  /**
   * @brief Compute the centroid of the facet in pixel coordinates.
   * Internally, pixel coordinates are converted to floats, to avoid rounding
   * issues in boost::geometry.
   */
  PixelPosition Centroid() const;

  /**
   * Get the ra+dec coordinates. This function is mainly for testing purposes.
   * @return The ra+dec coordinates of the facet.
   */
  const std::vector<Coord>& GetCoords() const { return coords_; }

  /**
   * Get the pixel coordinates. This function is mainly for testing purposes.
   */
  const std::vector<PixelPosition>& GetPixels() const { return pixels_; }

  const std::string& DirectionLabel() const { return direction_label_; }

  void SetRA(double dir_ra) { dir_.ra = dir_ra; }
  void SetDec(double dir_dec) { dir_.dec = dir_dec; }
  void SetDirectionLabel(const std::string& direction_label) {
    direction_label_ = direction_label;
  }

  /**
   * Get the trimmed bounding box for the facet.
   * The trimmed bounding box contains all pixels and is aligned.
   */
  const BoundingBox& GetTrimmedBoundingBox() const { return trimmed_box_; };

  /**
   * Get the untrimmed bounding box for the facet.
   * The untrimmed bounding box is calculated by applying padding to the trimmed
   * bounding box and then aligning the resulting box.
   */
  const BoundingBox& GetUntrimmedBoundingBox() const { return untrimmed_box_; };

  /**
   * Get a box that has requested extra space around all edges. This box can
   * be used to perform the convolution of the mask when feathering, and will
   * have enough space around it to avoid the feather kernel from wrapping
   * around.
   */
  const BoundingBox& GetConvolutionBox() const { return convolution_box_; };

  /**
   * @brief Calculate intersection between polygons \param poly1 and \param
   * poly2. Makes use of boost::geometry::intersection.
   *
   * @param poly1 Polygon 1 (the facet)
   * @param poly2 Polygon 2 (the full image)
   * @return A vector containing the vertices of the intersecting polygon, or an
   * empty vector if there is no overlap.
   */
  static std::vector<PixelPosition> PolygonIntersection(
      std::vector<PixelPosition> poly1, std::vector<PixelPosition> poly2);

  /**
   * Point-in-polygon: returns true if pixel in polygon (edges included), false
   * otherwise
   */
  bool Contains(const PixelPosition& pixel) const;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  std::vector<Coord> coords_;  ///< Ra+Dec coordinates of the Facet vertices.
  std::vector<PixelPosition>
      pixels_;  ///< Pixel coordinates of the Facet vertices.
  int min_y_;   ///< Minimum y coordinate for all pixels.
  int max_y_;   ///< Maximum y coordinate for all pixels.
  Coord dir_;   ///< (Custom) facet direction (ra, dec) in radians.
  std::string direction_label_;  ///< Description of the facet direction.
  BoundingBox trimmed_box_;      ///< Aligned bounding box for the pixels only.
  BoundingBox untrimmed_box_;    ///< Aligned bounding box including padding.
  /// Aligned bounding box that has extra zero-padding for the
  /// convolution wrap while feathering
  BoundingBox convolution_box_;
};

}  // namespace facets
}  // namespace schaapcommon
#endif
