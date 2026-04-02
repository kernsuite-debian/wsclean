// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FACETS_BOUNDING_BOX_H_
#define SCHAAPCOMMON_FACETS_BOUNDING_BOX_H_

#include <cstddef>
#include <vector>

#include <aocommon/io/serialstreamfwd.h>

namespace schaapcommon::facets {

/// Structure for holding pixel coordinates.
struct PixelPosition {
  constexpr PixelPosition() = default;

  constexpr PixelPosition(int _x, int _y) : x(_x), y(_y) {}

  constexpr friend PixelPosition operator+(const PixelPosition& a,
                                           const PixelPosition& b) {
    return {a.x + b.x, a.y + b.y};
  }

  constexpr friend PixelPosition operator-(const PixelPosition& a,
                                           const PixelPosition& b) {
    return {a.x - b.x, a.y - b.y};
  }

  constexpr friend bool operator==(const PixelPosition& a,
                                   const PixelPosition& b) {
    return (a.x == b.x) && (a.y == b.y);
  }

  constexpr friend bool operator!=(const PixelPosition& a,
                                   const PixelPosition& b) {
    return (a.x != b.x) || (a.y != b.y);
  }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

  int x{0};
  int y{0};
};

class BoundingBox {
 public:
  BoundingBox() : min_(0, 0), max_(0, 0) {}

  /**
   * Determine bounding box from a list of pixel positions. This constructor
   * allows adding extra space for feathering or convolution. When feathering,
   * it is not necessary to add extra space outside the full image, and thus
   * the bounding box will not extend beyond the full image to circumvent
   * increasing the cost of inversion more than necessary. If however the
   * padding should extend beyond the image, @c always_pad can be set to
   * @c true. This is e.g. useful when performing a convolution and padded
   * space is required to take wrapping into account.
   */
  explicit BoundingBox(const std::vector<PixelPosition>& pixels, size_t align,
                       bool make_square, size_t extra_space, size_t full_width,
                       size_t full_height, bool always_pad);

  explicit BoundingBox(const std::vector<PixelPosition>& pixels,
                       size_t align = 1, bool make_square = false)
      : BoundingBox(pixels, align, make_square, 0, 0, 0, false) {}

  constexpr friend bool operator==(const BoundingBox& lhs,
                                   const BoundingBox& rhs) {
    return lhs.min_ == rhs.min_ && lhs.max_ == rhs.max_;
  }

  constexpr friend bool operator!=(const BoundingBox& lhs,
                                   const BoundingBox& rhs) {
    return !(lhs == rhs);
  }

  /**
   * @return The minimum (x,y) coordinates of the bounding box. For a
   * Facet - with x-axis positive rightwards and y-axis positive
   * upward - this coordinate is the lower left point of the bounding box.
   *
   */
  const PixelPosition& Min() const { return min_; }

  /**
   * @return The maximum (x,y) coordinates of the bounding box. For a
   * Facet - with x-axis positive rightwards and y-axis positive
   * upward - this coordinate is the upper right point of the bounding box.
   */
  const PixelPosition& Max() const { return max_; }

  size_t Width() const { return static_cast<size_t>(max_.x - min_.x); }

  size_t Height() const { return static_cast<size_t>(max_.y - min_.y); }

  /**
   * Return the centre x and y of the bounding box.
   */
  PixelPosition Centre() const {
    return {(min_.x + max_.x) / 2, (min_.y + max_.y) / 2};
  }

  /**
   * Returns true if the bounding box contains the pixel, otherwise return
   * false
   */
  bool Contains(const PixelPosition& pixel) const {
    return (min_.x <= pixel.x) && (min_.y <= pixel.y) && (max_.x > pixel.x) &&
           (max_.y > pixel.y);
  }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  PixelPosition min_;  // Minimum x and y coordinates.
  PixelPosition max_;  // Maximum x and y coordinates.
};

}  // namespace schaapcommon::facets

#endif
