// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "boundingbox.h"

#include <cassert>

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

namespace {
constexpr int ApplyExtraSpaceToMin(int min_value, int extra_space,
                                   bool always_pad) {
  if (always_pad) {
    return min_value - extra_space;
  } else {
    // Extra space for feathering is not needed beyond the full image, so if the
    // min_value would become negative, it is set to zero. If the min_value was
    // already negative, it is left unchanged.
    return min_value > 0 ? std::max(0, min_value - extra_space) : min_value;
  }
}

constexpr int ApplyExtraSpaceToMax(int max_value, int extra_space,
                                   int full_size, bool always_pad) {
  if (always_pad) {
    return max_value + extra_space;
  } else {
    return max_value < full_size ? std::min(full_size, max_value + extra_space)
                                 : max_value;
  }
}

}  // namespace

namespace schaapcommon {
namespace facets {

void PixelPosition::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(x).UInt32(y);
}

void PixelPosition::Unserialize(aocommon::SerialIStream& stream) {
  stream.UInt32(x).UInt32(y);
}

BoundingBox::BoundingBox(const std::vector<PixelPosition>& pixels, size_t align,
                         bool make_square, size_t extra_space,
                         size_t full_width, size_t full_height,
                         bool always_pad) {
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

  if (extra_space) {
    assert(full_width && full_height);
    min_.x = ApplyExtraSpaceToMin(min_.x, extra_space, always_pad);
    min_.y = ApplyExtraSpaceToMin(min_.y, extra_space, always_pad);
    max_.x = ApplyExtraSpaceToMax(max_.x, extra_space, full_width, always_pad);
    max_.y = ApplyExtraSpaceToMax(max_.y, extra_space, full_height, always_pad);
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

}  // namespace facets
}  // namespace schaapcommon
