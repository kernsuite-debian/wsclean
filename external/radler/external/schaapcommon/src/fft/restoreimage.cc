// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "restoreimage.h"

#include <cmath>
#include <array>
#include <vector>

#include <aocommon/imagecoordinates.h>
#include <aocommon/uvector.h>

#include "convolution.h"

namespace schaapcommon {
namespace fft {
void RestoreImage(float* image_data, const float* model_data,
                  size_t image_width, size_t image_height,
                  long double beam_major_axis, long double beam_minor_axis,
                  long double beam_position_angle, long double pixel_scale_l,
                  long double pixel_scale_m, size_t thread_count) {
  if (beam_major_axis == 0.0 && beam_minor_axis == 0.0) {
    for (size_t j = 0; j != image_width * image_height; ++j) {
      image_data[j] += model_data[j];
    }
  } else {
    // Using the FWHM formula for a Gaussian:
    const long double sigma_major =
        beam_major_axis / (2.0L * std::sqrt(2.0L * std::log(2.0L)));
    const long double sigma_minor =
        beam_minor_axis / (2.0L * std::sqrt(2.0L * std::log(2.0L)));

    // Position angle is angle from North:
    const long double angle = beam_position_angle + M_PI_2;
    const long double cos_angle = std::cos(angle);
    const long double sin_angle = std::sin(angle);

    // Make rotation matrix
    std::array<long double, 4> transf;
    transf[0] = cos_angle;
    transf[1] = -sin_angle;
    transf[2] = sin_angle;
    transf[3] = cos_angle;

    const double sigma_max = std::max(std::fabs(sigma_major * transf[0]),
                                      std::fabs(sigma_major * transf[1]));
    // Multiply with scaling matrix to make variance 1.
    transf[0] /= sigma_major;
    transf[1] /= sigma_major;
    transf[2] /= sigma_minor;
    transf[3] /= sigma_minor;

    const size_t min_dimension = std::min(image_width, image_height);
    size_t bounding_box_size = std::min<size_t>(
        std::ceil(sigma_max * 40.0 / std::min(pixel_scale_l, pixel_scale_m)),
        min_dimension);
    if (bounding_box_size % 2 != 0) {
      ++bounding_box_size;
    }
    if (bounding_box_size > std::min(image_width, image_height)) {
      bounding_box_size = std::min(image_width, image_height);
    }
    aocommon::UVector<float> kernel(bounding_box_size * bounding_box_size);
    auto iter = kernel.begin();
    for (size_t y = 0; y != bounding_box_size; ++y) {
      for (size_t x = 0; x != bounding_box_size; ++x) {
        long double l;
        long double m;
        aocommon::ImageCoordinates::XYToLM<long double>(
            x, y, pixel_scale_l, pixel_scale_m, bounding_box_size,
            bounding_box_size, l, m);
        const long double l_transf = l * transf[0] + m * transf[1];
        const long double m_transf = l * transf[2] + m * transf[3];
        // Kernel value is evaluation of a Gaussian with unit-valued
        // coefficients
        const long double dist_squared =
            l_transf * l_transf + m_transf * m_transf;
        *iter = std::exp(-0.5 * dist_squared);
        ++iter;
      }
    }

    aocommon::UVector<float> convolved_model(
        model_data, model_data + image_width * image_height);

    schaapcommon::fft::ResizeAndConvolve(convolved_model.data(), image_width,
                                         image_height, kernel.data(),
                                         bounding_box_size, thread_count);
    for (size_t j = 0; j != image_width * image_height; ++j) {
      image_data[j] += convolved_model[j];
    }
  }
}
}  // namespace fft
}  // namespace schaapcommon