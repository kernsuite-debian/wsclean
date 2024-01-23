// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FFT_RESTORE_IMAGE_H_
#define SCHAAPCOMMON_FFT_RESTORE_IMAGE_H_

#include <cstddef>

namespace schaapcommon {
namespace fft {

/**
 * @brief Restore a diffuse image (e.g. produced with a multi-scale clean) using
 * FFT convolution.
 *
 * Caller is assumed to ensure thread safety of the FFT plan creation. This can
 * be achieved by running \c schaapcommon::fft::MakeFftwPlannerThreadSafe() at a
 * sufficiently high level in the calling code.
 *
 * @param image_data (Residual??) image data buffer
 * @param model_data Model image data buffer
 * @param image_width Image height in pixels
 * @param image_height Image height in pixels
 * @param beam_major_axis Length of beam major axis [rad??]
 * @param beam_minor_axis Length of beam minor axis [rad??]
 * @param beam_position_angle Beam position angle [rad]
 * @param pixel_scale_l Pixel scale in L-direction [rad]
 * @param pixel_scale_m Pixel scale in M-direction [rad]
 * @param thread_count Number of threads to use in convolution.
 */
void RestoreImage(float* image_data, const float* model_data,
                  size_t image_width, size_t image_height,
                  long double beam_major_axis, long double beam_minor_axis,
                  long double beam_position_angle, long double pixel_scale_l,
                  long double pixel_scale_m, size_t thread_count);

}  // namespace fft
}  // namespace schaapcommon
#endif