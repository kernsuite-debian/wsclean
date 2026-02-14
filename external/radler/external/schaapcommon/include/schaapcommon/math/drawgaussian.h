#ifndef SCHAAPCOMMON_MATH_DRAW_GAUSSIAN_H_
#define SCHAAPCOMMON_MATH_DRAW_GAUSSIAN_H_

#include <cstring>

#include "ellipse.h"

namespace schaapcommon::math {

/**
 * Draw a two-dimensional Gaussian onto a gridded image. Units are given in
 * pixels.
 */
void DrawGaussianToXy(float* image_data, size_t image_width,
                      size_t image_height, double source_x, double source_y,
                      const Ellipse& ellipse, long double integrated_flux);

/**
 * Draw a two-dimensional Gaussian onto an lm-space image.
 * As of yet, @p pixel_scale_l and @p pixel_scale_m must be equal.
 */
void DrawGaussianToLm(float* image_data, size_t image_width,
                      size_t image_height, long double phase_centre_ra,
                      long double phase_centre_dec, long double pixel_scale_l,
                      long double pixel_scale_m, long double l_shift,
                      long double m_shift, long double source_ra,
                      long double source_dec, const Ellipse& ellipse,
                      long double integrated_flux);

}  // namespace schaapcommon::math

#endif
