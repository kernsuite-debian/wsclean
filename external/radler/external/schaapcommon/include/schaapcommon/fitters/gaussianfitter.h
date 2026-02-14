// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_GAUSSIAN_FITTER_H_
#define SCHAAPCOMMON_FITTERS_GAUSSIAN_FITTER_H_

#include <cstring>

#include "../math/ellipse.h"

namespace schaapcommon {
namespace fitters {

/**
 * Deconvolve two Gaussian kernels. The result convolved with @p smallest
 * results again in @p largest. If @p smallest is larger than @p largest,
 * an ellipse is returned with all parameters set to NaN.
 *
 * As the output axes are returned with the same units as the input axes,
 * the input may be both in FWHM or both in sigma.
 */
math::Ellipse DeconvolveGaussian(const math::Ellipse& largest,
                                 const math::Ellipse& smallest);

/**
 * Fit Gaussian parameters to the shape in the center of the image.
 * @param beam_est provides the initial value for the beam's size. It must be
 * somewhat close and must be non-zero for the algorithm to converge (quickly).
 */
math::Ellipse Fit2DGaussianCentred(const float* image, size_t width,
                                   size_t height, double beam_est,
                                   double box_scale_factor = 10.0,
                                   bool verbose = false);

/**
 * Fits a circular (one parameter) Gaussian to the shape in the centre of the
 * image.
 * @param [in,out] beam_size on input, the initial value for beam size. On
 * output, the fit result.
 */
void Fit2DCircularGaussianCentred(const float* image, size_t width,
                                  size_t height, double& beam_size,
                                  double box_scale_factor = 10.0);

/**
 * Fit all Gaussian parameters to a shape in an image. Parameters
 * @p pos_x, @p pos_y and @p beam_major are input and output function
 * parameters: on input they provide initial values for the fitting,
 * and on output they provide the result. Parameters @p val,
 * @p beam_minor, @p beam_pa and @p floor_level are only used for outputting
 * the result.
 * @param [out] val amplitude of Gaussian.
 * @param [in,out] floor_level If not nullptr, will be used to store the floor
 * level parameter. If nullptr, this value is not fitted and assumed to be zero.
 */
void Fit2DGaussianFull(const float* image, size_t width, size_t height,
                       double& val, double& pos_x, double& pos_y,
                       double& beam_major, double& beam_minor, double& beam_pa,
                       double* floor_level = nullptr);

}  // namespace fitters
}  // namespace schaapcommon
#endif
