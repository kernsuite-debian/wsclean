// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FFT_CONVOLUTION_H_
#define SCHAAPCOMMON_FFT_CONVOLUTION_H_

#include <cstring>

namespace schaapcommon {
namespace fft {

/**
 * @brief Make the FFTW float planner thread safe.
 */
void MakeFftwfPlannerThreadSafe();

/**
 * Convolve an image with a smaller kernel. No preparation of either image is
 * needed.
 *
 * This function assumes that (n+1)/2 is the middle pixel for uneven image
 * sizes. In the case of even image sizes, the middle falls between two
 * pixels.
 *
 * To ensure thread safety make sure that \c MakeFftwPlannerThreadSafe() is
 * called at a sufficiently high level before calling this function.
 */
void ResizeAndConvolve(float* image, size_t image_width, size_t image_height,
                       const float* kernel, size_t kernel_size,
                       size_t thread_count);

/**
 * Prepare a smaller kernel for convolution with \c Convolve(). The kernel
 * is zero-padded and translated such that it is correctly centered for the
 * convolution with the larger image. When the kernel is used more often, it
 * is more efficient to call \c PrepareConvolutionKernel() or
 * \c PrepareSmallConvolutionKernel() once and multiple times \c Convolve(),
 * than calling \c Convolve() or \c ResizeAndConvolve() multiple times.
 */
void PrepareSmallConvolutionKernel(float* dest, size_t image_width,
                                   size_t image_height, const float* kernel,
                                   size_t kernel_size, size_t thread_count);
/**
 * Prepare a kernel for convolution with \c ConvolveSameSize(), by translating
 * the input buffer such that the center value ends-up at postion 0. The
 * kernel should be of the same size as the image to be convolved.
 * Otherwise, \c PrepareSmallConvolutionKernel() should be used.
 */
void PrepareConvolutionKernel(float* dest, const float* source,
                              size_t image_width, size_t image_height,
                              size_t thread_count);

/**
 * Convolve an image with an already prepared kernel of the same size.
 *
 * To ensure thread safety make sure that \c MakeFftwPlannerThreadSafe() is
 * called at a sufficiently high level before calling this function.
 */
void Convolve(float* image, const float* kernel, size_t image_width,
              size_t image_height, size_t thread_count);
}  // namespace fft
}  // namespace schaapcommon

#endif
