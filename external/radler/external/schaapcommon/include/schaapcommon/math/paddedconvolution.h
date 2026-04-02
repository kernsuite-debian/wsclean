#include "convolution.h"

#include <aocommon/image.h>

#ifndef SCHAAPCOMMON_MATH_PADDED_CONVOLUTION_H_
#define SCHAAPCOMMON_MATH_PADDED_CONVOLUTION_H_

namespace schaapcommon::math {

/**
 * Performs an in-place padded image convolution. The convolution kernel (psf)
 * is appropriately transformed so that its centre pixel is the centre of the
 * kernel. This overload requires two scratch images that have been allocated
 * with the padded dimensions.
 *
 * This operation is sometimes referred to as correlation, because the kernel is
 * not flipped as is formally done when convolving. The operation performed by
 * this function is typically used to convolve an interferometric model image
 * with the PSF.
 *
 * This overload allows specifying the scratch images, to minimize allocations
 * for use-cases that perform multiple padded convolutions. Both scratch images
 * should be allocated with the padded size on input.
 */
inline void PaddedConvolution(aocommon::Image& image,
                              const aocommon::Image& psf,
                              aocommon::Image& scratch_a,
                              aocommon::Image& scratch_b, size_t padded_width,
                              size_t padded_height) {
  assert(image.Width() == psf.Width());
  assert(image.Height() == psf.Height());
  assert(image.Size() > 0);
  assert(scratch_a.Width() * scratch_a.Height() >=
         padded_width * padded_height);
  assert(scratch_b.Width() * scratch_b.Height() >=
         padded_width * padded_height);
  assert(padded_width >= image.Width());
  assert(padded_height >= image.Height());
  using aocommon::Image;
  // scratch_a = padded psf
  Image::Untrim(scratch_a.Data(), padded_width, padded_height, psf.Data(),
                psf.Width(), psf.Height());
  // scratch_b = prepared padded psf
  PrepareConvolutionKernel(scratch_b.Data(), scratch_a.Data(), padded_width,
                           padded_height);
  // scratch_a = padded image
  Image::Untrim(scratch_a.Data(), padded_width, padded_height, image.Data(),
                image.Width(), image.Height());

  // Convolve and store in scratch_a
  Convolve(scratch_a.Data(), scratch_b.Data(), padded_width, padded_height);

  Image::Trim(image.Data(), image.Width(), image.Height(), scratch_a.Data(),
              padded_width, padded_height);
}

/**
 * Convenience overload of the above PaddedConvolution() function that does
 * not require scratch images. The required scratch images are allocated
 * inside the function. See other overload for help.
 */
inline void PaddedConvolution(aocommon::Image& image,
                              const aocommon::Image& psf, size_t padded_width,
                              size_t padded_height) {
  aocommon::Image scratch_a(padded_width, padded_height);
  aocommon::Image scratch_b(padded_width, padded_height);
  PaddedConvolution(image, psf, scratch_a, scratch_b, padded_width,
                    padded_height);
}

}  // namespace schaapcommon::math

#endif
