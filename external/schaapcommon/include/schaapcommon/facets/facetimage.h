// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FACETS_FACETIMAGE_H_
#define SCHAAPCOMMON_FACETS_FACETIMAGE_H_

#include "facet.h"

#include <aocommon/image.h>
#include <aocommon/uvector.h>

namespace schaapcommon {
namespace facets {

/**
 * This class allows copying image data between a main image and facets of that
 * main image. It has an internal data buffer for storing the facet image.
 */
class FacetImage {
 public:
  /**
   * @brief Construct a new FacetImage object.
   *
   * @param image_width Width of the main image, in pixels.
   * @param image_height Height of the main image, in pixels.
   * @param nr_spectral_terms Number of spectral terms to use.
   * @throw std::invalid_argument If image_width, image_height or
   *        nr_spectral_terms is zero.
   */
  FacetImage(size_t image_width, size_t image_height, size_t nr_spectral_terms);

  /**
   * Update the facet geometry.
   * This function initializes the FacetImage for using the given facet.
   * Call this function before calling any other functions.
   *
   * @param facet The facet geometry. This facet should have a properly padded
   *        and aligned bounding box. This bounding box may extend beyond the
   *        borders of the main image due to padding or alignment.
   *        All horizontal intersections should be within the mainimage, though.
   * @param trimmed True: Use trimmed facet bounding box. False: Use untrimmed
   *        facet bounding box.
   */
  void SetFacet(const Facet& facet, bool trimmed);

  /**
   * @brief Get the Facet object previously set in SetFacet()
   */
  const std::optional<Facet>& GetFacet() const { return facet_; }

  size_t Width() const { return box_.Width(); }
  size_t Height() const { return box_.Height(); }

  size_t FullImageWidth() const { return image_width_; }
  size_t FullImageHeight() const { return image_height_; }

  /**
   * @return The number of pixels that this facet is displaced in X direction,
   * relative to origin of main image.
   */
  int OffsetX() const { return box_.Min().x; }

  /**
   * @return The number of pixels that this facet is displaced in Y direction,
   * relative to origin of main image.
   */
  int OffsetY() const { return box_.Min().y; }

  const float* Data(size_t spectral_term) const {
    return data_[spectral_term].data();
  }

  float* Data(size_t spectral_term) { return data_[spectral_term].data(); }

  /**
   * @brief Multiply facet image data for all spectral terms with a scalar.
   */
  FacetImage& operator*=(float factor);

  /**
   * Fill the facet buffer with image data.
   *
   * @param images Vector of images.
   *        The outer vector length should match the number of spectral terms.
   *        The inner vector size should should match the image width and
   *        image height values passed to the constructor.
   * @throw std::invalid_argument if 'images' has an incorrect size.
   * @throw std::runtime_error if SetFacet was not called.
   */
  void CopyToFacet(const std::vector<aocommon::UVector<float>>& images);

  /**
   * Fill the facet buffer with the image data given as a buffer.
   *
   * @param images Vector of buffered images.
   *        The vector length should match the number of spectral terms.
   *        The buffer is assumed to match the image width and image height
   *        values passed to the constructor
   * @throw std::runtime_error if SetFacet was not called.
   */
  void CopyToFacet(const std::vector<const float*>& images);

  /**
   * Add data from the facet buffer to main images.
   * @param images Vector of images.
   *        The vector length should match the number of spectral terms.
   *        The width and height of each image should match the values passed
   *        to the constructor.
   */
  void AddToImage(const std::vector<float*>& images) const;

  /**
   * Returns an image that contains a value of 1 inside the facet
   * and a value of 0 outside the facet. The image will have the
   * size of the bounding box of the facet, and can be used in a call
   * to AddWithMask().
   */
  aocommon::Image MakeMask() const;

  /**
   * Add this facet to the provided full image and weight image,
   * using a specified weight mask. In case the facet extends beyond
   * the right (high x) or bottom (high y) edges, the facet image is
   * trimmed on this side to fall inside the full image. This
   * situation might arise because of the alignment requirements.
   * @param[in,out] image  Full image to which the data from this facet will be
   * added.
   * @param[in,out] weight Weight image (same size as the full image) that will
   * be updated.
   * @param[in] mask Image with values >= 0 that are used as weights for the
   * facet. This can typically be created with @ref MakeMask(), after
   * modifications (e.g. feathering).
   * @param Spectral index term.
   */
  void AddWithMask(float* image, float* weight, const aocommon::Image& mask,
                   size_t term) const;

  /**
   * @brief Multiply the image values enclosed by the facet with scalar valued
   * @param factor
   *
   * @param images Vector of images.
   *         The vector length should match the number of spectral terms.
   *        The width and height of each image should match the values passed
   *        to the constructor.
   * @param factor Scalar used for multiplying the image data inside the facet.
   */
  void MultiplyImageInsideFacet(std::vector<float*>& images,
                                float factor) const;

 private:
  /**
   * @brief Pre-calculate the horizontal intersections (xmin, xmax)
   * of for each y-coordinate row in the facet bounding box with the
   * main image.
   */
  void SetHorizontalIntersections(const Facet& facet);

  /**
   * A vector with one element for each spectral frequency term.
   * Each inner vector holds the data for a facet image.
   */
  std::vector<std::vector<float>> data_;
  std::optional<Facet> facet_;
  const size_t image_width_;   ///< Width of the main image.
  const size_t image_height_;  ///< Height of the main image.
  BoundingBox box_;  ///< Facet bounding box, including padding/alignment.

  /**
   * For each y value in box_, a vector with start and end x coordinate in the
   * main image.
   */
  std::vector<std::vector<std::pair<int, int>>> horizontal_intersections_;
};
}  // namespace facets
}  // namespace schaapcommon
#endif
