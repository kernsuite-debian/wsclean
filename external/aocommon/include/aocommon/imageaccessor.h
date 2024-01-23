#ifndef AOCOMMON_IMAGE_ACCESSOR_H_
#define AOCOMMON_IMAGE_ACCESSOR_H_

#include <cstddef>

namespace aocommon {

/**
 * @brief Abstract interface for loading and storing an image.
 *
 * An ImageAccessor object knows how to store a single image and load it back
 * later. When used in an interface, the caller can implement storing the image
 * in different ways (e.g. on disk, in memory, etc)."
 *
 * Load() and Store() use raw pointers to the image data instead of @ref Image
 * references. This approach allows avoiding memory reallocations.
 * An Image reallocates its data buffer when the size changes. Using raw
 * pointers allows using scratch buffers that support multiple image sizes,
 * which do not need memory reallocations when the image size changes.
 */
class ImageAccessor {
 public:
  virtual ~ImageAccessor() {}

  /**
   * @return The width of the image, in number of pixels.
   */
  virtual std::size_t Width() const = 0;

  /**
   * @return The height of the image, in number of pixels.
   */
  virtual std::size_t Height() const = 0;

  /**
   * @brief Load the image.
   *
   * @param data Location where the image will be loaded into. The location
   *             should be large enough to hold Width() * Height() pixels.
   *             The pointer is only valid during execution of this function.
   */
  virtual void Load(float* data) const = 0;

  /**
   * @brief Store the image, so it can be loaded back later.
   *
   * @param data Location of the image that must be stored. The location should
   *             contain Width() * Height() consecutive pixels.
   *             The pointer is only valid during execution of this function.
   */
  virtual void Store(const float* data) = 0;
};

}  // namespace aocommon

#endif