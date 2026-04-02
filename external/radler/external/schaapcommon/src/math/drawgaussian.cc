#include "drawgaussian.h"

#include <cassert>
#include <cmath>

#include <aocommon/imagecoordinates.h>

using aocommon::ImageCoordinates;

namespace {
inline long double Gaussian(long double x) {
  // Evaluate unnormalized Gaussian (unit valued peak, but not a unit integral).
  return std::exp(-0.5L * x * x);
}
}  // namespace

namespace schaapcommon::math {

void DrawGaussianToXy(float* image_data, size_t image_width,
                      size_t image_height, double source_x, double source_y,
                      const Ellipse& ellipse, long double integrated_flux) {
  // Using the FWHM formula for a Gaussian:
  const long double fwhm_constant = 2.0L * std::sqrt(2.0L * M_LN2);
  const long double sigma_major_axis = ellipse.major / fwhm_constant;
  const long double sigma_minor_axis = ellipse.minor / fwhm_constant;

  // Position angle is angle from North:
  const long double angle = ellipse.position_angle + M_PI_2;
  const long double cos_angle = std::cos(angle);
  const long double sin_angle = std::sin(angle);

  // Make rotation matrix
  long double transf[4];
  transf[0] = cos_angle;
  transf[1] = -sin_angle;
  transf[2] = sin_angle;
  transf[3] = cos_angle;

  const double sigma_max = std::max(std::fabs(sigma_major_axis * transf[0]),
                                    std::fabs(sigma_major_axis * transf[1]));
  // Multiply with scaling matrix to make variance 1.
  transf[0] /= sigma_major_axis;
  transf[1] /= sigma_major_axis;
  transf[2] /= sigma_minor_axis;
  transf[3] /= sigma_minor_axis;

  // Calculate the bounding box
  const int bounding_box_size = std::ceil(sigma_max * 20.0);
  const int x_left =
      std::clamp<int>(source_x - bounding_box_size, 0, image_width);
  const int x_right =
      std::clamp<int>(source_x + bounding_box_size, x_left, image_width);
  const int y_top =
      std::clamp<int>(source_y - bounding_box_size, 0, image_height);
  const int y_bottom =
      std::clamp<int>(source_y + bounding_box_size, y_top, image_height);

  std::vector<double> values;
  values.reserve((x_right - x_left) * (y_bottom - y_top));
  double flux_sum = 0.0;
  for (int y = y_top; y != y_bottom; ++y) {
    for (int x = x_left; x != x_right; ++x) {
      const long double x_transf =
          (x - source_x) * transf[0] + (y - source_y) * transf[1];
      const long double y_transf =
          (x - source_x) * transf[2] + (y - source_y) * transf[3];
      const long double dist =
          std::sqrt(x_transf * x_transf + y_transf * y_transf);
      long double v = Gaussian(dist);
      flux_sum += static_cast<double>(v);
      values.emplace_back(v);
    }
  }
  const double* iter = values.data();
  // While the integral of a continuous Gaussian is known, the Gaussian that is
  // drawn here is a sampled (gridded) Gaussian. Therefore, it is not accurate
  // to use the standard Gaussian formula to normalize the Gaussian.
  // TODO if the Gaussian cuts the side of the image, this leads to unexpected
  // results.
  const double factor = integrated_flux / flux_sum;
  for (int y = y_top; y != y_bottom; ++y) {
    float* image_data_ptr = image_data + y * image_width + x_left;
    for (int x = x_left; x != x_right; ++x) {
      *image_data_ptr += *iter * factor;
      ++image_data_ptr;
      ++iter;
    }
  }
}

void DrawGaussianToLm(float* image_data, size_t image_width,
                      size_t image_height, long double phase_centre_ra,
                      long double phase_centre_dec, long double pixel_scale_l,
                      long double pixel_scale_m, long double l_shift,
                      long double m_shift, long double source_ra,
                      long double source_dec, const Ellipse& ellipse,
                      long double integrated_flux) {
  assert(pixel_scale_l == pixel_scale_m);
  long double source_l;
  long double source_m;
  ImageCoordinates::RaDecToLM(source_ra, source_dec, phase_centre_ra,
                              phase_centre_dec, source_l, source_m);
  long double source_x;
  long double source_y;
  ImageCoordinates::LMToXYfloat<long double>(
      source_l - l_shift, source_m - m_shift, pixel_scale_l, pixel_scale_m,
      image_width, image_height, source_x, source_y);

  Ellipse xy_ellipse;
  xy_ellipse.major = ellipse.major / pixel_scale_l;
  xy_ellipse.minor = ellipse.minor / pixel_scale_l;

  // Position angle is North through East. Because l increases to the left
  // whereas x increases to the right (see e.g. ImageCoordinates::LMToXY()),
  // the sign is flipped here.
  xy_ellipse.position_angle = -ellipse.position_angle;
  DrawGaussianToXy(image_data, image_width, image_height, source_x, source_y,
                   xy_ellipse, integrated_flux);
}

}  // namespace schaapcommon::math
