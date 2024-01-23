// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "gaussianfitter.h"

#include <cmath>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

using aocommon::Matrix2x2;

namespace schaapcommon {
namespace fitters {

namespace {

const long double kSigmaToBeam = 2.0L * std::sqrt(2.0L * std::log(2.0L));

void ToAnglesAndFwhm(double sx, double sy, double beta, double& ellipse_major,
                     double& ellipse_minor, double& ellipse_phase_angle) {
  const double beta_factor = 1.0 - beta * beta;
  double cov[4];
  cov[0] = sx * sx / beta_factor;
  cov[1] = beta * sx * sy / beta_factor;
  cov[2] = cov[1];
  cov[3] = sy * sy / beta_factor;

  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
  if (std::isfinite(e1)) {
    ellipse_major = std::sqrt(std::fabs(e1)) * kSigmaToBeam;
    ellipse_minor = std::sqrt(std::fabs(e2)) * kSigmaToBeam;
    if (ellipse_major < ellipse_minor) {
      std::swap(ellipse_major, ellipse_minor);
      vec1[0] = vec2[0];
      vec1[1] = vec2[1];
    }
    ellipse_phase_angle = -std::atan2(vec1[0], vec1[1]);
  } else {
    ellipse_major = std::sqrt(std::fabs(sx)) * kSigmaToBeam;
    ellipse_minor = std::sqrt(std::fabs(sx)) * kSigmaToBeam;
    ellipse_phase_angle = 0.0;
  }
}

/**
 * Calculates a two dimensional gaussian with the specified parameters.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param sx Sigma value for the x direction.
 * @param sy Sigma value for the y direction.
 * @param beta Beta value.
 */
double GaussCentered(double x, double y, double sx, double sy, double beta) {
  return std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                  y * y / (2.0 * sy * sy));
}

/**
 * Calculates a circular two dimensional gaussian with the specified parameters.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param s Sigma value for both x and y directions.
 */
double GaussCircularCentered(double x, double y, double s) {
  return std::exp((-x * x - y * y) / (2.0 * s * s));
}

/**
 * Fitting function for SingleFit2DGaussianCentred(). Calculates the sum of the
 * squared errors(/residuals).
 */
int FittingCentered(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double e = GaussCentered(x, y, sx, sy, beta) - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingCircularCentered(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double e = GaussCircularCentered(x, y, s) - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Derivative function belong with SingleFit2DGaussianCentred().
 */
int FittingDerivativeCentered(const gsl_vector* xvec, void* data,
                              gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      gsl_matrix_set(J, data_index, 0, dsx);
      gsl_matrix_set(J, data_index, 1, dsy);
      gsl_matrix_set(J, data_index, 2, dbeta);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeCircularCentered(const gsl_vector* xvec, void* data,
                                      gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2, y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double exp_term = GaussCircularCentered(x, y, s);
      // derivative of exp((-x*x - y*y)/(2.0*s*s)) to s
      // = (-x*x - y*y)/2.0*-2/(s*s*s)
      // = (-x*x - y*y)/(-s*s*s)
      // = (x*x + y*y)/(s*s*s)
      double ds = ((x * x + y * y) / (s * s * s)) * exp_term;
      gsl_matrix_set(J, data_index, 0, ds);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Squared error and derivative function together.
 */
int FittingBothCentered(const gsl_vector* x, void* data, gsl_vector* f,
                        gsl_matrix* J) {
  FittingCentered(x, data, f);
  FittingDerivativeCentered(x, data, J);
  return GSL_SUCCESS;
}

int FittingBothCircularCentered(const gsl_vector* x, void* data, gsl_vector* f,
                                gsl_matrix* J) {
  FittingCircularCentered(x, data, f);
  FittingDerivativeCircularCentered(x, data, J);
  return GSL_SUCCESS;
}

int FittingWithAmplitude(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - x_mid) * scale;
      double e =
          GaussCentered(xS, yS, sx, sy, beta) * v - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeWithAmplitude(const gsl_vector* xvec, void* data,
                                   gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale)) {
    std::cout << "GSL_EDOM\n";
    return GSL_EDOM;
  }
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      // TODO I need to go over the signs -- ds, dy, dsx, dsy in particular
      double x = xc + (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dv = exp_term;
      exp_term *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * exp_term;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * exp_term;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      gsl_matrix_set(J, data_index, 0, dv);
      gsl_matrix_set(J, data_index, 1, dx);
      gsl_matrix_set(J, data_index, 2, dy);
      gsl_matrix_set(J, data_index, 3, dsx);
      gsl_matrix_set(J, data_index, 4, dsy);
      gsl_matrix_set(J, data_index, 5, dbeta);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingBothWithAmplitude(const gsl_vector* x, void* data, gsl_vector* f,
                             gsl_matrix* J) {
  FittingWithAmplitude(x, data, f);
  FittingDerivativeWithAmplitude(x, data, J);
  return GSL_SUCCESS;
}

int FittingWithAmplitudeAndFloor(const gsl_vector* xvec, void* data,
                                 gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const double fl = gsl_vector_get(xvec, 6);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale)) {
    return GSL_EDOM;
  }
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - x_mid) * scale;
      double e = GaussCentered(xS, yS, sx, sy, beta) * v -
                 fitter.Image()[data_index] + fl;
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeWithAmplitudeAndFloor(const gsl_vector* xvec, void* data,
                                           gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double x = xc + (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dv = exp_term;
      exp_term *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * exp_term;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * exp_term;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      double dfl = 1.0;
      gsl_matrix_set(J, data_index, 0, dv);
      gsl_matrix_set(J, data_index, 1, dx);
      gsl_matrix_set(J, data_index, 2, dy);
      gsl_matrix_set(J, data_index, 3, dsx);
      gsl_matrix_set(J, data_index, 4, dsy);
      gsl_matrix_set(J, data_index, 5, dbeta);
      gsl_matrix_set(J, data_index, 6, dfl);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingBothWithAmplitudeAndFloor(const gsl_vector* x, void* data,
                                     gsl_vector* f, gsl_matrix* J) {
  FittingWithAmplitudeAndFloor(x, data, f);
  FittingDerivativeWithAmplitudeAndFloor(x, data, J);
  return GSL_SUCCESS;
}

}  // namespace

void GaussianFitter::Fit2DGaussianCentred(
    const float* image, size_t width, size_t height, double beam_estimate,
    double& beam_major, double& beam_minor, double& beam_phase_angle,
    double box_scale_factor, bool verbose) {
  size_t preferred_size = std::max<size_t>(
      std::ceil(box_scale_factor), std::ceil(beam_estimate * box_scale_factor));
  if (preferred_size % 2 != 0) ++preferred_size;
  if (preferred_size < width || preferred_size < height) {
    size_t n_iterations = 0;
    bool box_was_large_enough;
    do {
      size_t box_width = std::min(preferred_size, width);
      size_t box_height = std::min(preferred_size, height);
      if (verbose) std::cout << "Fit initial value:" << beam_estimate << "\n";
      Fit2DGaussianCentredInBox(image, width, height, beam_estimate, beam_major,
                                beam_minor, beam_phase_angle, box_width,
                                box_height, verbose);
      if (verbose) {
        std::cout << "Fit result:" << beam_major << " x " << beam_minor
                  << " px, " << beam_phase_angle << " (box was " << box_width
                  << " x " << box_height << ")\n";
      }

      box_was_large_enough =
          (beam_major * box_scale_factor * 0.8 < box_width ||
           box_width >= width) &&
          (beam_major * box_scale_factor * 0.8 < box_height ||
           box_height >= height);
      if (!box_was_large_enough) {
        preferred_size =
            std::max<size_t>(std::ceil(box_scale_factor),
                             std::ceil(beam_major * box_scale_factor));
        if (preferred_size % 2 != 0) ++preferred_size;
        beam_estimate = std::max(beam_major, beam_estimate);
      }
      ++n_iterations;
    } while (!box_was_large_enough && n_iterations < 5);
  } else {
    if (verbose) std::cout << "Image is as large as the fitting box.\n";
    SingleFit2DGaussianCentred(image, width, height, beam_estimate, beam_major,
                               beam_minor, beam_phase_angle, verbose);
  }
}

void GaussianFitter::Fit2DCircularGaussianCentred(const float* image,
                                                  size_t width, size_t height,
                                                  double& beam_size,
                                                  double box_scale_factor) {
  double initial_value = beam_size;
  size_t preferred_size = std::max<size_t>(
      std::ceil(box_scale_factor), std::ceil(beam_size * box_scale_factor));
  if (preferred_size % 2 != 0) ++preferred_size;
  if (preferred_size < width || preferred_size < height) {
    size_t box_width = std::min(preferred_size, width);
    size_t box_height = std::min(preferred_size, height);
    size_t n_iterations = 0;
    bool box_was_large_enough;
    do {
      Fit2DCircularGaussianCentredInBox(image, width, height, beam_size,
                                        box_width, box_height);

      box_was_large_enough = (beam_size * box_scale_factor * 0.8 < box_width ||
                              width >= box_width) &&
                             (beam_size * box_scale_factor * 0.8 < box_height ||
                              height >= box_height);
      if (!box_was_large_enough) {
        preferred_size =
            std::max<size_t>(std::ceil(box_scale_factor),
                             std::ceil(beam_size * box_scale_factor));
        if (preferred_size % 2 != 0) ++preferred_size;
        beam_size = std::max(initial_value, beam_size);
      }
      ++n_iterations;
    } while (!box_was_large_enough && n_iterations < 5);
  } else {
    SingleFit2DCircularGaussianCentred(image, width, height, beam_size);
  }
}

void GaussianFitter::Fit2DGaussianFull(const float* image, size_t width,
                                       size_t height, double& val,
                                       double& pos_x, double& pos_y,
                                       double& beam_major, double& beam_minor,
                                       double& beam_phase_angle,
                                       double* floor_level) {
  size_t preferred_size = std::max<size_t>(10, std::ceil(beam_major * 10.0));
  if (preferred_size % 2 != 0) ++preferred_size;
  if (preferred_size < width || preferred_size < height) {
    size_t x_start =
        std::max<int>(0, int(std::round(pos_x)) - int(preferred_size) / 2);
    size_t x_end =
        std::min(width, size_t(std::round(pos_x)) + preferred_size / 2);
    size_t y_start =
        std::max<int>(0, int(std::round(pos_y)) - int(preferred_size) / 2);
    size_t y_end =
        std::min(height, size_t(std::round(pos_y)) + preferred_size / 2);
    size_t n_iterations = 0;
    bool box_was_large_enough;
    do {
      Fit2DGaussianWithAmplitudeInBox(
          image, width, height, val, pos_x, pos_y, beam_major, beam_minor,
          beam_phase_angle, floor_level, x_start, x_end, y_start, y_end);

      size_t box_width = x_end - x_start;
      size_t box_height = y_end - y_start;
      box_was_large_enough =
          (beam_major * 4.0 < box_width || width >= box_width) &&
          (beam_major * 4.0 < box_height || height >= box_height);
      if (!box_was_large_enough) {
        preferred_size = std::max<size_t>(10, std::ceil(beam_major * 10.0));
        if (preferred_size % 2 != 0) ++preferred_size;
      }
      ++n_iterations;
    } while (!box_was_large_enough && n_iterations < 5);
  } else {
    Fit2DGaussianWithAmplitude(image, width, height, val, pos_x, pos_y,
                               beam_major, beam_minor, beam_phase_angle,
                               floor_level);
  }
}

void GaussianFitter::Fit2DGaussianCentredInBox(
    const float* image, size_t width, size_t height, double beam_estimate,
    double& beam_major, double& beam_minor, double& beam_phase_angle,
    size_t box_width, size_t box_height, bool verbose) {
  const size_t x_start = (width - box_width) / 2;
  const size_t y_start = (height - box_height) / 2;
  aocommon::UVector<float> small_image(box_width * box_height);
  for (size_t y = y_start; y != (height + box_height) / 2; ++y) {
    std::copy_n(&image[y * width + x_start], box_width,
                &small_image[(y - y_start) * box_width]);
  }

  SingleFit2DGaussianCentred(small_image.data(), box_width, box_height,
                             beam_estimate, beam_major, beam_minor,
                             beam_phase_angle, verbose);
}

void GaussianFitter::Fit2DCircularGaussianCentredInBox(
    const float* image, size_t width, size_t height, double& beam_size,
    size_t box_width, size_t box_height) {
  const size_t x_start = (width - box_width) / 2;
  const size_t y_start = (height - box_height) / 2;
  aocommon::UVector<float> small_image(box_width * box_height);
  for (size_t y = y_start; y != (height + box_height) / 2; ++y) {
    std::copy_n(&image[y * width + x_start], box_width,
                &small_image[(y - y_start) * box_width]);
  }

  SingleFit2DCircularGaussianCentred(small_image.data(), box_width, box_height,
                                     beam_size);
}

void GaussianFitter::SingleFit2DGaussianCentred(
    const float* image, size_t width, size_t height, double beam_estimate,
    double& beam_major, double& beam_minor, double& beam_phase_angle,
    bool verbose) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 3);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingCentered;
  fdf.df = &FittingDerivativeCentered;
  fdf.fdf = &FittingBothCentered;
  fdf.n = width_ * height_;
  fdf.p = 3;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  double initial_values_array[3] = {
      beam_estimate / (scale_factor_ * double(kSigmaToBeam)),
      beam_estimate / (scale_factor_ * double(kSigmaToBeam)), 0.0};
  gsl_vector_view initial_values =
      gsl_vector_view_array(initial_values_array, 3);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initial_values.vector);

  int status;
  size_t iter = 0;
  do {
    if (verbose) std::cout << "Iteration " << iter << ": ";
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-7, 1.0e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  double sx = gsl_vector_get(solver->x, 0), sy = gsl_vector_get(solver->x, 1),
         beta = gsl_vector_get(solver->x, 2);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beam_major, beam_minor, beam_phase_angle);
  beam_major *= scale_factor_;
  beam_minor *= scale_factor_;
}

void GaussianFitter::SingleFit2DCircularGaussianCentred(const float* image,
                                                        size_t width,
                                                        size_t height,
                                                        double& beam_size) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 1);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingCircularCentered;
  fdf.df = &FittingDerivativeCircularCentered;
  fdf.fdf = &FittingBothCircularCentered;
  fdf.n = width_ * height_;
  fdf.p = 1;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  double initial_values_array[1] = {beam_size /
                                    (scale_factor_ * double(kSigmaToBeam))};
  gsl_vector_view initial_values =
      gsl_vector_view_array(initial_values_array, 1);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initial_values.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-7, 1.0e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  const double s = gsl_vector_get(solver->x, 0);
  gsl_multifit_fdfsolver_free(solver);

  beam_size = s * kSigmaToBeam * scale_factor_;
}

void GaussianFitter::Fit2DGaussianWithAmplitudeInBox(
    const float* image, size_t width, size_t /*height*/, double& val,
    double& pos_x, double& pos_y, double& beam_major, double& beam_minor,
    double& beam_phase_angle, double* floor_level, size_t x_start, size_t x_end,
    size_t y_start, size_t y_end) {
  const size_t box_width = x_end - x_start;
  const size_t box_height = y_end - y_start;
  aocommon::UVector<float> small_image(box_width * box_height);
  for (size_t y = y_start; y != y_end; ++y) {
    std::copy_n(&image[y * width + x_start], box_width,
                &small_image[(y - y_start) * box_width]);
  }

  pos_x -= x_start;
  pos_y -= y_start;
  Fit2DGaussianWithAmplitude(small_image.data(), box_width, box_height, val,
                             pos_x, pos_y, beam_major, beam_minor,
                             beam_phase_angle, floor_level);
  pos_x += x_start;
  pos_y += y_start;
}

/**
 * Fits the position, size and amplitude of a Gaussian. If floor_level is not
 * a nullptr, the floor (background level, or zero level) is fitted too.
 */
void GaussianFitter::Fit2DGaussianWithAmplitude(
    const float* image, size_t width, size_t height, double& val, double& pos_x,
    double& pos_y, double& beam_major, double& beam_minor,
    double& beam_phase_angle, double* floor_level) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  if (floor_level == nullptr) {
    Fit2DGaussianWithAmplitude(val, pos_x, pos_y, beam_major, beam_minor,
                               beam_phase_angle);
  } else {
    Fit2DGaussianWithAmplitudeWithFloor(val, pos_x, pos_y, beam_major,
                                        beam_minor, beam_phase_angle,
                                        *floor_level);
  }
}

void GaussianFitter::Fit2DGaussianWithAmplitude(double& val, double& pos_x,
                                                double& pos_y,
                                                double& beam_major,
                                                double& beam_minor,
                                                double& beam_phase_angle) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 6);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingWithAmplitude;
  fdf.df = &FittingDerivativeWithAmplitude;
  fdf.fdf = &FittingBothWithAmplitude;
  fdf.n = width_ * height_;
  fdf.p = 6;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  x_init_ = -(pos_x - width_ / 2) / scale_factor_;
  y_init_ = -(pos_y - height_ / 2) / scale_factor_;
  double initial_values_array[6] = {
      val,
      x_init_,
      y_init_,
      beam_major / (scale_factor_ * double(kSigmaToBeam)),
      beam_major / (scale_factor_ * double(kSigmaToBeam)),
      0.0};
  gsl_vector_view initial_values =
      gsl_vector_view_array(initial_values_array, 6);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initial_values.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-7, 1.0e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  pos_x = -1.0 * gsl_vector_get(solver->x, 1) * scale_factor_ + width_ / 2;
  pos_y = -1.0 * gsl_vector_get(solver->x, 2) * scale_factor_ + height_ / 2;
  double sx = gsl_vector_get(solver->x, 3), sy = gsl_vector_get(solver->x, 4),
         beta = gsl_vector_get(solver->x, 5);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beam_major, beam_minor, beam_phase_angle);
  beam_major *= scale_factor_;
  beam_minor *= scale_factor_;
}

void GaussianFitter::Fit2DGaussianWithAmplitudeWithFloor(
    double& val, double& pos_x, double& pos_y, double& beam_major,
    double& beam_minor, double& beam_phase_angle, double& floor_level) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 7);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingWithAmplitudeAndFloor;
  fdf.df = &FittingDerivativeWithAmplitudeAndFloor;
  fdf.fdf = &FittingBothWithAmplitudeAndFloor;
  fdf.n = width_ * height_;
  fdf.p = 7;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  x_init_ = -(pos_x - width_ / 2) / scale_factor_;
  y_init_ = -(pos_y - height_ / 2) / scale_factor_;
  double initial_values_array[7] = {
      val,
      x_init_,
      y_init_,
      beam_major / (scale_factor_ * double(kSigmaToBeam)),
      beam_major / (scale_factor_ * double(kSigmaToBeam)),
      0.0,
      0.0};
  gsl_vector_view initial_values =
      gsl_vector_view_array(initial_values_array, 7);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initial_values.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-7, 1.0e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  pos_x = -1.0 * gsl_vector_get(solver->x, 1) * scale_factor_ + width_ / 2;
  pos_y = -1.0 * gsl_vector_get(solver->x, 2) * scale_factor_ + height_ / 2;
  double sx = gsl_vector_get(solver->x, 3);
  double sy = gsl_vector_get(solver->x, 4);
  double beta = gsl_vector_get(solver->x, 5);
  floor_level = gsl_vector_get(solver->x, 6);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beam_major, beam_minor, beam_phase_angle);
  beam_major *= scale_factor_;
  beam_minor *= scale_factor_;
}
}  // namespace fitters
}  // namespace schaapcommon