// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "nlplfitter.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <aocommon/uvector.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

namespace schaapcommon {
namespace fitters {

class NLPLFitterData {
 public:
  using PointVec = aocommon::UVector<std::pair<double, double>>;
  PointVec points;
  size_t n_terms;
  gsl_multifit_fdfsolver* solver;

  static int Fitting(const gsl_vector* x_vector, void* data, gsl_vector* f) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);
    const double exponent = gsl_vector_get(x_vector, 0);
    const double factor = gsl_vector_get(x_vector, 1);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;
      const double y = fitterData.points[i].second;

      gsl_vector_set(f, i, factor * std::pow(x, exponent) - y);
    }

    return GSL_SUCCESS;
  }

  static int FittingDerivative(const gsl_vector* x_vector, void* data,
                               gsl_matrix* j_matrix) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);
    const double exponent = gsl_vector_get(x_vector, 0);
    const double factor = gsl_vector_get(x_vector, 1);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;

      const double x_to_the_e = std::pow(x, exponent);
      const double dfd_exp = factor * std::log(x) * x_to_the_e;
      const double dfd_fac = x_to_the_e;

      gsl_matrix_set(j_matrix, i, 0, dfd_exp);
      gsl_matrix_set(j_matrix, i, 1, dfd_fac);
    }

    return GSL_SUCCESS;
  }

  static int FittingBoth(const gsl_vector* x_vector, void* data, gsl_vector* f,
                         gsl_matrix* j_matrix) {
    Fitting(x_vector, data, f);
    FittingDerivative(x_vector, data, j_matrix);
    return GSL_SUCCESS;
  }

  static int FittingSecondOrder(const gsl_vector* x_vector, void* data,
                                gsl_vector* f) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);
    const double exponent = gsl_vector_get(x_vector, 0);
    const double b = gsl_vector_get(x_vector, 1);
    const double c = gsl_vector_get(x_vector, 2);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;
      const double y = fitterData.points[i].second;

      // f(x) = (bx + cx^2)^a
      gsl_vector_set(f, i, pow(b * x + c * x * x, exponent) - y);
    }

    return GSL_SUCCESS;
  }

  static int FittingSecondOrderDerivative(const gsl_vector* x_vector,
                                          void* data, gsl_matrix* j_matrix) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);
    const double a = gsl_vector_get(x_vector, 0);
    const double b = gsl_vector_get(x_vector, 1);
    const double c = gsl_vector_get(x_vector, 2);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;

      // f(x)    = (bx + cx^2)^a
      // f(x)/da = ln(bx + cx^2) (bx + cx^2)^a
      // f(x)/db =    ax (bx + cx^2)^(a-1)
      // f(x)/dc =  ax^2 (bx + cx^2)^(a-1)
      const double inner_term = b * x + c * x * x;
      const double to_the_a = std::pow(inner_term, a);
      const double dfd_exp = std::log(inner_term) * to_the_a;
      const double to_the_a_minus_one = to_the_a / inner_term;
      const double dfd_factor_b = a * x * to_the_a_minus_one;
      const double dfd_factor_c = a * x * x * to_the_a_minus_one;

      gsl_matrix_set(j_matrix, i, 0, dfd_exp);
      gsl_matrix_set(j_matrix, i, 1, dfd_factor_b);
      gsl_matrix_set(j_matrix, i, 2, dfd_factor_c);
    }

    return GSL_SUCCESS;
  }

  static int FittingSecondOrderBoth(const gsl_vector* x, void* data,
                                    gsl_vector* f, gsl_matrix* j_matrix) {
    FittingSecondOrder(x, data, f);
    FittingSecondOrderDerivative(x, data, j_matrix);
    return GSL_SUCCESS;
  }

  static int FittingMultiOrder(const gsl_vector* x_vector, void* data,
                               gsl_vector* f) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;
      const double y = fitterData.points[i].second;
      const double lg = std::log10(x);

      // Horner's method
      double fity = 0.0;
      for (size_t k = 0; k != fitterData.n_terms; ++k) {
        const size_t j = fitterData.n_terms - k - 1;
        const double a_j = gsl_vector_get(x_vector, j);
        fity = a_j + fity * lg;
      }
      gsl_vector_set(f, i, std::pow(10.0, fity) - y);
    }

    return GSL_SUCCESS;
  }

  static int FittingMultiOrderDerivative(const gsl_vector* x_vector, void* data,
                                         gsl_matrix* j_matrix) {
    const NLPLFitterData& fitterData = *reinterpret_cast<NLPLFitterData*>(data);

    for (size_t i = 0; i != fitterData.points.size(); ++i) {
      const double x = fitterData.points[i].first;
      const double lg = std::log10(x);

      // Horner's method
      double fity = 0.0;
      for (size_t k = 0; k != fitterData.n_terms; ++k) {
        const size_t j = fitterData.n_terms - k - 1;
        const double a_j = gsl_vector_get(x_vector, j);
        fity = a_j + fity * lg;
      }
      fity = std::pow(10.0, fity);
      // dY/da_i = e^[ a_0...a_i-1,a_i+1...a_n] * (e^[a_i {log x}^i]) {log x}^i
      gsl_matrix_set(j_matrix, i, 0, fity);

      double lgPower = lg;
      for (size_t j = 1; j != fitterData.n_terms; ++j) {
        gsl_matrix_set(j_matrix, i, j, fity * lgPower);
        lgPower *= lg;
      }
    }

    return GSL_SUCCESS;
  }

  static int FittingMultiOrderBoth(const gsl_vector* x, void* data,
                                   gsl_vector* f, gsl_matrix* j_matrix) {
    FittingMultiOrder(x, data, f);
    FittingMultiOrderDerivative(x, data, j_matrix);
    return GSL_SUCCESS;
  }
};

void NonLinearPowerLawFitter::Fit(NumT& exponent, NumT& factor) {
  if (data_->points.size() >= 2) {
    const gsl_multifit_fdfsolver_type* t = gsl_multifit_fdfsolver_lmsder;
    data_->solver = gsl_multifit_fdfsolver_alloc(t, data_->points.size(), 2);

    gsl_multifit_function_fdf fdf;
    fdf.f = &NLPLFitterData::Fitting;
    fdf.df = &NLPLFitterData::FittingDerivative;
    fdf.fdf = &NLPLFitterData::FittingBoth;
    fdf.n = data_->points.size();
    fdf.p = 2;
    fdf.params = data_.get();

    double initial_values_array[2] = {exponent, factor};
    gsl_vector_view initial_values =
        gsl_vector_view_array(initial_values_array, 2);
    gsl_multifit_fdfsolver_set(data_->solver, &fdf, &initial_values.vector);

    int status;
    size_t iter = 0;
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(data_->solver);

      if (status) break;

      status = gsl_multifit_test_delta(data_->solver->dx, data_->solver->x,
                                       1.0e-7, 1.0e-7);

    } while (status == GSL_CONTINUE && iter < 500);

    exponent = gsl_vector_get(data_->solver->x, 0);
    factor = gsl_vector_get(data_->solver->x, 1);

    gsl_multifit_fdfsolver_free(data_->solver);
  } else {
    exponent = 0.0;
    factor = 0.0;
    for (size_t i = 0; i != data_->points.size(); ++i) {
      factor += data_->points[i].second;
    }
    factor /= NumT(data_->points.size());
  }
}

void NonLinearPowerLawFitter::Fit(NumT& a, NumT& b, NumT& c) {
  Fit(a, b);
  b = pow(b, 1.0 / a);

  if (data_->points.size() >= 3) {
    const gsl_multifit_fdfsolver_type* t = gsl_multifit_fdfsolver_lmsder;
    data_->solver = gsl_multifit_fdfsolver_alloc(t, data_->points.size(), 3);

    gsl_multifit_function_fdf fdf;
    fdf.f = &NLPLFitterData::FittingSecondOrder;
    fdf.df = &NLPLFitterData::FittingSecondOrderDerivative;
    fdf.fdf = &NLPLFitterData::FittingSecondOrderBoth;
    fdf.n = data_->points.size();
    fdf.p = 3;
    fdf.params = data_.get();

    double initial_values_array[3] = {a, b, c};
    gsl_vector_view initial_values =
        gsl_vector_view_array(initial_values_array, 3);
    gsl_multifit_fdfsolver_set(data_->solver, &fdf, &initial_values.vector);

    int status;
    size_t iter = 0;
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(data_->solver);

      if (status) break;

      status = gsl_multifit_test_delta(data_->solver->dx, data_->solver->x,
                                       1.0e-7, 1.0e-7);

    } while (status == GSL_CONTINUE && iter < 500);

    a = gsl_vector_get(data_->solver->x, 0);
    b = gsl_vector_get(data_->solver->x, 1);
    c = gsl_vector_get(data_->solver->x, 2);

    gsl_multifit_fdfsolver_free(data_->solver);
  }
}

void NonLinearPowerLawFitter::FitImplementation(std::vector<NumT>& terms,
                                                size_t n_terms) {
  data_->n_terms = n_terms;
  const gsl_multifit_fdfsolver_type* t = gsl_multifit_fdfsolver_lmsder;
  data_->solver =
      gsl_multifit_fdfsolver_alloc(t, data_->points.size(), n_terms);

  gsl_multifit_function_fdf fdf;
  fdf.f = &NLPLFitterData::FittingMultiOrder;
  fdf.df = &NLPLFitterData::FittingMultiOrderDerivative;
  fdf.fdf = &NLPLFitterData::FittingMultiOrderBoth;
  fdf.n = data_->points.size();
  fdf.p = n_terms;
  fdf.params = data_.get();

  std::vector<double> termsView(terms.begin(), terms.end());
  gsl_vector_view initial_values =
      gsl_vector_view_array(termsView.data(), n_terms);
  gsl_multifit_fdfsolver_set(data_->solver, &fdf, &initial_values.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(data_->solver);

    if (status) break;

    status = gsl_multifit_test_delta(data_->solver->dx, data_->solver->x,
                                     1.0e-6, 1.0e-6);

  } while (status == GSL_CONTINUE && iter < 5000);

  if (status != GSL_SUCCESS) {
    std::cout << "Warning: not converged! (niter=" << iter
              << ", status=" << gsl_strerror(status) << ")\n";
  }
  for (size_t i = 0; i != n_terms; ++i) {
    terms[i] = gsl_vector_get(data_->solver->x, i);
  }

  gsl_multifit_fdfsolver_free(data_->solver);
}

NonLinearPowerLawFitter::NonLinearPowerLawFitter()
    : data_(new NLPLFitterData()) {}

NonLinearPowerLawFitter::~NonLinearPowerLawFitter() = default;

void NonLinearPowerLawFitter::AddDataPoint(NumT x, NumT y) {
  data_->points.emplace_back(x, y);
}

void NonLinearPowerLawFitter::Fit(std::vector<NumT>& terms, size_t n_terms) {
  terms.assign(n_terms, 0.0);

  if (data_->points.size() < n_terms) n_terms = data_->points.size();

  if (n_terms == 0) return;

  NumT a = 1.0;
  NumT b = 0.0;
  Fit(a, b);
  bool is_negative = b < 0.0;
  if (is_negative) {
    for (std::pair<double, double>& point : data_->points) {
      point.second = -point.second;
    }
    terms[0] = -std::log10(-b);
    a = -a;
  } else {
    terms[0] = std::log10(b);  // - a*log(NLPLFact);
  }

  if (b != 0.0) {
    if (n_terms > 1) terms[1] = a;

    FitImplementation(terms, n_terms);
  }

  if (is_negative) {
    terms[0] = -std::pow(10.0, terms[0]);
  } else {
    terms[0] = std::pow(10.0, terms[0]);
  }
}

void NonLinearPowerLawFitter::FitStable(std::vector<NumT>& terms,
                                        size_t n_terms) {
  terms.assign(n_terms, 0.0);
  if (n_terms == 0) return;

  NumT a = 1.0;
  NumT b = 0.0;
  Fit(a, b);

  bool is_negative = b < 0.0;
  if (is_negative) {
    for (std::pair<double, double>& point : data_->points) {
      point.second = -point.second;
    }
    terms[0] = -std::log10(-b);
    a = -a;
  } else {
    terms[0] = std::log10(b);  // - a*log(NLPLFact);
  }

  if (b != 0.0) {
    if (n_terms > 1) terms[1] = a;
    size_t n_terms_estimated = 2;
    while (n_terms_estimated < n_terms) {
      ++n_terms_estimated;
      FitImplementation(terms, n_terms_estimated);
    }
  }

  if (is_negative) {
    terms[0] = -std::pow(10.0, terms[0]);
  } else {
    terms[0] = std::pow(10.0, terms[0]);
  }
}

void NonLinearPowerLawFitter::FastFit(NumT& exponent, NumT& factor) {
  NumT sumxy = 0.0;
  NumT sumx = 0.0;
  NumT sumy = 0.0;
  NumT sumxx = 0.0;
  size_t n = 0;
  bool require_non_linear = false;

  for (const std::pair<double, double>& point : data_->points) {
    const NumT x = point.first;
    const NumT y = point.second;
    if (y <= 0) {
      require_non_linear = true;
      break;
    }
    if (x > 0 && y > 0) {
      const long double log_x = std::log(x);
      const long double log_y = std::log(y);
      sumxy += log_x * log_y;
      sumx += log_x;
      sumy += log_y;
      sumxx += log_x * log_x;
      ++n;
    }
  }
  if (require_non_linear) {
    exponent = 0.0;
    factor = 1.0;
    Fit(exponent, factor);
  } else {
    if (n == 0) {
      exponent = std::numeric_limits<NumT>::quiet_NaN();
      factor = std::numeric_limits<NumT>::quiet_NaN();
    } else {
      const NumT d = n * sumxx - sumx * sumx;
      if (d == 0.0) {
        exponent = 0.0;
      } else {
        exponent = (n * sumxy - sumx * sumy) / d;
      }
      factor = std::exp((sumy - exponent * sumx) / n);
    }
  }
}

NonLinearPowerLawFitter::NumT NonLinearPowerLawFitter::Evaluate(
    NumT x, const std::vector<NumT>& terms, NumT reference_frequency_hz) {
  if (terms.empty()) return 0.0;
  NumT y = 0.0;
  const NumT lg = std::log10(x / reference_frequency_hz);

  for (size_t k = 0; k != terms.size() - 1; ++k) {
    const size_t j = terms.size() - k - 1;
    y = y * lg + terms[j];
  }
  return std::pow(10.0, y * lg) * terms[0];
}

}  // namespace fitters
}  // namespace schaapcommon