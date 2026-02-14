#ifndef AOCOMMON_SCALAR_EIGENVALUES_H_
#define AOCOMMON_SCALAR_EIGENVALUES_H_

#include <cmath>
#include <complex>

namespace aocommon {

/**
 * Compute eigen values and vectors for real matrix. Assumes
 * the determinant > 0.
 */
inline void EigenValuesAndVectors(const double* matrix, double& e1, double& e2,
                                  double* vec1, double* vec2) {
  double tr = matrix[0] + matrix[3];
  double d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
  double term = std::sqrt(tr * tr * 0.25 - d);
  double trHalf = tr * 0.5;
  e1 = trHalf + term;
  e2 = trHalf - term;
  double limit = std::min(std::fabs(e1), std::fabs(e2)) * 1e-6;
  if (std::fabs(matrix[2]) > limit) {
    vec1[0] = matrix[3] - e1;
    vec1[1] = -matrix[2];
    vec2[0] = matrix[3] - e2;
    vec2[1] = -matrix[2];
  } else if (std::fabs(matrix[1]) > limit) {
    vec1[0] = -matrix[1];
    vec1[1] = matrix[0] - e1;
    vec2[0] = -matrix[1];
    vec2[1] = matrix[0] - e2;
  } else {
    // We know that A v = lambda v, and we know that v1 or v2 = [1, 0]:
    double
        // Evaluate for v = [1, 0] and see if the error is smaller for e1 than
        // for e2
        err1_0 = matrix[0] - e1,
        err2_0 = matrix[0] - e2;
    if (err1_0 * err1_0 < err2_0 * err2_0) {
      vec1[0] = 1.0;
      vec1[1] = 0.0;
      vec2[0] = 0.0;
      vec2[1] = 1.0;
    } else {
      vec1[0] = 0.0;
      vec1[1] = 1.0;
      vec2[0] = 1.0;
      vec2[1] = 0.0;
    }
  }
}

/**
 * Compute eigen values and vectors for complex-valued matrix. Assumes
 * the determinant > 0.
 *
 * TODO: can probably be merged with previous function
 */
inline void EigenValuesAndVectors(const std::complex<double>* matrix,
                                  std::complex<double>& e1,
                                  std::complex<double>& e2,
                                  std::complex<double>* vec1,
                                  std::complex<double>* vec2) {
  std::complex<double> tr = matrix[0] + matrix[3];
  std::complex<double> d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
  std::complex<double> term = std::sqrt(tr * tr * 0.25 - d);
  std::complex<double> trHalf = tr * 0.5;
  e1 = trHalf + term;
  e2 = trHalf - term;
  double limit = std::min(std::abs(e1), std::abs(e2)) * 1e-6;
  if (std::abs(matrix[2]) > limit) {
    vec1[0] = matrix[3] - e1;
    vec1[1] = -matrix[2];
    vec2[0] = matrix[3] - e2;
    vec2[1] = -matrix[2];
  } else if (std::abs(matrix[1]) > limit) {
    vec1[0] = -matrix[1];
    vec1[1] = matrix[0] - e1;
    vec2[0] = -matrix[1];
    vec2[1] = matrix[0] - e2;
  } else {
    // We know that A v = lambda v, and we know that v1 or v2 = [1, 0]:
    auto
        // Evaluate for v = [1, 0] and see if the error is smaller for e1 than
        // for e2
        err1_0 = std::norm(matrix[0] - e1),
        err2_0 = std::norm(matrix[0] - e2);
    if (err1_0 < err2_0) {
      vec1[0] = 1.0;
      vec1[1] = 0.0;
      vec2[0] = 0.0;
      vec2[1] = 1.0;
    } else {
      vec1[0] = 0.0;
      vec1[1] = 1.0;
      vec2[0] = 1.0;
      vec2[1] = 0.0;
    }
  }
}

}  // namespace aocommon

#endif
