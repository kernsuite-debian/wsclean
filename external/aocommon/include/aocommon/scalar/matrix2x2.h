#ifndef AOCOMMON_SCALAR_MATRIX_2X2_H_
#define AOCOMMON_SCALAR_MATRIX_2X2_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>
#include <sstream>

#include "eigenvalues.h"
#include "vector4.h"

namespace aocommon::scalar {

template <typename T>
class MC2x2DiagBase;

/**
 * Class wraps functionality around a size 4 pointer
 * as if it were a 2x2 matrix.
 *
 */
class Matrix2x2 {
 public:
  /**
   * Copy complex-valued source buffer to complex-valued dest buffer.
   *
   * TODO: seems redundant?
   */
  template <typename LHS_T, typename RHS_T>
  static void Assign(std::complex<LHS_T>* dest,
                     const std::complex<RHS_T>* source) {
    for (size_t p = 0; p != 4; ++p) dest[p] = source[p];
  }

  /**
   * Copy source buffer to dest buffer.
   *
   */
  template <typename LHS_T, typename RHS_T>
  static void Assign(LHS_T* dest, const RHS_T* source) {
    for (size_t p = 0; p != 4; ++p) dest[p] = source[p];
  }

  /**
   * Add assign rhs buffer to complex-valued dest buffer.
   */
  template <typename T, typename RHS_T>
  static void Add(std::complex<T>* dest, const RHS_T* rhs) {
    for (size_t p = 0; p != 4; ++p) dest[p] += rhs[p];
  }

  /**
   * Subtract assign complex-valued rhs buffer to complex-valued dest buffer.
   * Assumes that T and RHS_T admit an implicit conversion.
   *
   */
  template <typename T>
  static void Subtract(std::complex<T>* dest, const std::complex<T>* rhs) {
    for (size_t p = 0; p != 4; ++p) dest[p] -= rhs[p];
  }

  /**
   * Check if all entries in matrix are finite
   */
  template <typename T>
  static bool IsFinite(const std::complex<T>* matrix) {
    return std::isfinite(matrix[0].real()) && std::isfinite(matrix[0].imag()) &&
           std::isfinite(matrix[1].real()) && std::isfinite(matrix[1].imag()) &&
           std::isfinite(matrix[2].real()) && std::isfinite(matrix[2].imag()) &&
           std::isfinite(matrix[3].real()) && std::isfinite(matrix[3].imag());
  }

  /**
   * Scalar multiplication of matrix.
   */
  template <typename LHS_T, typename RHS_T>
  static void ScalarMultiply(LHS_T* dest, RHS_T factor) {
    for (size_t p = 0; p != 4; ++p) dest[p] *= factor;
  }

  /**
   * Multiply rhs matrix with factor, then add assign to lhs matrix
   */
  template <typename T, typename RHS, typename FactorType>
  static void MultiplyAdd(std::complex<T>* dest, const RHS* rhs,
                          FactorType factor) {
    for (size_t p = 0; p != 4; ++p) dest[p] += rhs[p] * factor;
  }

  /**
   * Matrix multiplication
   */
  template <typename ComplType, typename LHS_T, typename RHS_T>
  static void ATimesB(std::complex<ComplType>* dest, const LHS_T* lhs,
                      const RHS_T* rhs) {
    dest[0] = lhs[0] * rhs[0] + lhs[1] * rhs[2];
    dest[1] = lhs[0] * rhs[1] + lhs[1] * rhs[3];
    dest[2] = lhs[2] * rhs[0] + lhs[3] * rhs[2];
    dest[3] = lhs[2] * rhs[1] + lhs[3] * rhs[3];
  }

  /**
   * Add assign matrix multiplication to destination buffer
   *
   * TODO: use templated type?
   */
  static void PlusATimesB(std::complex<double>* dest,
                          const std::complex<double>* lhs,
                          const std::complex<double>* rhs) {
    dest[0] += lhs[0] * rhs[0] + lhs[1] * rhs[2];
    dest[1] += lhs[0] * rhs[1] + lhs[1] * rhs[3];
    dest[2] += lhs[2] * rhs[0] + lhs[3] * rhs[2];
    dest[3] += lhs[2] * rhs[1] + lhs[3] * rhs[3];
  }

  /**
   * Matrix multiplication of matrix A with the Hermitian transpose of matrix B,
   * i.e. result = A * B^H
   *
   * TODO: seems unnecessary to use three templated types
   */
  template <typename ComplType, typename LHS_T, typename RHS_T>
  static void ATimesHermB(std::complex<ComplType>* dest, const LHS_T* lhs,
                          const RHS_T* rhs) {
    dest[0] = lhs[0] * std::conj(rhs[0]) + lhs[1] * std::conj(rhs[1]);
    dest[1] = lhs[0] * std::conj(rhs[2]) + lhs[1] * std::conj(rhs[3]);
    dest[2] = lhs[2] * std::conj(rhs[0]) + lhs[3] * std::conj(rhs[1]);
    dest[3] = lhs[2] * std::conj(rhs[2]) + lhs[3] * std::conj(rhs[3]);
  }

  /**
   * Add assign matrix multiplication of matrix A with the
   * Hermitian transpose of matrix B, i.e. result += A * B^H
   *
   * TODO: seems unnecessary to use three templated types
   */
  template <typename ComplType, typename LHS_T, typename RHS_T>
  static void PlusATimesHermB(std::complex<ComplType>* dest, const LHS_T* lhs,
                              const RHS_T* rhs) {
    dest[0] += lhs[0] * std::conj(rhs[0]) + lhs[1] * std::conj(rhs[1]);
    dest[1] += lhs[0] * std::conj(rhs[2]) + lhs[1] * std::conj(rhs[3]);
    dest[2] += lhs[2] * std::conj(rhs[0]) + lhs[3] * std::conj(rhs[1]);
    dest[3] += lhs[2] * std::conj(rhs[2]) + lhs[3] * std::conj(rhs[3]);
  }

  /**
   * Matrix multiplication of the Hermitian transpose of matrix A with matrix B,
   * i.e. A^H * B
   */
  template <typename ComplType, typename LHS_T, typename RHS_T>
  static void HermATimesB(std::complex<ComplType>* dest, const LHS_T* lhs,
                          const RHS_T* rhs) {
    dest[0] = std::conj(lhs[0]) * rhs[0] + std::conj(lhs[2]) * rhs[2];
    dest[1] = std::conj(lhs[0]) * rhs[1] + std::conj(lhs[2]) * rhs[3];
    dest[2] = std::conj(lhs[1]) * rhs[0] + std::conj(lhs[3]) * rhs[2];
    dest[3] = std::conj(lhs[1]) * rhs[1] + std::conj(lhs[3]) * rhs[3];
  }

  /**
   * Matrix multiplication of the Hermitian transpose of matrix A with Hermitian
   * transpose of B, i.e. A^H * B^H
   */
  static void HermATimesHermB(std::complex<double>* dest,
                              const std::complex<double>* lhs,
                              const std::complex<double>* rhs) {
    dest[0] = std::conj(lhs[0]) * std::conj(rhs[0]) +
              std::conj(lhs[2]) * std::conj(rhs[1]);
    dest[1] = std::conj(lhs[0]) * std::conj(rhs[2]) +
              std::conj(lhs[2]) * std::conj(rhs[3]);
    dest[2] = std::conj(lhs[1]) * std::conj(rhs[0]) +
              std::conj(lhs[3]) * std::conj(rhs[1]);
    dest[3] = std::conj(lhs[1]) * std::conj(rhs[2]) +
              std::conj(lhs[3]) * std::conj(rhs[3]);
  }

  /**
   * Add assign matrix multiplication A^H * B to the destination buffer.
   *
   * TODO: seems redundant to template three types
   */
  template <typename ComplType, typename LHS_T, typename RHS_T>
  static void PlusHermATimesB(std::complex<ComplType>* dest, const LHS_T* lhs,
                              const RHS_T* rhs) {
    dest[0] += std::conj(lhs[0]) * rhs[0] + std::conj(lhs[2]) * rhs[2];
    dest[1] += std::conj(lhs[0]) * rhs[1] + std::conj(lhs[2]) * rhs[3];
    dest[2] += std::conj(lhs[1]) * rhs[0] + std::conj(lhs[3]) * rhs[2];
    dest[3] += std::conj(lhs[1]) * rhs[1] + std::conj(lhs[3]) * rhs[3];
  }

  /**
   * Compute matrix inverse
   */
  template <typename T>
  static bool Invert(T* matrix) {
    T d = ((matrix[0] * matrix[3]) - (matrix[1] * matrix[2]));
    if (d == T(0.0)) return false;
    T determinant_reciprocal = T(1.0) / d;
    T temp;
    temp = matrix[3] * determinant_reciprocal;
    matrix[1] = -matrix[1] * determinant_reciprocal;
    matrix[2] = -matrix[2] * determinant_reciprocal;
    matrix[3] = matrix[0] * determinant_reciprocal;
    matrix[0] = temp;
    return true;
  }

  /**
   * Compute conjugate transpose (a.k.a. Hermitian transpose) of matrix
   */
  template <typename T>
  static void ConjugateTranspose(T* matrix) {
    matrix[0] = std::conj(matrix[0]);
    T temp = matrix[1];
    matrix[1] = std::conj(matrix[2]);
    matrix[2] = std::conj(temp);
    matrix[3] = std::conj(matrix[3]);
  }

  /**
   * Multiply lhs buffer with inverse of rhs buffer. Returns false if
   * rhs not invertible.
   */
  static bool MultiplyWithInverse(std::complex<double>* lhs,
                                  const std::complex<double>* rhs) {
    std::complex<double> d = ((rhs[0] * rhs[3]) - (rhs[1] * rhs[2]));
    if (d == 0.0) return false;
    std::complex<double> determinant_reciprocal = 1.0 / d;
    std::complex<double> temp[4];
    temp[0] = rhs[3] * determinant_reciprocal;
    temp[1] = -rhs[1] * determinant_reciprocal;
    temp[2] = -rhs[2] * determinant_reciprocal;
    temp[3] = rhs[0] * determinant_reciprocal;

    std::complex<double> temp2 = lhs[0];
    lhs[0] = lhs[0] * temp[0] + lhs[1] * temp[2];
    lhs[1] = temp2 * temp[1] + lhs[1] * temp[3];

    temp2 = lhs[2];
    lhs[2] = lhs[2] * temp[0] + lhs[3] * temp[2];
    lhs[3] = temp2 * temp[1] + lhs[3] * temp[3];
    return true;
  }

  /**
   * Compute singular values of the matrix buffer
   */
  static void SingularValues(const std::complex<double>* matrix, double& e1,
                             double& e2) {
    // This is not the ultimate fastest method, since we
    // don't need to calculate the imaginary values of b,c at all.
    // Calculate M M^H
    std::complex<double> temp[4] = {
        matrix[0] * std::conj(matrix[0]) + matrix[1] * std::conj(matrix[1]),
        matrix[0] * std::conj(matrix[2]) + matrix[1] * std::conj(matrix[3]),
        matrix[2] * std::conj(matrix[0]) + matrix[3] * std::conj(matrix[1]),
        matrix[2] * std::conj(matrix[2]) + matrix[3] * std::conj(matrix[3])};
    // Use quadratic formula, with a=1.
    double b = -temp[0].real() - temp[3].real(),
           c = temp[0].real() * temp[3].real() - (temp[1] * temp[2]).real(),
           d = b * b - (4.0 * 1.0) * c, sqrtd = std::sqrt(d);

    e1 = std::sqrt((-b + sqrtd) * 0.5);
    e2 = std::sqrt((-b - sqrtd) * 0.5);
  }

  /**
   * Compute eigen values of input matrix buffer. It assumes that the
   * determinant > 0, so that eigenvalues are real.
   */
  static void EigenValues(const double* matrix, double& e1, double& e2) {
    double tr = matrix[0] + matrix[3];
    double d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    double term = std::sqrt(tr * tr * 0.25 - d);
    double trHalf = tr * 0.5;
    e1 = trHalf + term;
    e2 = trHalf - term;
  }

  /**
   * Compute the eigen values of a complex matrix.
   *
   * TODO: can probably be merged with previous method.
   */
  template <typename ValType>
  static void EigenValues(const std::complex<ValType>* matrix,
                          std::complex<ValType>& e1,
                          std::complex<ValType>& e2) {
    std::complex<ValType> tr = matrix[0] + matrix[3];
    std::complex<ValType> d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    std::complex<ValType> term = std::sqrt(tr * tr * ValType(0.25) - d);
    std::complex<ValType> trHalf = tr * ValType(0.5);
    e1 = trHalf + term;
    e2 = trHalf - term;
  }

  /**
   * Compute eigen values and vectors for real matrix. Assumes
   * the determinant > 0.
   */
  static void EigenValuesAndVectors(const double* matrix, double& e1,
                                    double& e2, double* vec1, double* vec2) {
    aocommon::EigenValuesAndVectors(matrix, e1, e2, vec1, vec2);
  }

  /**
   * Compute eigen values and vectors for complex-valued matrix. Assumes
   * the determinant > 0.
   *
   * TODO: can probably be merged with previous method
   */
  static void EigenValuesAndVectors(const std::complex<double>* matrix,
                                    std::complex<double>& e1,
                                    std::complex<double>& e2,
                                    std::complex<double>* vec1,
                                    std::complex<double>* vec2) {
    aocommon::EigenValuesAndVectors(matrix, e1, e2, vec1, vec2);
  }

  /**
   * Computes the positive square root of a real-valued matrix buffer such that
   * M = R * R. Assumes that determinant > 0. Note that matrix M might have more
   * square roots.
   */
  static void SquareRoot(double* matrix) {
    double tr = matrix[0] + matrix[3];
    double d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    double s = /*+/-*/ std::sqrt(d);
    double t = /*+/-*/ std::sqrt(tr + 2.0 * s);
    if (t != 0.0) {
      matrix[0] = (matrix[0] + s) / t;
      matrix[1] = (matrix[1] / t);
      matrix[2] = (matrix[2] / t);
      matrix[3] = (matrix[3] + s) / t;
    } else {
      if (matrix[0] == 0.0 && matrix[1] == 0.0 && matrix[2] == 0.0 &&
          matrix[3] == 0.0) {
        // done: it's the zero matrix
      } else {
        for (size_t i = 0; i != 4; ++i)
          matrix[i] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  /**
   * Computes the positive square root of a complex-valued matrix buffer,
   * such that M = R * R. Assumes that determinant > 0.
   * Note that matrix M might have more square roots.
   */
  static void SquareRoot(std::complex<double>* matrix) {
    std::complex<double> tr = matrix[0] + matrix[3];
    std::complex<double> d = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    std::complex<double> s = /*+/-*/ std::sqrt(d);
    std::complex<double> t = /*+/-*/ std::sqrt(tr + 2.0 * s);
    if (t != 0.0) {
      matrix[0] = (matrix[0] + s) / t;
      matrix[1] = (matrix[1] / t);
      matrix[2] = (matrix[2] / t);
      matrix[3] = (matrix[3] + s) / t;
    } else {
      if (matrix[0] == 0.0 && matrix[1] == 0.0 && matrix[2] == 0.0 &&
          matrix[3] == 0.0) {
        // done: it's the zero matrix
      } else {
        for (size_t i = 0; i != 4; ++i)
          matrix[i] =
              std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN());
      }
    }
  }

  /**
   * Calculates L, the lower triangle of the Cholesky decomposition, such that
   * L L^H = M. The result is undefined when the matrix is not positive
   * definite.
   */
  static void UncheckedCholesky(std::complex<double>* matrix) {
    // solve:
    // ( a 0 ) ( a* b* ) = ( aa* ;    ab*    )
    // ( b c ) ( 0  c* )   ( a*b ; bb* + cc* )
    // With a and c necessarily real.
    double a = std::sqrt(matrix[0].real());
    std::complex<double> b = std::conj(matrix[1] / a);
    double bbConj = b.real() * b.real() + b.imag() * b.imag();
    double c = std::sqrt(matrix[3].real() - bbConj);
    matrix[0] = a;
    matrix[1] = 0.0;
    matrix[2] = b;
    matrix[3] = c;
  }

  /**
   * Calculates L, the lower triangle of the Cholesky decomposition, such that
   * L L^H = M. Return false when the result would not be finite.
   */
  static bool Cholesky(std::complex<double>* matrix) {
    if (matrix[0].real() < 0.0) return false;
    double a = std::sqrt(matrix[0].real());
    std::complex<double> b = std::conj(matrix[1] / a);
    double bbConj = b.real() * b.real() + b.imag() * b.imag();
    double cc = matrix[3].real() - bbConj;
    if (cc < 0.0) return false;
    double c = std::sqrt(cc);
    matrix[0] = a;
    matrix[1] = 0.0;
    matrix[2] = b;
    matrix[3] = c;
    return true;
  }

  /**
   * Calculates L, the lower triangle of the Cholesky decomposition, such that
   * L L^H = M. Return false when the matrix was not positive semi-definite.
   */
  static bool CheckedCholesky(std::complex<double>* matrix) {
    if (matrix[0].real() <= 0.0 || matrix[0].imag() != 0.0 ||
        matrix[3].real() <= 0.0 || matrix[3].imag() != 0.0 ||
        matrix[1] != std::conj(matrix[2]))
      return false;
    UncheckedCholesky(matrix);
    return true;
  }

  /**
   * Calculates the rotation angle of a complex-valued matrix.
   */
  template <typename T>
  static T RotationAngle(const std::complex<T>* matrix) {
    return std::atan2((matrix[2].real() - matrix[1].real()) * 0.5,
                      (matrix[0].real() + matrix[3].real()) * 0.5);
  }

  /**
   * Calculates the rotation matrix, given a rotation angle \p alpha.
   */
  template <typename T>
  static void RotationMatrix(std::complex<T>* matrix, double alpha) {
    T cos_alpha = std::cos(alpha), sin_alpha = std::sin(alpha);
    matrix[0] = cos_alpha;
    matrix[1] = -sin_alpha;
    matrix[2] = sin_alpha;
    matrix[3] = cos_alpha;
  }
};

/**
 * Class implements a 2x2 complex-valued matrix.
 */
template <typename ValType>
class MC2x2Base {
 public:
  /**
   * Creates an uninitialized 2x2 matrix. Be aware that the complex
   * values are not initialized to zero (use Zero() for that).
   */
  MC2x2Base() {}

  /**
   * Copy constructor. Even though the template copy constructor below covers
   * this case, the compiler declares this copy constructor implicitly, which
   * is deprecated in C++11 -> Declare the copy constructor explicitly.
   */
  constexpr MC2x2Base(const MC2x2Base& source) = default;

  template <typename OtherValType>
  MC2x2Base(const MC2x2Base<OtherValType>& source)
      : values_{source.Get(0).real(), source.Get(0).imag(),
                source.Get(1).real(), source.Get(1).imag(),
                source.Get(2).real(), source.Get(2).imag(),
                source.Get(3).real(), source.Get(3).imag()} {}

  /**
   * Construct MC2x2Base object from (length 4) data buffer with real values.
   */
  explicit constexpr MC2x2Base(const float source[4])
      : values_{source[0], 0.0, source[1], 0.0,
                source[2], 0.0, source[3], 0.0} {}
  /**
   * Construct MC2x2Base object from (length 4) data buffer with complex values.
   */
  explicit constexpr MC2x2Base(const std::complex<float> source[4])
      : values_{source[0].real(), source[0].imag(), source[1].real(),
                source[1].imag(), source[2].real(), source[2].imag(),
                source[3].real(), source[3].imag()} {}
  /**
   * Construct MC2x2Base object from (length 4) data buffer with real values.
   */
  explicit constexpr MC2x2Base(const double source[4])
      : values_{ValType(source[0]), ValType(0.0),       ValType(source[1]),
                ValType(0.0),       ValType(source[2]), ValType(0.0),
                ValType(source[3]), ValType(0.0)} {}
  /**
   * Construct MC2x2Base object from (length 4) data buffer with complex values.
   */
  explicit constexpr MC2x2Base(const std::complex<double> source[4])
      : values_{ValType(source[0].real()), ValType(source[0].imag()),
                ValType(source[1].real()), ValType(source[1].imag()),
                ValType(source[2].real()), ValType(source[2].imag()),
                ValType(source[3].real()), ValType(source[3].imag())} {}

  /**
   * Construct MC2x2Base object from four real values. Internally, values are
   * converted to complex type.
   */
  constexpr MC2x2Base(ValType m00, ValType m01, ValType m10, ValType m11)
      : values_{m00, ValType(0.0), m01, ValType(0.0),
                m10, ValType(0.0), m11, ValType(0.0)} {}

  /**
   * Construct MC2x2Base object from four complex-valued input values.
   */
  constexpr MC2x2Base(std::complex<ValType> m00, std::complex<ValType> m01,
                      std::complex<ValType> m10, std::complex<ValType> m11)
      : values_{m00.real(), m00.imag(), m01.real(), m01.imag(),
                m10.real(), m10.imag(), m11.real(), m11.imag()} {}

  /**
   * Construct from a diagonal matrix
   */
  constexpr explicit MC2x2Base(const MC2x2DiagBase<ValType>& diag);

  /**
   * Construct from initializer list, values are internally converted
   * to complex type. Assumes that list has size four.
   */
  MC2x2Base(std::initializer_list<ValType> list) {
    assert(list.size() == 4);
    typename std::initializer_list<ValType>::const_iterator i = list.begin();
    values_[0] = *i;
    values_[1] = 0.0;
    ++i;
    values_[2] = *i;
    values_[3] = 0.0;
    ++i;
    values_[4] = *i;
    values_[5] = 0.0;
    ++i;
    values_[6] = *i;
    values_[7] = 0.0;
  }

  /**
   * Construct from initializer list. Assumes that list has size four.
   */
  MC2x2Base(std::initializer_list<std::complex<ValType>> list) {
    assert(list.size() == 4);
    typename std::initializer_list<std::complex<ValType>>::const_iterator i =
        list.begin();
    values_[0] = i->real();
    values_[1] = i->imag();
    ++i;
    values_[2] = i->real();
    values_[3] = i->imag();
    ++i;
    values_[4] = i->real();
    values_[5] = i->imag();
    ++i;
    values_[6] = i->real();
    values_[7] = i->imag();
  }

  MC2x2Base<ValType>& operator=(const MC2x2Base<ValType>& source) = default;

  MC2x2Base<ValType> operator+(const MC2x2Base<ValType>& rhs) const {
    MC2x2Base<ValType> result;
    std::transform(std::begin(values_), std::end(values_),
                   std::begin(rhs.values_), std::begin(result.values_),
                   std::plus<ValType>());
    return result;
  }

  MC2x2Base<ValType> operator+(const MC2x2DiagBase<ValType>& rhs) const {
    return MC2x2Base<ValType>(Get(0) + rhs.Get(0), Get(1), Get(2),
                              Get(3) + rhs.Get(1));
  }

  friend MC2x2Base<ValType> operator+(const MC2x2DiagBase<ValType>& lhs,
                                      const MC2x2Base<ValType>& rhs) {
    return MC2x2Base<ValType>(lhs.Get(0) + rhs.Get(0), rhs.Get(1), rhs.Get(2),
                              lhs.Get(1) + rhs.Get(3));
  }

  MC2x2Base<ValType>& operator+=(const MC2x2Base<ValType>& rhs) {
    std::transform(std::begin(values_), std::end(values_),
                   std::begin(rhs.values_), std::begin(values_),
                   std::plus<ValType>());
    return *this;
  }

  MC2x2Base<ValType>& operator+=(const MC2x2DiagBase<ValType>& rhs) {
    values_[0] += rhs.Get(0).real();
    values_[1] += rhs.Get(0).imag();
    values_[6] += rhs.Get(1).real();
    values_[7] += rhs.Get(1).imag();
    return *this;
  }

  MC2x2Base<ValType> operator-(const MC2x2Base<ValType>& rhs) const {
    MC2x2Base<ValType> result;
    std::transform(std::begin(values_), std::end(values_),
                   std::begin(rhs.values_), std::begin(result.values_),
                   std::minus<ValType>());
    return result;
  }

  MC2x2Base<ValType> operator-(const MC2x2DiagBase<ValType>& rhs) const {
    return MC2x2Base<ValType>(Get(0) - rhs.Get(0), Get(1), Get(2),
                              Get(3) - rhs.Get(1));
  }

  friend MC2x2Base<ValType> operator-(const MC2x2DiagBase<ValType>& lhs,
                                      const MC2x2Base<ValType>& rhs) {
    return MC2x2Base<ValType>(lhs.Get(0) - rhs.Get(0), -rhs.Get(1), -rhs.Get(2),
                              lhs.Get(1) - rhs.Get(3));
  }

  MC2x2Base<ValType>& operator-=(const MC2x2Base<ValType>& rhs) {
    std::transform(std::begin(values_), std::end(values_),
                   std::begin(rhs.values_), std::begin(values_),
                   std::minus<ValType>());
    return *this;
  }

  MC2x2Base<ValType>& operator-=(const MC2x2DiagBase<ValType>& rhs) {
    values_[0] -= rhs.Get(0).real();
    values_[1] -= rhs.Get(0).imag();
    values_[6] -= rhs.Get(1).real();
    values_[7] -= rhs.Get(1).imag();
    return *this;
  }

  MC2x2Base<ValType>& operator*=(const MC2x2Base<ValType>& rhs) {
    *this = *this * rhs;
    return *this;
  }

  MC2x2Base<ValType> operator*(const MC2x2Base<ValType>& rhs) const {
    return MC2x2Base<ValType>{Get(0) * rhs.Get(0) + Get(1) * rhs.Get(2),
                              Get(0) * rhs.Get(1) + Get(1) * rhs.Get(3),
                              Get(2) * rhs.Get(0) + Get(3) * rhs.Get(2),
                              Get(2) * rhs.Get(1) + Get(3) * rhs.Get(3)};
  }

  bool operator==(const MC2x2Base<ValType>& rhs) const {
    return values_ == rhs.values_;
  }

  /**
   * Matrix multiplication assignment operator given a length 4 rhs buffer
   * of possibly different type
   */
  template <typename T>
  MC2x2Base<ValType>& operator*=(const T* rhs) {
    *this = *this * rhs;
    return *this;
  }

  /**
   * Matrix multiplication given a length 4 rhs buffer of possibly different
   * type
   */
  template <typename T>
  MC2x2Base<ValType> operator*(const T* rhs) const {
    return *this * MC2x2Base<ValType>(rhs);
  }

  /**
   * Scalar multiplication assignment operator
   */
  MC2x2Base<ValType>& operator*=(ValType rhs) {
    *this = *this * rhs;
    return *this;
  }

  /**
   * Scalar multiplication operator
   */
  MC2x2Base<ValType> operator*(ValType rhs) const {
    return MC2x2Base<ValType>(Get(0) * rhs, Get(1) * rhs, Get(2) * rhs,
                              Get(3) * rhs);
  }

  friend MC2x2Base<ValType> operator*(ValType lhs,
                                      const MC2x2Base<ValType>& rhs) {
    return MC2x2Base<ValType>{rhs.Get(0) * lhs, rhs.Get(1) * lhs,
                              rhs.Get(2) * lhs, rhs.Get(3) * lhs};
  }

  /**
   * Complex scalar multiplication
   */
  MC2x2Base<ValType> operator*(std::complex<ValType> rhs) const {
    return MC2x2Base<ValType>{Get(0) * rhs, Get(1) * rhs, Get(2) * rhs,
                              Get(3) * rhs};
  }

  friend MC2x2Base<ValType> operator*(std::complex<ValType> lhs,
                                      const MC2x2Base<ValType>& rhs) {
    return MC2x2Base<ValType>{rhs.Get(0) * lhs, rhs.Get(1) * lhs,
                              rhs.Get(2) * lhs, rhs.Get(3) * lhs};
  }

  /**
   * Scalar division assignment operator
   */
  MC2x2Base<ValType>& operator/=(ValType rhs) {
    *this *= ValType(1.0) / rhs;
    return *this;
  }

  std::complex<ValType> Get(size_t index) const {
    return std::complex<ValType>(values_[index * 2], values_[index * 2 + 1]);
  }
  void Set(size_t index, std::complex<ValType> value) {
    values_[index * 2] = value.real();
    values_[index * 2 + 1] = value.imag();
  }

  /**
   * Indexes the values as real values instead of complex values. This means
   * that IndexReal(0) returns the real value of the first complex value,
   * IndexReal(1) returns the imaginary value of the first complex value, and
   * IndexReal(7) returns the imaginary value of the last complex value.
   * @param index should be such that 0 <= index < 8.
   */
  const ValType& IndexReal(size_t index) const { return values_[index]; }
  ValType& IndexReal(size_t index) { return values_[index]; }

  /**
   * Return MC2x2Base matrix filled with zeros
   */
  static MC2x2Base<ValType> Zero() {
    return MC2x2Base<ValType>(0.0, 0.0, 0.0, 0.0);
  }

  /**
   * Return 2x2 identity matrix
   */
  static MC2x2Base<ValType> Unity() { return MC2x2Base(1.0, 0.0, 0.0, 1.0); }

  /**
   * Return 2x2 matrix filled with NaN values
   */
  static MC2x2Base<ValType> NaN() {
    return MC2x2Base<ValType>(
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()),
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()),
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()),
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()));
  }

  /**
   * Assign data stored by 2x2 matrix to destination buffer
   */
  template <typename T>
  void AssignTo(std::complex<T>* destination) const {
    destination[0] = Get(0);
    destination[1] = Get(1);
    destination[2] = Get(2);
    destination[3] = Get(3);
  }

  /**
   * Flatten 2x2 matrix to length 4 vector
   */
  Vector4 Vec() const { return Vector4(Get(0), Get(2), Get(1), Get(3)); }

  /**
   * Matrix multiplication, alias for the overloaded * operator and thus equally
   * computationally efficient.
   */
  MC2x2Base<ValType> Multiply(const MC2x2Base<ValType>& rhs) const {
    return *this * rhs;
  }

  /**
   * Matrix multiplication of internal matrix with Hermitian transpose of input
   * matrix, i.e. returns A * B^H. Note that this is preferred over combining
   * operator* with HermTranspose() for computational efficiency.
   */
  MC2x2Base<ValType> MultiplyHerm(const MC2x2Base<ValType>& rhs) const {
    return MC2x2Base<ValType>{
        Get(0) * std::conj(rhs.Get(0)) + Get(1) * std::conj(rhs.Get(1)),
        Get(0) * std::conj(rhs.Get(2)) + Get(1) * std::conj(rhs.Get(3)),
        Get(2) * std::conj(rhs.Get(0)) + Get(3) * std::conj(rhs.Get(1)),
        Get(2) * std::conj(rhs.Get(2)) + Get(3) * std::conj(rhs.Get(3))};
  }

  /**
   * Matrix multiplication Hermitian transpose of internal matrix with input
   * matrix, i.e. returns A^H * B. Note that this is preferred over combining
   * operator* with HermTranspose() for computational efficiency.
   */
  MC2x2Base<ValType> HermThenMultiply(const MC2x2Base<ValType>& rhs) const {
    return MC2x2Base<ValType>{
        std::conj(Get(0)) * rhs.Get(0) + std::conj(Get(2)) * rhs.Get(2),
        std::conj(Get(0)) * rhs.Get(1) + std::conj(Get(2)) * rhs.Get(3),
        std::conj(Get(1)) * rhs.Get(0) + std::conj(Get(3)) * rhs.Get(2),
        std::conj(Get(1)) * rhs.Get(1) + std::conj(Get(3)) * rhs.Get(3)};
  }

  /**
   * Matrix multiplication of Hermitian transposes of both the internal matrix
   * and the input matrix, i.e. returns A^H * B^H. Note that this is preferred
   * over combining operator* with HermTranspose() for computational efficiency.
   */
  MC2x2Base<ValType> HermThenMultiplyHerm(const MC2x2Base<ValType>& rhs) const {
    return MC2x2Base<ValType>{std::conj(Get(0)) * std::conj(rhs(0)) +
                                  std::conj(Get(2)) * std::conj(rhs(1)),
                              std::conj(Get(0)) * std::conj(rhs(2)) +
                                  std::conj(Get(2)) * std::conj(rhs(3)),
                              std::conj(Get(1)) * std::conj(rhs(0)) +
                                  std::conj(Get(3)) * std::conj(rhs(1)),
                              std::conj(Get(1)) * std::conj(rhs(2)) +
                                  std::conj(Get(3)) * std::conj(rhs(3))};
  }

  /**
   * Computes the double dot, i.e. A:B (A_ij Bij)
   * See https://en.wikipedia.org/wiki/Dyadics#Double-dot_product
   */
  std::complex<ValType> DoubleDot(const MC2x2Base<ValType>& rhs) const {
    return Get(0) * rhs.Get(0) + Get(1) * rhs.Get(1) + Get(2) * rhs.Get(2) +
           Get(3) * rhs.Get(3);
  }

  /**
   * Multiply input matrix with factor, then add assign to stored matrix
   */
  template <typename FactorType>
  void AddWithFactorAndAssign(const MC2x2Base<ValType>& rhs,
                              FactorType factor) {
    *this += MC2x2Base<ValType>{rhs.Get(0) * factor, rhs.Get(1) * factor,
                                rhs.Get(2) * factor, rhs.Get(3) * factor};
  }

  /**
   * Compute (regular) transpose of matrix
   */
  MC2x2Base<ValType> Transpose() const {
    return MC2x2Base(Get(0), Get(2), Get(1), Get(3));
  }

  /**
   * Compute Hermitian transpose of matrix
   */
  MC2x2Base<ValType> HermTranspose() const {
    return MC2x2Base(std::conj(Get(0)), std::conj(Get(2)), std::conj(Get(1)),
                     std::conj(Get(3)));
  }

  /**
   * Compute the elementwise conjugate of the matrix (without transposing!)
   */
  MC2x2Base<ValType> Conjugate() const {
    return MC2x2Base(std::conj(Get(0)), std::conj(Get(1)), std::conj(Get(2)),
                     std::conj(Get(3)));
  }

  /**
   * @returns A^H * A.
   */
  MC2x2Base<ValType> HermitianSquare() const {
    const std::complex<ValType> m01 =
        std::conj(Get(0)) * Get(1) + std::conj(Get(2)) * Get(3);
    return MC2x2Base(std::norm(Get(0)) + std::norm(Get(2)), m01, std::conj(m01),
                     std::norm(Get(1)) + std::norm(Get(3)));
  }

  /**
   * Invert 2x2 matrix, returns false if matrix is not invertible
   */
  bool Invert() {
    std::complex<ValType> data[4]{Get(0), Get(1), Get(2), Get(3)};
    const bool result = Matrix2x2::Invert(data);
    *this = MC2x2Base<ValType>(data);
    return result;
  }

  /**
   * Convert matrix to pretty string
   */
  std::string ToString() const {
    std::stringstream str;
    str << Get(0) << ", " << Get(1) << "; " << Get(2) << ", " << Get(3);
    return str.str();
  }

  [[deprecated("Use AssignTo()")]] void CopyValues(
      std::complex<ValType>* values) const {
    AssignTo(values);
  }

  /**
   * Calculate eigen values
   */
  void EigenValues(std::complex<ValType>& e1, std::complex<ValType>& e2) const {
    const std::complex<double> data[4]{Get(0), Get(1), Get(2), Get(3)};
    Matrix2x2::EigenValues(data, e1, e2);
  }

  /**
   * Check if matrix entries are finite
   */
  bool IsFinite() const {
    return std::isfinite(values_[0]) && std::isfinite(values_[1]) &&
           std::isfinite(values_[2]) && std::isfinite(values_[3]) &&
           std::isfinite(values_[4]) && std::isfinite(values_[5]) &&
           std::isfinite(values_[6]) && std::isfinite(values_[7]);
  }
  /**
   * Calculates L, the lower triangle of the Cholesky decomposition, such that
   * L L^H = M.
   */
  bool Cholesky() {
    std::complex<ValType> data[4]{Get(0), Get(1), Get(2), Get(3)};
    const bool result = Matrix2x2::Cholesky(data);
    *this = MC2x2Base<ValType>(data);
    return result;
  }

  /**
   * See Matrix2x2::CheckedCholesky
   */
  bool CheckedCholesky() {
    std::complex<ValType> data[4]{Get(0), Get(1), Get(2), Get(3)};
    const bool result = Matrix2x2::CheckedCholesky(data);
    *this = MC2x2Base<ValType>(data);
    return result;
  }

  /**
   * See Matrix2x2::UncheckedCholesky
   */
  void UncheckedCholesky() {
    std::complex<ValType> data[4]{Get(0), Get(1), Get(2), Get(3)};
    Matrix2x2::UncheckedCholesky(data);
    *this = MC2x2Base<ValType>(data);
  }

  /**
   * Decompose a Hermitian matrix X into A A^H such that
   *   X = A A^H = U D D^H U^H
   *   with A = U D
   * where D D^H = E is a diagonal matrix
   *       with the eigen values of X, and U contains the eigen vectors.
   */
  MC2x2Base<ValType> DecomposeHermitianEigenvalue() const {
    const std::complex<ValType> data[4]{Get(0), Get(1), Get(2), Get(3)};
    std::complex<ValType> e1, e2, vec1[2], vec2[2];
    Matrix2x2::EigenValuesAndVectors(data, e1, e2, vec1, vec2);
    ValType v1norm = std::norm(vec1[0]) + std::norm(vec1[1]);
    vec1[0] /= std::sqrt(v1norm);
    vec1[1] /= std::sqrt(v1norm);
    ValType v2norm = std::norm(vec2[0]) + std::norm(vec2[1]);
    vec2[0] /= std::sqrt(v2norm);
    vec2[1] /= std::sqrt(v2norm);

    return MC2x2Base<ValType>(
        vec1[0] * std::sqrt(e1.real()), vec2[0] * std::sqrt(e2.real()),
        vec1[1] * std::sqrt(e1.real()), vec2[1] * std::sqrt(e2.real()));
  }

 private:
  /**
   * This is actually a std::complex<ValType>[4], but because std::complex
   * initializes its values by default, we use std::array<ValType, 8> as type.
   *
   * Be aware that casting a ValType to a std::complex<ValType> is not allowed
   * because of type-punning rules (the other direction is fine though: it
   * is allowed to cast a std::complex<ValType> to a ValType).
   * Operations that need a complex can use @ref Get(size_t).
   */
  std::array<ValType, 8> values_;
};

/**
 * Left shift operator to write the matrix to ostream
 */
template <typename ValType>
std::ostream& operator<<(std::ostream& output,
                         const MC2x2Base<ValType>& value) {
  output << "[{" << value.Get(0) << ", " << value.Get(1) << "}, {"
         << value.Get(2) << ", " << value.Get(3) << "}]";
  return output;
}

/**
 * Calculate the Hermite transpose of a 2x2 matrix.
 */
template <typename ValType>
MC2x2Base<ValType> HermTranspose(const MC2x2Base<ValType>& matrix) {
  return MC2x2Base<ValType>(std::conj(matrix.Get(0)), std::conj(matrix.Get(2)),
                            std::conj(matrix.Get(1)), std::conj(matrix.Get(3)));
}

/**
 * Calculate the sum of the diagonal elements.
 */
template <typename ValType>
std::complex<ValType> Trace(const MC2x2Base<ValType>& matrix) {
  return matrix.Get(0) + matrix.Get(3);
}

/**
 * Calculate the Frobenius norm of a matrix. This
 * is the sum of squares over all the real and imaginary values
 * in the matrix.
 */
template <typename ValType>
ValType Norm(const MC2x2Base<ValType>& matrix) {
  return std::norm(matrix.Get(0)) + std::norm(matrix.Get(1)) +
         std::norm(matrix.Get(2)) + std::norm(matrix.Get(3));
}

/**
 * Calculate the L1 norm of a matrix: the sum of absolute values.
 */
template <typename ValType>
ValType SumOfAbsolute(const MC2x2Base<ValType>& matrix) {
  return std::abs(matrix.Get(0)) + std::abs(matrix.Get(1)) +
         std::abs(matrix.Get(2)) + std::abs(matrix.Get(3));
}

/**
 * Element-wise product of two 2x2 matrices.
 */
template <typename ValType>
MC2x2Base<ValType> ElementProduct(const MC2x2Base<ValType>& lhs,
                                  const MC2x2Base<ValType>& rhs) {
  return MC2x2Base<ValType>(lhs.Get(0) * rhs.Get(0), lhs.Get(1) * rhs.Get(1),
                            lhs.Get(2) * rhs.Get(2), lhs.Get(3) * rhs.Get(3));
}

}  // namespace aocommon::scalar

#include "../matrix2x2diag.h"

template <typename ValType>
constexpr aocommon::scalar::MC2x2Base<ValType>::MC2x2Base(
    const aocommon::scalar::MC2x2DiagBase<ValType>& diag)
    : values_{diag.Get(0).real(), diag.Get(0).imag(), 0.0, 0.0, 0.0, 0.0,
              diag.Get(1).real(), diag.Get(1).imag()} {}

#endif
