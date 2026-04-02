#ifndef AOCOMMON_SCALAR_MATRIX_2X2_DIAG_H_
#define AOCOMMON_SCALAR_MATRIX_2X2_DIAG_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>
#include <sstream>

#include "matrix2x2.h"

namespace aocommon::scalar {

/**
 * Class implements a 2x2 complex-valued diagonal matrix.
 */
template <typename ValType>
class MC2x2DiagBase {
 public:
  MC2x2DiagBase() {}

  MC2x2DiagBase(const MC2x2DiagBase<ValType>& source) noexcept {
    _values[0] = source._values[0];
    _values[1] = source._values[1];
  }

  template <typename OtherValType>
  explicit MC2x2DiagBase(const MC2x2DiagBase<OtherValType>& source) noexcept {
    _values[0] = static_cast<std::complex<OtherValType>>(source.Get(0));
    _values[1] = static_cast<std::complex<OtherValType>>(source.Get(1));
  }

  /**
   * Construct MC2x2Base object from (length 2) data buffer
   */
  template <typename T>
  explicit MC2x2DiagBase(const T* source) {
    _values[0] = source[0];
    _values[1] = source[1];
  }

  /**
   * Construct object from diagonal input values. Values are
   * converted to complex type.
   */
  MC2x2DiagBase(ValType m00, ValType m11) {
    _values[0] = m00;
    _values[1] = m11;
  }

  /**
   * Construct MC2x2Base object from two complex-valued input values.
   */
  MC2x2DiagBase(std::complex<ValType> m00, std::complex<ValType> m11) {
    _values[0] = m00;
    _values[1] = m11;
  }

  /**
   * Construct via initializer list, will be converted to complex type.
   * Assumes that list has size two.
   */
  MC2x2DiagBase(std::initializer_list<ValType> list) {
    assert(list.size() == 2);
    std::copy_n(list.begin(), 2, &_values[0]);
  }

  /**
   * Construct via initializer list. Assumes that list has size two.
   */
  MC2x2DiagBase(std::initializer_list<std::complex<ValType>> list) {
    assert(list.size() == 2);
    std::copy_n(list.begin(), 2, &_values[0]);
  }

  MC2x2DiagBase<ValType>& operator=(const MC2x2DiagBase<ValType>& source) {
    _values[0] = source._values[0];
    _values[1] = source._values[1];
    return *this;
  }

  bool operator==(const MC2x2DiagBase<ValType>& rhs) const {
    return _values[0] == rhs._values[0] && _values[1] == rhs._values[1];
  }

  MC2x2DiagBase<ValType> operator*(const MC2x2DiagBase<ValType>& source) const {
    return Multiply(source);
  }

  MC2x2DiagBase<ValType>& operator*=(const MC2x2DiagBase<ValType>& source) {
    (*this) = (*this) * source;
    return *this;
  }

  /**
   * Scalar multiplication operator, real valued rhs
   */
  MC2x2DiagBase<ValType> operator*(ValType rhs) const {
    MC2x2DiagBase<ValType> dest(*this);
    MC2x2DiagBase::ScalarMultiply(dest._values, rhs);
    return dest;
  }

  /**
   * Scalar multiplication assignment operator, real valued rhs
   */
  MC2x2DiagBase<ValType>& operator*=(ValType rhs) {
    MC2x2DiagBase::ScalarMultiply(_values, rhs);
    return *this;
  }

  /**
   * Scalar multiplication operator, complex valued rhs
   */
  MC2x2DiagBase<ValType> operator*(std::complex<ValType> rhs) const {
    MC2x2DiagBase<ValType> dest(*this);
    MC2x2DiagBase::ScalarMultiply(dest._values, rhs);
    return dest;
  }

  /**
   * Scalar multiplication assignment operator, complex valued rhs
   */
  MC2x2DiagBase<ValType>& operator*=(std::complex<ValType> rhs) {
    MC2x2DiagBase::ScalarMultiply(_values, rhs);
    return *this;
  }

  /**
   * Scalar division operator
   */
  MC2x2DiagBase<ValType> operator/(ValType rhs) const {
    MC2x2DiagBase<ValType> dest(*this);
    MC2x2DiagBase::ScalarMultiply(dest._values, 1.0 / rhs);
    return dest;
  }

  /**
   * Scalar division assignment operator
   */
  MC2x2DiagBase<ValType>& operator/=(ValType rhs) {
    MC2x2DiagBase::ScalarMultiply(_values, 1.0 / rhs);
    return *this;
  }

  MC2x2DiagBase<ValType> operator+(const MC2x2DiagBase<ValType>& source) {
    return MC2x2DiagBase<ValType>(_values[0] + source._values[0],
                                  _values[1] + source._values[1]);
  }

  MC2x2DiagBase<ValType>& operator+=(const MC2x2DiagBase<ValType>& source) {
    _values[0] += source._values[0];
    _values[1] += source._values[1];
    return *this;
  }

  MC2x2DiagBase<ValType>& operator-=(const MC2x2DiagBase<ValType>& source) {
    _values[0] -= source._values[0];
    _values[1] -= source._values[1];
    return *this;
  }

  std::complex<ValType> Get(size_t index) const { return _values[index]; }
  void Set(size_t index, std::complex<ValType> value) {
    _values[index] = value;
  }

  /**
   * Return MC2x2Base matrix filled with zeros
   */
  static MC2x2DiagBase<ValType> Zero() {
    return MC2x2DiagBase<ValType>(0.0, 0.0);
  }

  /**
   * Return 2x2 identity matrix
   */
  static MC2x2DiagBase<ValType> Unity() { return MC2x2DiagBase(1.0, 1.0); }

  /**
   * Return 2x2 matrix filled with NaN values
   */
  static MC2x2DiagBase<ValType> NaN() {
    return MC2x2DiagBase<ValType>(
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()),
        std::complex<ValType>(std::numeric_limits<ValType>::quiet_NaN(),
                              std::numeric_limits<ValType>::quiet_NaN()));
  }

  /**
   * Get pointer to underlying data
   */
  std::complex<ValType>* Data() { return _values; }
  const std::complex<ValType>* Data() const { return _values; }

  /**
   * Matrix multiplication, alias for the overloaded * operator
   */
  MC2x2DiagBase<ValType> Multiply(const MC2x2DiagBase<ValType>& rhs) const {
    return MC2x2DiagBase<ValType>(_values[0] * rhs.Get(0),
                                  _values[1] * rhs.Get(1));
  }

  /**
   * Scalar multiplication of diagonal matrix.
   */
  template <typename T>
  static void ScalarMultiply(std::complex<T>* dest, T factor) {
    for (size_t p = 0; p != 2; ++p) dest[p] *= factor;
  }

  /**
   * Scalar multiplication of matrix.
   */
  template <typename T>
  static void ScalarMultiply(T* dest, T factor) {
    for (size_t p = 0; p != 2; ++p) dest[p] *= factor;
  }

  /**
   * Compute Hermitian transpose of matrix
   */
  MC2x2DiagBase<ValType> HermTranspose() const {
    return MC2x2DiagBase(std::conj(_values[0]), std::conj(_values[1]));
  }

  /**
   * Compute conjugate of matrix (for a Diagonal matrix, this is the same as
   * HermTranspose()).
   */
  MC2x2DiagBase<ValType> Conjugate() const { return HermTranspose(); }

  /**
   * Computation is short-circuited to avoid unnecessary computation of
   * imaginary parts
   * @returns real(A^H * A)
   */
  std::array<ValType, 2> RealHermitianSquare() const {
    const MC2x2DiagBase<ValType> conjugate = Conjugate();
    return {(std::real(conjugate.Get(0)) * std::real(Get(0))) -
                (std::imag(conjugate.Get(0)) * std::imag(Get(0))),
            (std::real(conjugate.Get(1)) * std::real(Get(1))) -
                (std::imag(conjugate.Get(1)) * std::imag(Get(1)))};
  }

 private:
  std::complex<ValType> _values[2];
};

/**
 * Diagonal - non-diagonal Matrix multiplication operator
 */
template <typename ValType>
scalar::MC2x2Base<ValType> operator*(const MC2x2DiagBase<ValType>& lhs,
                                     const scalar::MC2x2Base<ValType>& rhs) {
  return scalar::MC2x2Base<ValType>(
      lhs.Get(0) * rhs.Get(0), lhs.Get(0) * rhs.Get(1), lhs.Get(1) * rhs.Get(2),
      lhs.Get(1) * rhs.Get(3));
}

/**
 * Non-diagonal - diagonal Matrix multiplication operator
 */
template <typename ValType>
scalar::MC2x2Base<ValType> operator*(const scalar::MC2x2Base<ValType>& lhs,
                                     const MC2x2DiagBase<ValType>& rhs) {
  return scalar::MC2x2Base<ValType>(
      lhs.Get(0) * rhs.Get(0), lhs.Get(1) * rhs.Get(1), lhs.Get(2) * rhs.Get(0),
      lhs.Get(3) * rhs.Get(1));
}

/**
 * Left shift operator to write the matrix to ostream
 */
template <typename ValType>
std::ostream& operator<<(std::ostream& output,
                         const MC2x2DiagBase<ValType>& value) {
  output << "[{" << value.Get(0) << ", " << std::complex<float>{0.0} << "}, {"
         << std::complex<float>{0.0} << ", " << value.Get(1) << "}]";
  return output;
}

/**
 * Obtain the diagonal of a 2x2 matrix
 */
template <typename ValType>
MC2x2DiagBase<ValType> Diagonal(const scalar::MC2x2Base<ValType>& matrix) {
  return MC2x2DiagBase<ValType>(matrix.Get(0), matrix.Get(3));
}

template <typename ValType>
inline MC2x2DiagBase<ValType> HermTranspose(
    MC2x2DiagBase<ValType> matrix) noexcept {
  return matrix.HermTranspose();
}

/// @returns the sum of the diagonal elements.
template <typename ValType>
inline std::complex<float> Trace(MC2x2DiagBase<ValType> matrix) noexcept {
  return matrix.Get(0) + matrix.Get(1);
}

/// @returns the Frobenius norm of the matrix.
template <typename ValType>
inline float Norm(MC2x2DiagBase<ValType> matrix) noexcept {
  return std::norm(matrix.Get(0)) + std::norm(matrix.Get(1));
}

/// Returns the original diagonal matrix. Can be useful for templated code so
/// that it can be used for both the diagonal and full matrix variants.
template <typename ValType>
inline MC2x2DiagBase<ValType> Diagonal(MC2x2DiagBase<ValType> matrix) noexcept {
  return matrix;
}

}  // namespace aocommon::scalar

#endif
