// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
#define AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H

#include "VectorComplexDouble2.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>
#include <iostream>  // DEBUG

namespace aocommon::avx {

/**
 * Implements a Diagonal 2x2 Matrix with complex double values.
 * The matrix is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is based on @ref aocommon::MC2x2Diag but uses AVX-256
 * instructions.
 */
class DiagonalMatrixComplexDouble2x2 {
 public:
  AVX_TARGET DiagonalMatrixComplexDouble2x2() noexcept = default;

  AVX_TARGET DiagonalMatrixComplexDouble2x2(
      const DiagonalMatrixComplexDouble2x2&) noexcept = default;

  AVX_TARGET /* implicit */
  DiagonalMatrixComplexDouble2x2(VectorComplexDouble2 data) noexcept
      : data_{data} {}

  /**
   * Construct from (length 2) data buffer
   */
  AVX_TARGET explicit DiagonalMatrixComplexDouble2x2(
      const std::complex<double> matrix[2]) noexcept
      : data_{VectorComplexDouble2(&matrix[0])} {}

  AVX_TARGET explicit DiagonalMatrixComplexDouble2x2(
      const double* data) noexcept
      : data_(data[0], data[1]) {}

  /**
   * Construct from initializer list, values are internally converted
   * to complex type. Assumes that list has size two.
   */
  template <typename ValType>
  AVX_TARGET DiagonalMatrixComplexDouble2x2(std::initializer_list<ValType> data)
      : data_(*data.begin(), *(data.begin() + 1)) {
    assert(data.size() == 2);
  }

  AVX_TARGET DiagonalMatrixComplexDouble2x2(std::complex<double> a,
                                            std::complex<double> b) noexcept
      : data_{a, b} {}

  AVX_TARGET VectorComplexDouble2 Data() const noexcept { return data_; }

  AVX_TARGET std::complex<double> operator[](size_t index) const noexcept {
    assert(index < 2);
    return data_[index];
  }

  AVX_TARGET DiagonalMatrixComplexDouble2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  AVX_TARGET DiagonalMatrixComplexDouble2x2 HermTranspose() const noexcept {
    // The transpose has no effect for a diagonal matrix.
    return Conjugate();
  }

  AVX_TARGET DiagonalMatrixComplexDouble2x2& operator=(
      const DiagonalMatrixComplexDouble2x2&) = default;

  AVX_TARGET DiagonalMatrixComplexDouble2x2& operator+=(
      DiagonalMatrixComplexDouble2x2 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  template <typename T>
  AVX_TARGET DiagonalMatrixComplexDouble2x2& operator*=(T value) noexcept {
    // TODO could use avx
    *this = DiagonalMatrixComplexDouble2x2(data_[0] * value, data_[1] * value);
    return *this;
  }

  AVX_TARGET DiagonalMatrixComplexDouble2x2& operator*=(
      DiagonalMatrixComplexDouble2x2 value) noexcept {
    *this = *this * value;
    return *this;
  }

  template <typename T>
  AVX_TARGET DiagonalMatrixComplexDouble2x2& operator/=(T value) noexcept {
    // TODO could use avx
    *this = DiagonalMatrixComplexDouble2x2(data_[0] / value, data_[1] / value);
    return *this;
  }

  AVX_TARGET friend bool operator==(
      DiagonalMatrixComplexDouble2x2 lhs,
      DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  AVX_TARGET friend std::ostream& operator<<(
      std::ostream& output, DiagonalMatrixComplexDouble2x2 value) {
    output << "[{" << value[0] << ", " << std::complex<double>{} << "}, {"
           << std::complex<double>{} << ", " << value[1] << "}]";
    return output;
  }

  AVX_TARGET static DiagonalMatrixComplexDouble2x2 Zero() noexcept {
    return DiagonalMatrixComplexDouble2x2(VectorComplexDouble2::Zero());
  }

  AVX_TARGET static DiagonalMatrixComplexDouble2x2 Unity() {
    return DiagonalMatrixComplexDouble2x2(1.0, 1.0);
  }

  AVX_TARGET friend DiagonalMatrixComplexDouble2x2 operator+(
      DiagonalMatrixComplexDouble2x2 lhs,
      DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return lhs += rhs;
  }

  AVX_TARGET friend DiagonalMatrixComplexDouble2x2 operator*(
      DiagonalMatrixComplexDouble2x2 lhs,
      DiagonalMatrixComplexDouble2x2 rhs) noexcept {
    return lhs.data_ * rhs.data_;
  }

  template <typename T>
  AVX_TARGET friend DiagonalMatrixComplexDouble2x2 operator*(
      DiagonalMatrixComplexDouble2x2 lhs, T rhs) noexcept {
    return VectorComplexDouble2(lhs.data_[0] * rhs, lhs.data_[1] * rhs);
  }

  template <typename T>
  AVX_TARGET friend DiagonalMatrixComplexDouble2x2 operator*(
      DiagonalMatrixComplexDouble2x2 lhs, const T* rhs) noexcept {
    return lhs * DiagonalMatrixComplexDouble2x2(rhs);
  }

  template <typename T>
  AVX_TARGET friend DiagonalMatrixComplexDouble2x2 operator/(
      DiagonalMatrixComplexDouble2x2 lhs, T rhs) noexcept {
    return VectorComplexDouble2(lhs.data_[0] / rhs, lhs.data_[1] / rhs);
  }

 private:
  VectorComplexDouble2 data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_DOUBLE_2X2_H
