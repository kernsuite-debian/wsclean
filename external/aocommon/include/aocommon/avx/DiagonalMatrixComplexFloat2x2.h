// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
#define AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H

#include "AvxMacros.h"
#include "VectorComplexFloat2.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>

namespace aocommon::avx {

/**
 * Implements a Diagonal 2x2 Matrix with complex float values.
 * The matrix is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is based on @ref aocommon::MC2x2FDiag but uses AVX-128
 * instructions.
 */
class DiagonalMatrixComplexFloat2x2 {
 public:
  AVX_TARGET DiagonalMatrixComplexFloat2x2() noexcept = default;

  AVX_TARGET /* implicit */
  DiagonalMatrixComplexFloat2x2(VectorComplexFloat2 data) noexcept
      : data_{data} {}

  AVX_TARGET explicit DiagonalMatrixComplexFloat2x2(
      const std::complex<float> a, const std::complex<float> b) noexcept
      : data_{a, b} {}

  AVX_TARGET explicit DiagonalMatrixComplexFloat2x2(
      const std::complex<float> matrix[2]) noexcept
      : data_{VectorComplexFloat2{std::addressof(matrix[0])}} {}

  AVX_TARGET const std::complex<float> Get(size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return data_.Get(index);
  }

  AVX_TARGET void Set(size_t index, std::complex<float> value) noexcept {
    assert(index < 2 && "Index out of bounds.");
    data_.Set(index, value);
  }

  AVX_TARGET explicit operator __m128() const noexcept {
    return static_cast<__m128>(data_);
  }

  template <typename T>
  AVX_TARGET friend DiagonalMatrixComplexFloat2x2 operator*(
      DiagonalMatrixComplexFloat2x2 lhs, T rhs) noexcept {
    return lhs.data_ * rhs;
  }

  AVX_TARGET friend DiagonalMatrixComplexFloat2x2 operator*(
      DiagonalMatrixComplexFloat2x2 lhs,
      DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ * rhs.data_;
  }

  AVX_TARGET DiagonalMatrixComplexFloat2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  AVX_TARGET DiagonalMatrixComplexFloat2x2 HermTranspose() const noexcept {
    // The transpose has no effect for a diagonal matrix.
    return Conjugate();
  }

  AVX_TARGET static DiagonalMatrixComplexFloat2x2 Zero() noexcept {
    return DiagonalMatrixComplexFloat2x2{VectorComplexFloat2::Zero()};
  }

  AVX_TARGET DiagonalMatrixComplexFloat2x2
  operator+=(DiagonalMatrixComplexFloat2x2 value) noexcept {
    return data_ += value.data_;
  }

  AVX_TARGET DiagonalMatrixComplexFloat2x2
  operator-=(DiagonalMatrixComplexFloat2x2 value) noexcept {
    return data_ -= value.data_;
  }

  AVX_TARGET friend bool operator==(
      DiagonalMatrixComplexFloat2x2 lhs,
      DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  AVX_TARGET friend std::ostream& operator<<(
      std::ostream& output, DiagonalMatrixComplexFloat2x2 value) {
    output << "[{" << value.Get(0) << ", " << std::complex<float>{} << "}, {"
           << std::complex<float>{} << ", " << value.Get(1) << "}]";
    return output;
  }

 private:
  VectorComplexFloat2 data_;
};

AVX_TARGET inline DiagonalMatrixComplexFloat2x2 HermTranspose(
    DiagonalMatrixComplexFloat2x2 matrix) noexcept {
  return matrix.HermTranspose();
}

/// @returns the sum of the diagonal elements.
AVX_TARGET inline std::complex<float> Trace(
    DiagonalMatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Get(0) + matrix.Get(1);
}

/// @returns the Frobenius norm of the matrix.
AVX_TARGET inline float Norm(DiagonalMatrixComplexFloat2x2 matrix) noexcept {
  return std::norm(matrix.Get(0)) + std::norm(matrix.Get(1));
}

/// Returns the original diagonal matrix. Can be useful for templated code so
/// that it can be used for both the diagonal and full matrix variants.
AVX_TARGET inline DiagonalMatrixComplexFloat2x2 Diagonal(
    DiagonalMatrixComplexFloat2x2 matrix) noexcept {
  return matrix;
}

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
