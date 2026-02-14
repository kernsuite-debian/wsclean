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

  AVX_TARGET std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return data_[index];
  }

  AVX_TARGET explicit operator __m128() const noexcept {
    return static_cast<__m128>(data_);
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

  AVX_TARGET friend bool operator==(
      DiagonalMatrixComplexFloat2x2 lhs,
      DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  AVX_TARGET friend std::ostream& operator<<(
      std::ostream& output, DiagonalMatrixComplexFloat2x2 value) {
    output << "[{" << value[0] << ", " << std::complex<float>{} << "}, {"
           << std::complex<float>{} << ", " << value[1] << "}]";
    return output;
  }

 private:
  VectorComplexFloat2 data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_DIAGONAL_MATRIX_COMPLEX_FLOAT_2X2_H
