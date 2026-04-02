// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX_VECTOR_FLOAT_8_H_
#define AOCOMMON_AVX_VECTOR_FLOAT_8_H_

#include "AvxMacros.h"

#include <cassert>
#include <immintrin.h>
#include <memory>

namespace aocommon::avx {

/**
 * Implements a Vector with 8 float values.
 * The vector is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexFloat4, but can be used by itself.
 */
class VectorFloat8 {
 public:
  AVX_TARGET explicit VectorFloat8() noexcept = default;

  AVX_TARGET /* implicit */ VectorFloat8(__m256 data) noexcept : data_{data} {}

  AVX_TARGET explicit VectorFloat8(float value) noexcept
      : data_{_mm256_set1_ps(value)} {}

  AVX_TARGET explicit VectorFloat8(float a, float b, float c, float d, float e,
                                   float f, float g, float h) noexcept
      : data_{_mm256_setr_ps(a, b, c, d, e, f, g, h)} {}

  AVX_TARGET explicit VectorFloat8(const float vector[8]) noexcept
      : data_{_mm256_loadu_ps(vector)} {}

  AVX_TARGET explicit VectorFloat8(const double vector[8]) noexcept {
    __m256 lo = _mm256_castps128_ps256(
        _mm256_cvtpd_ps(_mm256_loadu_pd(std::addressof(vector[0]))));
    __m128 hi = _mm256_cvtpd_ps(_mm256_loadu_pd(std::addressof(vector[4])));

    data_ = _mm256_insertf128_ps(lo, hi, 1);
  }

  AVX_TARGET float operator[](size_t index) const noexcept {
    assert(index < 8 && "Index out of bounds.");
    return data_[index];
  }

  AVX_TARGET __m256 Value() const noexcept { return data_; }

  /** Assign data stored by 8 element vector to destination buffer */
  AVX_TARGET void AssignTo(float* destination) const noexcept {
    _mm256_storeu_ps(destination, data_);
  }

  AVX_TARGET void AssignTo(double* destination) const noexcept {
    _mm256_storeu_pd(destination,
                     _mm256_cvtps_pd(_mm256_castps256_ps128(data_)));
    _mm256_storeu_pd(destination + 4,
                     _mm256_cvtps_pd(_mm256_extractf128_ps(data_, 1)));
  }

  AVX_TARGET static VectorFloat8 Zero() noexcept {
    return VectorFloat8{_mm256_setzero_ps()};
  }

  AVX_TARGET VectorFloat8& operator+=(VectorFloat8 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  AVX_TARGET VectorFloat8& operator-=(VectorFloat8 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  AVX_TARGET friend VectorFloat8 operator+(VectorFloat8 lhs,
                                           VectorFloat8 rhs) noexcept {
    return lhs += rhs;
  }

  AVX_TARGET friend VectorFloat8 operator-(VectorFloat8 lhs,
                                           VectorFloat8 rhs) noexcept {
    return lhs -= rhs;
  }

  AVX_TARGET friend VectorFloat8 operator*(VectorFloat8 lhs,
                                           VectorFloat8 rhs) noexcept {
    return _mm256_mul_ps(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend VectorFloat8 operator^(VectorFloat8 lhs,
                                           VectorFloat8 rhs) noexcept {
    return _mm256_xor_ps(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend bool operator==(VectorFloat8 lhs,
                                    VectorFloat8 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3] &&
           lhs.data_[4] == rhs.data_[4] && lhs.data_[5] == rhs.data_[5] &&
           lhs.data_[6] == rhs.data_[6] && lhs.data_[7] == rhs.data_[7];
  }

 private:
  __m256 data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_VECTOR_FLOAT_8_H
