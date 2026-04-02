// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX_VECTOR_FLOAT_4_H_
#define AOCOMMON_AVX_VECTOR_FLOAT_4_H_

#include "AvxMacros.h"

#include <cassert>
#include <immintrin.h>
#include <memory>

namespace aocommon::avx {

/**
 * Implements a Vector with 4 float values.
 * The vector is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexFloat2, but can be used by itself.
 *
 * @note since 4 floats use 128-bits this vector uses AVX-128 vectors, despite
 * the 256 name.
 */
class VectorFloat4 {
 public:
  AVX_TARGET explicit VectorFloat4() noexcept = default;

  AVX_TARGET /* implicit */ VectorFloat4(__m128 data) noexcept : data_{data} {}

  AVX_TARGET explicit VectorFloat4(float value) noexcept
      : data_{_mm_set1_ps(value)} {}

  AVX_TARGET explicit VectorFloat4(float a, float b, float c, float d) noexcept
      : data_{_mm_setr_ps(a, b, c, d)} {}

  AVX_TARGET explicit VectorFloat4(const float vector[4]) noexcept
      : data_{_mm_loadu_ps(vector)} {}

  AVX_TARGET float operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  AVX_TARGET __m128 Value() const noexcept { return data_; }

  AVX_TARGET static VectorFloat4 Zero() noexcept {
    return VectorFloat4{_mm_setzero_ps()};
  }

  AVX_TARGET VectorFloat4& operator+=(VectorFloat4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  AVX_TARGET VectorFloat4& operator-=(VectorFloat4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  AVX_TARGET friend VectorFloat4 operator+(VectorFloat4 lhs,
                                           VectorFloat4 rhs) noexcept {
    return lhs += rhs;
  }

  AVX_TARGET friend VectorFloat4 operator-(VectorFloat4 lhs,
                                           VectorFloat4 rhs) noexcept {
    return lhs -= rhs;
  }

  AVX_TARGET friend VectorFloat4 operator*(VectorFloat4 lhs,
                                           VectorFloat4 rhs) noexcept {
    return _mm_mul_ps(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend VectorFloat4 operator^(VectorFloat4 lhs,
                                           VectorFloat4 rhs) noexcept {
    return _mm_xor_ps(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend bool operator==(VectorFloat4 lhs,
                                    VectorFloat4 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3];
  }

 private:
  __m128 data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_VECTOR_FLOAT_4_H
