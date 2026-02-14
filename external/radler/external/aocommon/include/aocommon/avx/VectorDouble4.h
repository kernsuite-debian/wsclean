// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
#define AOCOMMON_AVX256_VECTOR_DOUBLE_4_H

#include <cassert>
#include <immintrin.h>
#include <memory>

#include "AvxMacros.h"

namespace aocommon::avx {

/**
 * Implements a Vector with 4 double values.
 * The vector is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::VectorComplexDouble2, but can be used by itself.
 */
class VectorDouble4 {
 public:
  AVX_TARGET explicit VectorDouble4() noexcept = default;

  AVX_TARGET /* implicit */ VectorDouble4(__m256d data) noexcept
      : data_{data} {}

  AVX_TARGET explicit VectorDouble4(double value) noexcept
      : data_{_mm256_set1_pd(value)} {}

  AVX_TARGET explicit VectorDouble4(double a, double b, double c,
                                    double d) noexcept
      : data_{_mm256_setr_pd(a, b, c, d)} {}

  AVX_TARGET explicit VectorDouble4(const double vector[4]) noexcept
      : data_{_mm256_loadu_pd(vector)} {}

  AVX_TARGET explicit VectorDouble4(const float vector[4]) noexcept
      : data_(_mm256_cvtps_pd(_mm_loadu_ps(std::addressof(vector[0])))) {}

  AVX_TARGET VectorDouble4(const VectorDouble4&) = default;
  AVX_TARGET VectorDouble4& operator=(const VectorDouble4&) = default;

  AVX_TARGET double operator[](size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_[index];
  }

  AVX_TARGET __m256d Value() const noexcept { return data_; }

  /** Assign data stored by 4 element vector to destination buffer. */
  AVX_TARGET void AssignTo(double* destination) const noexcept {
    _mm256_storeu_pd(destination, data_);
  }

  AVX_TARGET static VectorDouble4 Zero() noexcept {
    return VectorDouble4{_mm256_setzero_pd()};
  }

  AVX_TARGET VectorDouble4& operator+=(VectorDouble4 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  AVX_TARGET VectorDouble4& operator-=(VectorDouble4 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  AVX_TARGET friend VectorDouble4 operator+(VectorDouble4 lhs,
                                            VectorDouble4 rhs) noexcept {
    return VectorDouble4(lhs.data_ + rhs.data_);
  }

  AVX_TARGET friend VectorDouble4 operator-(VectorDouble4 lhs,
                                            VectorDouble4 rhs) noexcept {
    return lhs -= rhs;
  }

  AVX_TARGET friend VectorDouble4 operator*(VectorDouble4 lhs,
                                            VectorDouble4 rhs) noexcept {
    return _mm256_mul_pd(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend VectorDouble4 operator^(VectorDouble4 lhs,
                                            VectorDouble4 rhs) noexcept {
    return _mm256_xor_pd(lhs.data_, rhs.data_);
  }

  AVX_TARGET friend bool operator==(VectorDouble4 lhs,
                                    VectorDouble4 rhs) noexcept {
    return lhs.data_[0] == rhs.data_[0] && lhs.data_[1] == rhs.data_[1] &&
           lhs.data_[2] == rhs.data_[2] && lhs.data_[3] == rhs.data_[3];
  }

 private:
  __m256d data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_VECTOR_DOUBLE_4_H
