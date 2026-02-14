// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX_VECTOR_COMPLEX_FLOAT_2_H_
#define AOCOMMON_AVX_VECTOR_COMPLEX_FLOAT_2_H_

#include "VectorFloat4.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>

namespace aocommon::avx {

/**
 * Implements a Vector with 2 complex float values.
 * The vector is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is an implementation detail of
 * @ref aocommon::Avx256::MatrixComplexFloat2x2, but can be used by itself.
 */
class VectorComplexFloat2 {
 public:
  AVX_TARGET VectorComplexFloat2() noexcept = default;

  AVX_TARGET /* implicit */ VectorComplexFloat2(VectorFloat4 data) noexcept
      : data_{data} {}

  AVX_TARGET explicit VectorComplexFloat2(std::complex<float> a,
                                          std::complex<float> b) noexcept
      : data_{VectorFloat4{a.real(), a.imag(), b.real(), b.imag()}} {}

  AVX_TARGET explicit VectorComplexFloat2(
      const std::complex<float> vector[2]) noexcept
      // reinterpret_cast explicitly allowed per [complex.numbers.general]/4.
      // (http://www.eelis.net/c++draft/complex.numbers#general-4)
      : data_{VectorFloat4{
            reinterpret_cast<const float*>(std::addressof(vector[0]))}} {}

  AVX_TARGET std::complex<float> operator[](size_t index) const noexcept {
    assert(index < 2 && "Index out of bounds.");
    return {data_[2 * index], data_[2 * index + 1]};
  }

  AVX_TARGET explicit operator __m128() const noexcept { return data_.Value(); }

  AVX_TARGET VectorComplexFloat2 Conjugate() const noexcept {
    // Xor-ing a float with  0.0 will not change the value.
    // Xor-ing a float with -0.0 will change the sign of the value.
    __m128 mask = _mm_setr_ps(0.0, -0.0, 0.0, -0.0);
    return data_ ^ mask;
  }

  AVX_TARGET static VectorComplexFloat2 Zero() noexcept {
    return VectorComplexFloat2{VectorFloat4::Zero()};
  }

  AVX_TARGET VectorComplexFloat2& operator+=(
      VectorComplexFloat2 value) noexcept {
    data_ += value.data_;
    return *this;
  }

  AVX_TARGET VectorComplexFloat2& operator-=(
      VectorComplexFloat2 value) noexcept {
    data_ -= value.data_;
    return *this;
  }

  AVX_TARGET friend VectorComplexFloat2 operator+(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
    return lhs += rhs;
  }

  AVX_TARGET friend VectorComplexFloat2 operator-(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
    return lhs -= rhs;
  }

  /// Multiplies the elements of 2 vectors in parallel.
  ///
  /// r[0] = lhs[0] * rhs[0]
  /// r[1] = lhs[1] * rhs[1]
  AVX_TARGET friend VectorComplexFloat2 operator*(
      VectorComplexFloat2 lhs, VectorComplexFloat2 rhs) noexcept {
    // The complex multiplication L * R is performed by:
    // L * R = (L Re * R Re) - (L Im * R Im) +
    //        ((L Re * R Im) + (L Im * R Re)) i

    // Transformed for AVX this becomes
    //
    // mul 1:
    //
    // Lv1 | L Re        | L Re        |
    // Lv2 | R Re        | R Im        | Note this is rhs
    //     ------------- * ----------- *
    // Lv3 | L Re * R Re | L Re * R Im |
    //
    // mul 2:
    //
    // Rv1 | L Im        | L Im        |
    // Rv2 | R Im        | R Re        |
    //     ------------- * ----------- *
    // Rv3 | L Im * R Im | L Im * R Re |
    //
    // add sub
    // Lv3 | L Re * R Re               | L Re * R Im               |
    // Rv3 | L Im * R Im               | L Im * R Re               |
    //     --------------------------- - ------------------------- +
    //     | L Re * R Re - L Im * R Im | L Re * R Im + L Im * R Re |
    //
    // It's also possible to do an fmul add sub
    // Which does (Lv1 fmul Lv2) add/sub Rv3

    // The algorithm "uses" 256 bit vectors in two 128 bit vectors. Testing
    // with 1 256 bit vector didn't improve much. The main issue is:
    //   128 bit: mul + fmaddsub
    //   256 bit: mul + addsub
    // The performance of fmaddsub and addsub seem to be the same according to
    // Intel. The additional register manipulation actually makes it more
    // expense.
    //
    // Not there can be an improvement by doing two diagonal matrices in
    // parallel. However that requires more code changes.

    __m128 Rv1 =
        _mm_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b11'11'01'01);
    __m128 Rv2 =
        _mm_shuffle_ps(rhs.data_.Value(), rhs.data_.Value(), 0b10'11'00'01);
    __m128 Lv1 =
        _mm_shuffle_ps(lhs.data_.Value(), lhs.data_.Value(), 0b10'10'00'00);
    __m128 Rv3 = _mm_mul_ps(Rv1, Rv2);
    return VectorFloat4{_mm_fmaddsub_ps(Lv1, rhs.data_.Value(), Rv3)};
  }

  AVX_TARGET friend VectorComplexFloat2 operator*(
      VectorComplexFloat2 lhs, std::complex<float> rhs) noexcept {
    return lhs * VectorComplexFloat2{rhs, rhs};
  }

  AVX_TARGET friend VectorComplexFloat2 operator*(
      std::complex<float> lhs, VectorComplexFloat2 rhs) noexcept {
    return rhs * lhs;
  }

  AVX_TARGET friend bool operator==(VectorComplexFloat2 lhs,
                                    VectorComplexFloat2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  AVX_TARGET friend std::ostream& operator<<(std::ostream& output,
                                             VectorComplexFloat2 value) {
    output << '[' << value[0] << ", " << value[1] << ']';
    return output;
  }

 private:
  VectorFloat4 data_;
};

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_VECTOR_COMPLEX_FLOAT_2_H
