// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_AVX_MATRIX_COMPLEX_FLOAT_2X2_H_
#define AOCOMMON_AVX_MATRIX_COMPLEX_FLOAT_2X2_H_

#include "AvxMacros.h"
#include "DiagonalMatrixComplexFloat2x2.h"
#include "VectorComplexFloat4.h"

#include <cassert>
#include <complex>
#include <immintrin.h>
#include <ostream>
#include <iostream>  // DEBUG

namespace aocommon::avx {

class MatrixComplexDouble2x2;

/**
 * Implements a 2x2 Matrix with complex float values.
 * The matrix is not initialized by default for performance considerations,
 * a call to Zero() can be used to do that when needed.
 *
 * This class is based on @ref aocommon::MC2x2F but uses AVX-256 instructions.
 */
class MatrixComplexFloat2x2 {
 public:
  AVX_TARGET MatrixComplexFloat2x2() noexcept = default;

  AVX_TARGET /* implicit */ MatrixComplexFloat2x2(
      VectorComplexFloat4 data) noexcept
      : data_{data} {}

  AVX_TARGET explicit MatrixComplexFloat2x2(std::complex<float> a,
                                            std::complex<float> b,
                                            std::complex<float> c,
                                            std::complex<float> d) noexcept
      : data_{a, b, c, d} {}

  AVX_TARGET explicit MatrixComplexFloat2x2(
      const std::complex<float> matrix[4]) noexcept
      : data_(matrix) {}

  AVX_TARGET explicit MatrixComplexFloat2x2(
      const std::complex<double> matrix[4]) noexcept
      : data_(matrix) {}

  // Supplied as a const ref argument implemented in
  // common/avx256/MatrixComplexDouble2x2.h. This avoids circular dependencies
  // in the headers.
  AVX_TARGET explicit MatrixComplexFloat2x2(
      const MatrixComplexDouble2x2& matrix) noexcept;

  AVX_TARGET const std::complex<float> Get(size_t index) const noexcept {
    assert(index < 4 && "Index out of bounds.");
    return data_.Get(index);
  }

  AVX_TARGET void Set(size_t index, std::complex<float> value) noexcept {
    assert(index < 4 && "Index out of bounds.");
    data_.Set(index, value);
  }

  AVX_TARGET MatrixComplexFloat2x2 Conjugate() const noexcept {
    return data_.Conjugate();
  }

  AVX_TARGET MatrixComplexFloat2x2 Transpose() const noexcept {
    // Note the compiler uses intrinsics without assistance.
    return MatrixComplexFloat2x2{data_.Get(0), data_.Get(2), data_.Get(1),
                                 data_.Get(3)};
  }

  AVX_TARGET explicit operator __m256() const noexcept {
    return static_cast<__m256>(data_);
  }

  /** Returns the sum of the diagonal elements. */
  AVX_TARGET std::complex<float> Trace() const noexcept {
    // Trace = M[0] + M[3]

    __m256 ret = static_cast<__m256>(data_);
    ret +=
        // Moves M[3] to the location of M[0] and adds it to ret.
        _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(ret), 0b11));
    return {ret[0], ret[1]};
  }

  /** Assign data stored by 2x2 complex matrix to destination buffer */
  AVX_TARGET void AssignTo(std::complex<float>* destination) const noexcept {
    data_.AssignTo(destination);
  }

  AVX_TARGET void AssignTo(std::complex<double>* destination) const noexcept {
    data_.AssignTo(destination);
  }

  AVX_TARGET static MatrixComplexFloat2x2 Zero() noexcept {
    return MatrixComplexFloat2x2{VectorComplexFloat4::Zero()};
  }

  AVX_TARGET static MatrixComplexFloat2x2 Unity() noexcept {
    return MatrixComplexFloat2x2{
        std::complex<float>(1.0f, 0.0f), std::complex<float>(0.0f, 0.0f),
        std::complex<float>(0.0f, 0.0f), std::complex<float>(1.0f, 0.0f)};
  }

  AVX_TARGET static MatrixComplexFloat2x2 NaN() noexcept {
    return MatrixComplexFloat2x2{
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()},
        std::complex<float>{std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN()}};
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator+(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return lhs += rhs;
  }

  AVX_TARGET MatrixComplexFloat2x2
  operator+(const DiagonalMatrixComplexFloat2x2& rhs) const {
    return MatrixComplexFloat2x2(Get(0) + rhs.Get(0), Get(1), Get(2),
                                 Get(3) + rhs.Get(1));
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator+(
      const DiagonalMatrixComplexFloat2x2& lhs,
      const MatrixComplexFloat2x2& rhs) {
    return MatrixComplexFloat2x2(lhs.Get(0) + rhs.Get(0), rhs.Get(1),
                                 rhs.Get(2), lhs.Get(1) + rhs.Get(3));
  }

  AVX_TARGET MatrixComplexFloat2x2& operator+=(
      MatrixComplexFloat2x2 rhs) noexcept {
    data_ += rhs.data_;
    return *this;
  }

  AVX_TARGET MatrixComplexFloat2x2& operator+=(
      const DiagonalMatrixComplexFloat2x2& rhs) {
    *this = *this + rhs;
    return *this;
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator-(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return lhs -= rhs;
  }

  AVX_TARGET MatrixComplexFloat2x2
  operator-(const DiagonalMatrixComplexFloat2x2& rhs) const {
    return MatrixComplexFloat2x2(Get(0) - rhs.Get(0), Get(1), Get(2),
                                 Get(3) - rhs.Get(1));
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator-(
      const DiagonalMatrixComplexFloat2x2& lhs,
      const MatrixComplexFloat2x2& rhs) {
    return MatrixComplexFloat2x2(lhs.Get(0) - rhs.Get(0), -rhs.Get(1),
                                 -rhs.Get(2), lhs.Get(1) - rhs.Get(3));
  }

  AVX_TARGET MatrixComplexFloat2x2& operator-=(
      MatrixComplexFloat2x2 rhs) noexcept {
    data_ -= rhs.data_;
    return *this;
  }

  AVX_TARGET MatrixComplexFloat2x2& operator-=(
      const DiagonalMatrixComplexFloat2x2& rhs) {
    *this = *this - rhs;
    return *this;
  }

  AVX_TARGET MatrixComplexFloat2x2& operator*=(
      MatrixComplexFloat2x2 rhs) noexcept {
    data_ = data_ * rhs.data_;
    return *this;
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    // The 2x2 matrix multiplication is done using the following algorithm.
    // ret.a = lhs.a * rhs.a + lhs.b * rhs.c
    // ret.b = lhs.a * rhs.b + lhs.b * rhs.d
    // ret.c = lhs.c * rhs.a + lhs.d * rhs.c
    // ret.d = lhs.c * rhs.b + lhs.d * rhs.d
    //       | c1    | c2    | c3    | c4    |
    //       | s1            | s2            |
    //

    VectorComplexFloat4 c1{lhs.Get(0), lhs.Get(0), lhs.Get(2), lhs.Get(2)};
    VectorComplexFloat4 c2{rhs.Get(0), rhs.Get(1), rhs.Get(0), rhs.Get(1)};
    VectorComplexFloat4 s1 = c1 * c2;

    VectorComplexFloat4 c3{lhs.Get(1), lhs.Get(1), lhs.Get(3), lhs.Get(3)};
    VectorComplexFloat4 c4{rhs.Get(2), rhs.Get(3), rhs.Get(2), rhs.Get(3)};
    VectorComplexFloat4 s2 = c3 * c4;

    return s1 + s2;
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, std::complex<float> rhs) noexcept {
    return lhs.data_ * rhs;
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator*(
      std::complex<float> lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return rhs * lhs;
  }

  AVX_TARGET DiagonalMatrixComplexFloat2x2 Diagonal() const noexcept {
    return VectorComplexFloat2{data_.Get(0), data_.Get(3)};
  }

  AVX_TARGET bool Invert() noexcept {
    // Note it would be possible to optimize further.
    VectorComplexFloat2 a{data_.Get(0), data_.Get(1)};
    VectorComplexFloat2 b{data_.Get(3), data_.Get(2)};
    VectorComplexFloat2 c = a * b;

    std::complex<float> d = c.Get(0) - c.Get(1);
    if (d == std::complex<float>{}) return false;

    const float n = std::norm(d);
    d.imag(-d.imag());
    __m256 reciprocal = _mm256_setr_ps(d.real(), d.imag(), d.real(), d.imag(),
                                       d.real(), d.imag(), d.real(), d.imag());
    reciprocal = _mm256_div_ps(reciprocal, _mm256_set1_ps(n));

    // std::swap(data[0],data[3]);
    // Using the fact that extracting as a double, the value has the number of
    // bits is the same as a complex float.
    __m256d data = _mm256_castps_pd(static_cast<__m256>(data_));
    data = _mm256_permute4x64_pd(data, 0b00'10'01'11);
    __m256 result = _mm256_castpd_ps(data);

    // data[0] = data[0]
    // data[1] = -data[1]
    // data[2] = -data[2]
    // data[3] = data[3]
    __m256 mask = _mm256_setr_ps(0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0);
    result = _mm256_xor_ps(result, mask);

    data_ = VectorComplexFloat4{result} * VectorComplexFloat4{reciprocal};

    return true;
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, DiagonalMatrixComplexFloat2x2 rhs) noexcept {
    // The multiplication of a 2x2 matrix M with a diagonal 2x2 matrix D is
    // not commutative and R = M * D can be written as
    // r[0][0] = m[0][0] * d[0][0] -> lhs[0] * rhs[0]
    // r[0][1] = m[0][1] * d[1][1] -> lhs[1] * rhs[1]
    // r[1][0] = m[1][0] * d[0][0] -> lhs[2] * rhs[0]
    // r[1][1] = m[1][1] * d[1][1] -> lhs[3] * rhs[1]

    __m128 lo = static_cast<__m128>(rhs);
    return lhs.data_ * VectorComplexFloat4{_mm256_set_m128(lo, lo)};
  }

  AVX_TARGET friend MatrixComplexFloat2x2 operator*(
      DiagonalMatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    // The multiplication of a diagonal 2x2 matrix D with a 2x2 matrix M is
    // not commutative and R = D * M can be written as:
    // r[0][0] = d[0][0] * m[0][0] -> lhs[0] * rhs[0]
    // r[0][1] = d[0][0] * m[0][1] -> lhs[0] * rhs[1]
    // r[1][0] = d[1][1] * m[1][0] -> lhs[1] * rhs[2]
    // r[1][1] = d[1][1] * m[1][1] -> lhs[1] * rhs[3]

    __m128 lo = static_cast<__m128>(lhs);
    __m256 l = _mm256_set_m128(lo, lo);

    // The vector contains the element 0 1 0 1 but it needs 0 0 1 1.  The
    // intrinsic _mm256_permute_ps does the same permutation for both 128-bit
    // lanes, which makes it impossible to get the desired output.  Since the
    // two floats are a complex pair it's possible to pretend they are 64-bit
    // entities. This makes it possible to use the _mm256_permute_pd intrinsic.
    // This intrinsic has separate control bits for both 128-bit lanes. (It
    // does not allow lane crossing, but that's not needed.)
    l = _mm256_castpd_ps(_mm256_permute_pd(_mm256_castps_pd(l), 0b00'00'11'00));
    return VectorComplexFloat4{l} * rhs.data_;
  }

  AVX_TARGET friend bool operator==(MatrixComplexFloat2x2 lhs,
                                    MatrixComplexFloat2x2 rhs) noexcept {
    return lhs.data_ == rhs.data_;
  }

  AVX_TARGET friend std::ostream& operator<<(std::ostream& output,
                                             MatrixComplexFloat2x2 value) {
    output << "[{" << value.Get(0) << ", " << value.Get(1) << "}, {"
           << value.Get(2) << ", " << value.Get(3) << "}]";
    return output;
  }

  AVX_TARGET MatrixComplexFloat2x2 HermTranspose() const noexcept {
    return Transpose().Conjugate();
  }

  template <typename Scalar>
  AVX_TARGET MatrixComplexFloat2x2& operator/=(Scalar scalar) {
    // TODO could be explicitly vectorized
    data_ = VectorComplexFloat4(data_.Get(0) / scalar, data_.Get(1) / scalar,
                                data_.Get(2) / scalar, data_.Get(3) / scalar);
    return *this;
  }

  /** @returns A^H * A. */
  AVX_TARGET MatrixComplexFloat2x2 HermitianSquare() const noexcept {
    // An explicit AVX version could be more efficient, but the use of this
    // method is rare so we have used a simpler implementation so far.
    return HermTranspose() * (*this);
  }

  AVX_TARGET std::complex<float> DoubleDot(MatrixComplexFloat2x2 rhs) const {
    return (*this).Get(0) * rhs.Get(0) + (*this).Get(1) * rhs.Get(1) +
           (*this).Get(2) * rhs.Get(2) + (*this).Get(3) * rhs.Get(3);
  }

  /**
   * Matrix multiplication, alias for the overloaded * operator and thus equally
   * computationally efficient.
   */
  AVX_TARGET MatrixComplexFloat2x2
  Multiply(MatrixComplexFloat2x2 rhs) const noexcept {
    return *this * rhs;
  }

  /**
   * Matrix multiplication of internal matrix with Hermitian transpose of input
   * matrix, i.e. returns A * B^H. Note that this is preferred over combining
   * operator* with HermTranspose() for computational efficiency in the non-AVX
   * case.
   */
  AVX_TARGET MatrixComplexFloat2x2
  MultiplyHerm(MatrixComplexFloat2x2 rhs) const noexcept {
    return *this * rhs.HermTranspose();
  }

  /**
   * Matrix multiplication Hermitian transpose of internal matrix with input
   * matrix, i.e. returns A^H * B. Note that this is preferred over combining
   * operator* with HermTranspose() for computational efficiency in the non-AVX
   * case.
   */
  AVX_TARGET MatrixComplexFloat2x2
  HermThenMultiply(MatrixComplexFloat2x2 rhs) const noexcept {
    return HermTranspose() * rhs;
  }

  /**
   * Matrix multiplication of Hermitian transposes of both the internal matrix
   * and the input matrix, i.e. returns A^H * B^H. Note that this is preferred
   * over combining operator* with HermTranspose() for computational efficiency
   * in the non-AVX case.
   */
  AVX_TARGET MatrixComplexFloat2x2
  HermThenMultiplyHerm(MatrixComplexFloat2x2 rhs) const noexcept {
    return HermTranspose() * rhs.HermTranspose();
  }

  /// Element-wise product
  AVX_TARGET friend MatrixComplexFloat2x2 ElementProduct(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) {
    return MatrixComplexFloat2x2(lhs.data_ * rhs.data_);
  }

 private:
  VectorComplexFloat4 data_;
  AVX_TARGET friend float Norm(MatrixComplexFloat2x2 matrix) noexcept;
};

/// MC2x2Base compatibility wrapper.
AVX_TARGET inline MatrixComplexFloat2x2 HermTranspose(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.HermTranspose();
}

/** Returns the sum of the diagonal elements. */
AVX_TARGET inline std::complex<float> Trace(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Trace();
}

AVX_TARGET inline DiagonalMatrixComplexFloat2x2 Diagonal(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Diagonal();
}

/// @returns the Frobenius norm of the matrix.
AVX_TARGET inline float Norm(MatrixComplexFloat2x2 matrix) noexcept {
  // This uses the same basic idea as MatrixComplexDouble2x2::Norm except
  // that the underlying data is stored in one __m256d value.

  // Note this function seems slower than expected.
  // MatrixComplexDouble2x2::Norm is faster than this function. It is
  // still faster than the scalar version. It would be nice to improve
  // this in the future.
  //
  //
  __m256 tmp = static_cast<__m256>(matrix.data_);
  tmp *= tmp;

  // For the addition we deviate from the double version.
  // | a     | b     | c     | d     |
  // | e     | f     | g     | h     |
  // ---------------------------------+
  // | a + e | b + f | c + g | d + h |

  __m128 ret = _mm256_castps256_ps128(tmp);
  ret += _mm256_extractf128_ps(tmp, 1);

  // | a + e         | b + f         |
  // | c + g         | d + h         |
  // ---------------------------------+
  // | a + e + c + g | b + f + d + h |

  // ret = _mm_hadd_ps(ret, ret);

  ret += _mm_movehl_ps(ret, ret);

  // | a + e + c + g                 |
  // | b + f + d + h                 |
  // ---------------------------------+
  // | a + e + c + g + b + f + d + h |

  return ret[0] + ret[1];
}

AVX_TARGET inline float SumOfAbsolute(MatrixComplexFloat2x2 matrix) {
  return std::abs(matrix.Get(0)) + std::abs(matrix.Get(1)) +
         std::abs(matrix.Get(2)) + std::abs(matrix.Get(3));
}

}  // namespace aocommon::avx

#endif  // AOCOMMON_AVX256_MATRIX_COMPLEX_FLOAT_2X2_H
