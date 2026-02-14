#ifndef AOCOMMON_MATRIX_2X2_H_
#define AOCOMMON_MATRIX_2X2_H_

#include "avx/AvxMacros.h"

#ifdef USE_AVX_MATRIX
#include "avx/MatrixComplexFloat2x2.h"
#include "avx/MatrixComplexDouble2x2.h"
#endif
#include "scalar/matrix2x2.h"

namespace aocommon {

#ifdef USE_AVX_MATRIX
using MC2x2 = avx::MatrixComplexDouble2x2;
using MC2x2F = avx::MatrixComplexFloat2x2;
#else
using MC2x2 = scalar::MC2x2Base<double>;
using MC2x2F = scalar::MC2x2Base<float>;
#endif
using Matrix2x2 = scalar::Matrix2x2;

/// These asserts are necessary for the casts below
static_assert(sizeof(MC2x2) == 4 * sizeof(std::complex<double>));
static_assert(sizeof(MC2x2F) == 4 * sizeof(std::complex<float>));
static_assert(std::is_standard_layout_v<MC2x2>);
static_assert(std::is_standard_layout_v<MC2x2F>);

/**
 * Cast the double-precision MC2x2 to a std::complex<double> pointer.
 * This makes use of the fact that:
 * - A trivial class can be cast to its first member if the class is standard
 * layout. This is statically asserted above. See also
 * https://stackoverflow.com/questions/49333703/is-cast-a-pointer-points-to-class-to-its-first-member-illegal
 * - The first member is either a std::complex<double> (without AVX) or __m256d
 * (with AVX). The first case may be cast because of the argument under point 1.
 * The second case seems to be a somewhat gcc special case: __m256d has the
 * property that strict aliasing does not apply to these avx types, and it is
 * therefore allowed. Because this cast is somewhat dubious and undesirable, the
 * use of this function should be minimized, e.g. only to interface the matrix
 * class with external libraries to avoid having to make an unnecessary copy.
 */
inline std::complex<double>* DubiousDComplexPointerCast(MC2x2& matrix) {
  return reinterpret_cast<std::complex<double>*>(&matrix);
}
inline const std::complex<double>* DubiousDComplexPointerCast(
    const MC2x2& matrix) {
  return reinterpret_cast<const std::complex<double>*>(&matrix);
}
/** Same as @ref DubiousDComplexPointerCast, but for MC2x2F. */
inline std::complex<float>* DubiousComplexPointerCast(MC2x2F& matrix) {
  return reinterpret_cast<std::complex<float>*>(&matrix);
}
inline const std::complex<float>* DubiousComplexPointerCast(
    const MC2x2F& matrix) {
  return reinterpret_cast<const std::complex<float>*>(&matrix);
}

}  // namespace aocommon

#endif
