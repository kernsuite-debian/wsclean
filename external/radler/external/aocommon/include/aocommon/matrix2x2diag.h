#ifndef AOCOMMON_MATRIX_2X2_DIAG_H_
#define AOCOMMON_MATRIX_2X2_DIAG_H_

#include "avx/AvxMacros.h"

#ifdef USE_AVX_MATRIX
#include "avx/DiagonalMatrixComplexDouble2x2.h"
#include "avx/DiagonalMatrixComplexFloat2x2.h"
#else
#include "scalar/matrix2x2diag.h"
#endif

namespace aocommon {

#ifdef USE_AVX_MATRIX
using MC2x2Diag = avx::DiagonalMatrixComplexDouble2x2;
using MC2x2FDiag = avx::DiagonalMatrixComplexFloat2x2;
#else
using MC2x2Diag = scalar::MC2x2DiagBase<double>;
using MC2x2FDiag = scalar::MC2x2DiagBase<float>;
#endif

}  // namespace aocommon

#endif
