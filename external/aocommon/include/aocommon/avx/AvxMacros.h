#ifndef AOCOMMON_AVX_MACROS_H_
#define AOCOMMON_AVX_MACROS_H_

#if defined(__AVX2__) && defined(__FMA__)

#define USE_AVX_MATRIX

/**
 * Functions in the avx namespace need to use a target attribute
 * like @c [[gnu::target("avx2,fma")]] for function multi-versioning.
 *
 * This is currently not yet enabled, because the current solution requires
 * __AVX2__ to be defined, but because the PORTABLE option disabled
 * -march=native, that macro does not get defined in portable builds. The
 * aoflagger code has an example (by Marcel Loose) that avoids this; it may be
 * possible to implement it similarly here.
 *
 * For now, it's disabled until FMV is properly implemented.
 */
//#define AVX_TARGET [[gnu::target("avx2,fma")]]
#define AVX_TARGET

#endif

#endif  // AOCOMMON_AVX_MACROS_H_
