// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "aocommon/avx/AvxMacros.h"

#ifdef USE_AVX_MATRIX

#include "aocommon/avx/DiagonalMatrixComplexFloat2x2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

using AvxDiagMC2x2F = aocommon::avx::DiagonalMatrixComplexFloat2x2;

BOOST_AUTO_TEST_SUITE(DiagonalMatrixComplexFloat2x2)

static_assert(std::is_default_constructible_v<AvxDiagMC2x2F>);
static_assert(std::is_nothrow_destructible_v<AvxDiagMC2x2F>);
static_assert(std::is_nothrow_copy_constructible_v<AvxDiagMC2x2F>);
static_assert(std::is_nothrow_move_constructible_v<AvxDiagMC2x2F>);
static_assert(std::is_nothrow_copy_assignable_v<AvxDiagMC2x2F>);
static_assert(std::is_nothrow_move_assignable_v<AvxDiagMC2x2F>);

BOOST_AUTO_TEST_CASE(construct_zero_initialized) {
  static_assert(std::is_nothrow_default_constructible_v<AvxDiagMC2x2F>);
  static_assert(noexcept(AvxDiagMC2x2F::Zero()));
  BOOST_TEST(AvxDiagMC2x2F::Zero()[0] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(AvxDiagMC2x2F::Zero()[1] == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(constructor_2_complex_float) {
  static_assert(
      std::is_nothrow_constructible_v<AvxDiagMC2x2F, std::complex<float>,
                                      std::complex<float>>);

  const AvxDiagMC2x2F result{std::complex<float>{-1.0, 1.0},
                             std::complex<float>{3.75, -3.75}};

  BOOST_TEST(result[0] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[1] == (std::complex<float>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_float_pointer) {
  static_assert(std::is_nothrow_constructible_v<AvxDiagMC2x2F,
                                                const std::complex<float>[4]>);
  const std::complex<float> input[] = {std::complex<float>{-1.0, 1.0},
                                       std::complex<float>{3.75, -3.75}};
  const AvxDiagMC2x2F result{input};

  BOOST_TEST((result[0] == std::complex<float>{-1.0, 1.0}));
  BOOST_TEST((result[1] == std::complex<float>{3.75, -3.75}));
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(
      noexcept(AvxDiagMC2x2F{static_cast<const std::complex<float>*>(nullptr)}
                   .Conjugate()));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{1.0, 2.0}, {10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{1.0, -2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, 2.0}, {10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, -2.0}, {10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.Conjugate()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(hermitian_transpose) {
  static_assert(
      noexcept(AvxDiagMC2x2F{static_cast<const std::complex<float>*>(nullptr)}
                   .HermTranspose()));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{1.0, 2.0}, {10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{1.0, -2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, 2.0}, {10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, -2.0}, {10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {10, -11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));

  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, 11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, -11}}));
  BOOST_CHECK_EQUAL((AvxDiagMC2x2F{{-1.0, -2.0}, {-10, -11}}.HermTranspose()),
                    (AvxDiagMC2x2F{{-1.0, 2.0}, {-10, 11}}));
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  AvxDiagMC2x2F r{{1.0, 2.0}, {10, 11}};

  const AvxDiagMC2x2F value{{4, 8}, {40, 44}};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r[0] == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r[1] == (std::complex<float>{50.0, 55.0}));
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(noexcept(
      AvxDiagMC2x2F{static_cast<const std::complex<float>*>(nullptr)} ==
      AvxDiagMC2x2F{static_cast<const std::complex<float>*>(nullptr)}));

  BOOST_TEST((AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}} ==
              AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxDiagMC2x2F{{42.0, 0.0}, {0.0, 0.0}} ==
               AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxDiagMC2x2F{{0.0, 42.0}, {0.0, 0.0}} ==
               AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxDiagMC2x2F{{0.0, 0.0}, {42.0, 0.0}} ==
               AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 42.0}} ==
               AvxDiagMC2x2F{{0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_AVX_MATRIX
