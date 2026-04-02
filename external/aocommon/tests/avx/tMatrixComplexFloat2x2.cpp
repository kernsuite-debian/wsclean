// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "aocommon/avx/AvxMacros.h"

#ifdef USE_AVX_MATRIX

#include "aocommon/avx/MatrixComplexFloat2x2.h"
#include "aocommon/avx/MatrixComplexDouble2x2.h"
#include "aocommon/scalar/matrix2x2.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

using AvxMC2x2F = aocommon::avx::MatrixComplexFloat2x2;
using ScalarMC2x2F = aocommon::scalar::MC2x2Base<float>;

inline AvxMC2x2F ToAvx(const ScalarMC2x2F& input) {
  return AvxMC2x2F(input.Get(0), input.Get(1), input.Get(2), input.Get(3));
}
inline ScalarMC2x2F ToScalar(AvxMC2x2F input) {
  return ScalarMC2x2F(input.Get(0), input.Get(1), input.Get(2), input.Get(3));
}

BOOST_AUTO_TEST_SUITE(matrix_complex_float_2x2)

static_assert(std::is_default_constructible_v<AvxMC2x2F>);
static_assert(std::is_nothrow_destructible_v<AvxMC2x2F>);
static_assert(std::is_nothrow_copy_constructible_v<AvxMC2x2F>);
static_assert(std::is_nothrow_move_constructible_v<AvxMC2x2F>);
static_assert(std::is_nothrow_copy_assignable_v<AvxMC2x2F>);
static_assert(std::is_nothrow_move_assignable_v<AvxMC2x2F>);

BOOST_AUTO_TEST_CASE(construct_zero_initialized) {
  static_assert(std::is_nothrow_default_constructible_v<AvxMC2x2F>);
  const AvxMC2x2F result = AvxMC2x2F::Zero();

  BOOST_TEST(result.Get(0) == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(constructor_vector_complex_float_4) {
  static_assert(
      std::is_nothrow_constructible_v<AvxMC2x2F,
                                      aocommon::avx::VectorComplexFloat4>);

  const AvxMC2x2F result{aocommon::avx::VectorComplexFloat4{
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}}};

  BOOST_TEST(result.Get(0) == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_4_complex_floats) {
  static_assert(
      std::is_nothrow_constructible_v<AvxMC2x2F, std::complex<float>,
                                      std::complex<float>, std::complex<float>,
                                      std::complex<float>>);

  const AvxMC2x2F result{
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};

  BOOST_TEST(result.Get(0) == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_float_pointer) {
  static_assert(
      std::is_nothrow_constructible_v<AvxMC2x2F, const std::complex<float>[4]>);
  const std::complex<float> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const AvxMC2x2F result{input};

  BOOST_TEST(result.Get(0) == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_complex_double_pointer) {
  static_assert(std::is_nothrow_constructible_v<AvxMC2x2F,
                                                const std::complex<double>[4]>);
  const std::complex<double> input[] = {
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const AvxMC2x2F result{input};

  BOOST_TEST(result.Get(0) == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(constructor_MatrixComplexDouble2x2) {
  static_assert(std::is_nothrow_constructible_v<
                AvxMC2x2F, const aocommon::avx::MatrixComplexDouble2x2&>);
  const aocommon::avx::MatrixComplexDouble2x2 input{
      {-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};
  const AvxMC2x2F result{input};

  BOOST_TEST(result.Get(0) == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result.Get(1) == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result.Get(2) == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result.Get(3) == (std::complex<float>{1.5, -1.5}));
}

BOOST_AUTO_TEST_CASE(to_avx) {
  const ScalarMC2x2F input({-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0},
                           {1.5, -1.5});
  const AvxMC2x2F result = ToAvx(input);

  BOOST_TEST(input.Get(0) == ToScalar(result).Get(0));
  BOOST_TEST(input.Get(1) == ToScalar(result).Get(1));
  BOOST_TEST(input.Get(2) == ToScalar(result).Get(2));
  BOOST_TEST(input.Get(3) == ToScalar(result).Get(3));
}

BOOST_AUTO_TEST_CASE(conjugate) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.Conjugate()));

  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}.Conjugate()),
      (AvxMC2x2F{{1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}}.Conjugate()),
      (AvxMC2x2F{{-1.0, -2.0}, {10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}}.Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {10, -11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, -11}, {100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, 11}, {100, -101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {-100, 101}, {1000, 1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, -11}, {-100, -101}, {1000, -1001}}));
  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, 11}, {-100, 101}, {1000, -1001}}));

  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}));
  BOOST_CHECK_EQUAL(
      (AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}
           .Conjugate()),
      (AvxMC2x2F{{-1.0, 2.0}, {-10, 11}, {-100, 101}, {-1000, 1001}}));
}

BOOST_AUTO_TEST_CASE(transpose) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.Transpose()));
  const ScalarMC2x2F input{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const AvxMC2x2F expected = ToAvx(input.Transpose());

  const AvxMC2x2F result = ToAvx(input).Transpose();

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(invert) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.Invert()));
  const ScalarMC2x2F input{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  ScalarMC2x2F inverted(input);
  BOOST_REQUIRE(inverted.Invert());
  const AvxMC2x2F expected = ToAvx(inverted);

  AvxMC2x2F result = ToAvx(input);
  BOOST_REQUIRE(result.Invert());

  BOOST_CHECK_CLOSE(result.Get(0).real(), expected.Get(0).real(), 3e-3);
  BOOST_CHECK_CLOSE(result.Get(0).imag(), expected.Get(0).imag(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(1).real(), expected.Get(1).real(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(1).imag(), expected.Get(1).imag(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(2).real(), expected.Get(2).real(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(2).imag(), expected.Get(2).imag(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(3).real(), expected.Get(3).real(), 1e-3);
  BOOST_CHECK_CLOSE(result.Get(3).imag(), expected.Get(3).imag(), 1e-3);
}

BOOST_AUTO_TEST_CASE(hermitian_transpose) {
  static_assert(
      noexcept(AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}
                   .HermTranspose()));
  const std::vector<ScalarMC2x2F> inputs{
      ScalarMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};

  for (const ScalarMC2x2F& input : inputs) {
    const AvxMC2x2F expected = ToAvx(input.HermTranspose());

    const AvxMC2x2F result = ToAvx(input).HermTranspose();

    BOOST_CHECK_EQUAL(result, expected);
  }
}

BOOST_AUTO_TEST_CASE(norm) {
  static_assert(noexcept(
      Norm(AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)})));
  const std::vector<AvxMC2x2F> inputs{
      AvxMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};

  for (const AvxMC2x2F& input : inputs)
    BOOST_CHECK_EQUAL(Norm(input), 2022428.0f);
}

BOOST_AUTO_TEST_CASE(trace) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.Trace()));
  const std::vector<ScalarMC2x2F> inputs{
      ScalarMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      ScalarMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};

  for (const ScalarMC2x2F& input : inputs) {
    const std::complex<float> expected = Trace(input);
    const std::complex<float> result = ToAvx(input).Trace();

    BOOST_CHECK_EQUAL(result, expected);
  }
}

BOOST_AUTO_TEST_CASE(assign_to_float) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.AssignTo(
          static_cast<std::complex<float>*>(nullptr))));
  const AvxMC2x2F input{{-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};

  std::vector<std::complex<float>> result(6);
  input.AssignTo(std::addressof(result[1]));

  BOOST_TEST(result[0] == (std::complex<float>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<float>{-1.0, 1.0}));
  BOOST_TEST(result[2] == (std::complex<float>{3.75, -3.75}));
  BOOST_TEST(result[3] == (std::complex<float>{99.0, -99.0}));
  BOOST_TEST(result[4] == (std::complex<float>{1.5, -1.5}));
  BOOST_TEST(result[5] == (std::complex<float>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(assign_to_double) {
  static_assert(noexcept(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}.AssignTo(
          static_cast<std::complex<double>*>(nullptr))));
  const AvxMC2x2F input{{-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};

  std::vector<std::complex<double>> result(6);
  input.AssignTo(std::addressof(result[1]));

  BOOST_TEST(result[0] == (std::complex<double>{0.0, 0.0}));
  BOOST_TEST(result[1] == (std::complex<double>{-1.0, 1.0}));
  BOOST_TEST(result[2] == (std::complex<double>{3.75, -3.75}));
  BOOST_TEST(result[3] == (std::complex<double>{99.0, -99.0}));
  BOOST_TEST(result[4] == (std::complex<double>{1.5, -1.5}));
  BOOST_TEST(result[5] == (std::complex<double>{0.0, 0.0}));
}

BOOST_AUTO_TEST_CASE(Unity) {
  static_assert(noexcept(AvxMC2x2F::Unity()));

  BOOST_CHECK_EQUAL(ToAvx(ScalarMC2x2F::Unity()), AvxMC2x2F::Unity());
}

BOOST_AUTO_TEST_CASE(NaN) {
  static_assert(noexcept(AvxMC2x2F::NaN()));

  const AvxMC2x2F value = AvxMC2x2F::NaN();

  BOOST_TEST(std::isnan(value.Get(0).real()));
  BOOST_TEST(std::isnan(value.Get(0).imag()));
  BOOST_TEST(std::isnan(value.Get(1).real()));
  BOOST_TEST(std::isnan(value.Get(1).imag()));
  BOOST_TEST(std::isnan(value.Get(2).real()));
  BOOST_TEST(std::isnan(value.Get(2).imag()));
  BOOST_TEST(std::isnan(value.Get(3).real()));
  BOOST_TEST(std::isnan(value.Get(3).imag()));
}

BOOST_AUTO_TEST_CASE(operator_plus_equal) {
  AvxMC2x2F r{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const AvxMC2x2F value{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  r += value;
  static_assert(noexcept(r += value));

  BOOST_TEST(r.Get(0) == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r.Get(1) == (std::complex<float>{50.0, 55.0}));
  BOOST_TEST(r.Get(2) == (std::complex<float>{500., 505.0}));
  BOOST_TEST(r.Get(3) == (std::complex<float>{5000., 5005.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus_equal) {
  AvxMC2x2F r{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const AvxMC2x2F value{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  r -= value;
  static_assert(noexcept(r -= value));

  BOOST_TEST(r.Get(0) == (std::complex<float>{-3, -6}));
  BOOST_TEST(r.Get(1) == (std::complex<float>{-30, -33}));
  BOOST_TEST(r.Get(2) == (std::complex<float>{-300, -303}));
  BOOST_TEST(r.Get(3) == (std::complex<float>{-3000, -3003}));
}

BOOST_AUTO_TEST_CASE(operator_plus) {
  const AvxMC2x2F lhs{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const AvxMC2x2F r = lhs + rhs;
  static_assert(noexcept(lhs + rhs));

  BOOST_TEST(r.Get(0) == (std::complex<float>{5.0, 10.0}));
  BOOST_TEST(r.Get(1) == (std::complex<float>{50.0, 55.0}));
  BOOST_TEST(r.Get(2) == (std::complex<float>{500., 505.0}));
  BOOST_TEST(r.Get(3) == (std::complex<float>{5000., 5005.0}));
}

BOOST_AUTO_TEST_CASE(operator_minus) {
  const AvxMC2x2F lhs{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const AvxMC2x2F r = lhs - rhs;
  static_assert(noexcept(lhs - rhs));

  BOOST_TEST(r.Get(0) == (std::complex<float>{-3, -6}));
  BOOST_TEST(r.Get(1) == (std::complex<float>{-30, -33}));
  BOOST_TEST(r.Get(2) == (std::complex<float>{-300, -303}));
  BOOST_TEST(r.Get(3) == (std::complex<float>{-3000, -3003}));
}

BOOST_AUTO_TEST_CASE(multiply) {
  {
    const AvxMC2x2F lhs{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

    const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

    static_assert(noexcept(lhs * rhs));
    const AvxMC2x2F r = lhs * rhs;

    BOOST_CHECK_CLOSE(r.Get(0).real(), -456, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(0).imag(), 8456, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).real(), -4092, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).imag(), 84164, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).real(), -4812, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).imag(), 805604, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).real(), -8448, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).imag(), 8016440, 1e-6);
  }

  // Emulate Matrix * Diagonal Matrix
  {
    const AvxMC2x2F lhs{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

    const AvxMC2x2F rhs{{4, 8}, {0, 0}, {0, 0}, {4000, 4004}};

    const AvxMC2x2F r = lhs * rhs;

    BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).real(), -4044, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).imag(), 84040, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).real(), -408, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).imag(), 1204, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).real(), -8004, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).imag(), 8008000, 1e-6);
  }

  // Emulate Diagonal Matrix * Matrix
  {
    const AvxMC2x2F lhs{{1.0, 2.0}, {0, 0}, {0, 0}, {1000, 1001}};

    const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

    static_assert(noexcept(lhs * rhs));
    const AvxMC2x2F r = lhs * rhs;

    BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).real(), -48, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(1).imag(), 124, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).real(), -4404, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(2).imag(), 804400, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).real(), -8004, 1e-6);
    BOOST_CHECK_CLOSE(r.Get(3).imag(), 8008000, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(multiply_matrix_and_value) {
  const AvxMC2x2F lhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  const std::complex<float> rhs{1.0, 2.0};

  static_assert(noexcept(lhs * rhs));
  const AvxMC2x2F r = lhs * rhs;

  BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).real(), -48, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).imag(), 124, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).real(), -408, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).imag(), 1204, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).real(), -4008, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).imag(), 12004, 1e-6);
}

BOOST_AUTO_TEST_CASE(multiply_value_and_matrix) {
  const std::complex<float> lhs{1.0, 2.0};

  const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  static_assert(noexcept(lhs * rhs));
  const AvxMC2x2F r = lhs * rhs;

  BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).real(), -48, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).imag(), 124, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).real(), -408, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).imag(), 1204, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).real(), -4008, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).imag(), 12004, 1e-6);
}

BOOST_AUTO_TEST_CASE(multiply_diagonal_matrix_and_matrix) {
  const aocommon::avx::DiagonalMatrixComplexFloat2x2 lhs{{1.0, 2.0},
                                                         {1000, 1001}};

  const AvxMC2x2F rhs{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  static_assert(noexcept(lhs * rhs));
  const AvxMC2x2F r = lhs * rhs;

  BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).real(), -48, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).imag(), 124, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).real(), -4404, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).imag(), 804400, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).real(), -8004, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).imag(), 8008000, 1e-6);
}

BOOST_AUTO_TEST_CASE(multiply_matrix_and_diagonal_matrix) {
  const AvxMC2x2F lhs{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}};

  const aocommon::avx::DiagonalMatrixComplexFloat2x2 rhs{{4, 8}, {4000, 4004}};

  static_assert(noexcept(lhs * rhs));
  const AvxMC2x2F r = lhs * rhs;

  BOOST_CHECK_CLOSE(r.Get(0).real(), -12, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(0).imag(), 16, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).real(), -4044, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(1).imag(), 84040, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).real(), -408, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(2).imag(), 1204, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).real(), -8004, 1e-6);
  BOOST_CHECK_CLOSE(r.Get(3).imag(), 8008000, 1e-6);
}

BOOST_AUTO_TEST_CASE(equal) {
  static_assert(
      noexcept(AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)} ==
               AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)}));

  BOOST_TEST((AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
              AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}, {0.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {42.0, 0.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
  BOOST_TEST(!(AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 42.0}} ==
               AvxMC2x2F{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}));
}

BOOST_AUTO_TEST_CASE(output) {
  const AvxMC2x2F input{{-1.0, 1.0}, {3.75, -3.75}, {99.0, -99.0}, {1.5, -1.5}};

  std::stringstream result;
  result << input;

  BOOST_CHECK_EQUAL(result.str(),
                    "[{(-1,1), (3.75,-3.75)}, {(99,-99), (1.5,-1.5)}]");
}

BOOST_AUTO_TEST_CASE(herm_transpose) {
  static_assert(noexcept(HermTranspose(
      AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)})));
  const std::vector<AvxMC2x2F> inputs{
      AvxMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};

  for (const AvxMC2x2F& input : inputs)
    BOOST_CHECK_EQUAL(input.HermTranspose(), HermTranspose(input));
}

BOOST_AUTO_TEST_CASE(non_member_norm) {
  static_assert(noexcept(
      Norm(AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)})));
  const std::vector<AvxMC2x2F> inputs{
      AvxMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};
}

BOOST_AUTO_TEST_CASE(non_member_trace) {
  static_assert(noexcept(
      Trace(AvxMC2x2F{static_cast<const std::complex<float>*>(nullptr)})));
  const std::vector<AvxMC2x2F> inputs{
      AvxMC2x2F{{1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, 2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, 11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, 101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, 1001}},
      AvxMC2x2F{{-1.0, -2.0}, {-10, -11}, {-100, -101}, {-1000, -1001}}};

  for (const AvxMC2x2F& input : inputs)
    BOOST_CHECK_EQUAL(Trace(input), input.Trace());
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_AVX_MATRIX
