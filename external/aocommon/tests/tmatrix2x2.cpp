#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

#include <boost/test/unit_test.hpp>

#include <iostream>

using aocommon::Matrix2x2;
using aocommon::MC2x2;
using aocommon::MC2x2Diag;
using aocommon::MC2x2F;
using aocommon::MC2x2FDiag;

namespace {
template <typename MType>
void CheckClose(const MType& a, const MType& b, double tolerance = 1e-6) {
  // Writing this out makes it easier to debug compared to using a for loop.
  BOOST_CHECK_CLOSE(a.Get(0).real(), b.Get(0).real(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(0).imag(), b.Get(0).imag(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(1).real(), b.Get(1).real(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(1).imag(), b.Get(1).imag(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(2).real(), b.Get(2).real(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(2).imag(), b.Get(2).imag(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(3).real(), b.Get(3).real(), tolerance);
  BOOST_CHECK_CLOSE(a.Get(3).imag(), b.Get(3).imag(), tolerance);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(matrix2x2)

BOOST_AUTO_TEST_CASE(construct_zero_initialized) {
  const MC2x2 m = MC2x2::Zero();
  BOOST_CHECK_EQUAL(m.Get(0).real(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(0).imag(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(1).real(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(1).imag(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(2).real(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(2).imag(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(3).real(), 0.0);
  BOOST_CHECK_EQUAL(m.Get(3).imag(), 0.0);
}

BOOST_AUTO_TEST_CASE(construct_initializer_list) {
  const MC2x2 a = {1.0, 2.0, 3.0, 4.0};
  BOOST_CHECK_EQUAL(a.Get(0).real(), 1.0);
  BOOST_CHECK_EQUAL(a.Get(0).imag(), 0.0);
  BOOST_CHECK_EQUAL(a.Get(1).real(), 2.0);
  BOOST_CHECK_EQUAL(a.Get(1).imag(), 0.0);
  BOOST_CHECK_EQUAL(a.Get(2).real(), 3.0);
  BOOST_CHECK_EQUAL(a.Get(2).imag(), 0.0);
  BOOST_CHECK_EQUAL(a.Get(3).real(), 4.0);
  BOOST_CHECK_EQUAL(a.Get(3).imag(), 0.0);
  const MC2x2 b = {
      std::complex<double>{1.0, 2.0}, std::complex<double>{3.0, 4.0},
      std::complex<double>{5.0, 6.0}, std::complex<double>{7.0, 8.0}};
  BOOST_CHECK_EQUAL(b.Get(0).real(), 1.0);
  BOOST_CHECK_EQUAL(b.Get(0).imag(), 2.0);
  BOOST_CHECK_EQUAL(b.Get(1).real(), 3.0);
  BOOST_CHECK_EQUAL(b.Get(1).imag(), 4.0);
  BOOST_CHECK_EQUAL(b.Get(2).real(), 5.0);
  BOOST_CHECK_EQUAL(b.Get(2).imag(), 6.0);
  BOOST_CHECK_EQUAL(b.Get(3).real(), 7.0);
  BOOST_CHECK_EQUAL(b.Get(3).imag(), 8.0);
}

BOOST_AUTO_TEST_CASE(copy_construct) {
  const MC2x2F source({3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0});
  const MC2x2F copy(source);
  constexpr size_t m[] = {0, 2, 1, 3};
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_EQUAL(copy.Get(m[i]), source.Get(m[i]));
  }
  BOOST_CHECK_EQUAL(copy, source);
}

BOOST_AUTO_TEST_CASE(move_construct) {
  const MC2x2F original({3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0});
  MC2x2F source(original);
  const MC2x2F dest(std::move(source));
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_EQUAL(dest.Get(i), original.Get(i));
  }
  BOOST_CHECK_EQUAL(dest, original);
}

BOOST_AUTO_TEST_CASE(move_assign) {
  MC2x2F dest;
  MC2x2F source({3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0});
  dest = std::move(source);
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_EQUAL(dest.Get(i), source.Get(i));
  }
  BOOST_CHECK_EQUAL(dest, source);
}

BOOST_AUTO_TEST_CASE(copy_assign) {
  const MC2x2F source({3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0});
  MC2x2F copy(MC2x2F::Zero());
  copy = source;
  BOOST_CHECK_EQUAL(copy, source);
}

BOOST_AUTO_TEST_CASE(from_diagonal) {
  const aocommon::MC2x2Diag a({1.0, 2.0}, {3.0, 4.0});
  MC2x2 b(a);
  CheckClose(b, MC2x2{{1.0, 2.0}, {0.0}, {0.0}, {3.0, 4.0}});
}

BOOST_AUTO_TEST_CASE(add) {
  const MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const MC2x2 b({-5.0, 7.0}, {9.0, -11.0}, {-13.0, -15.0}, {17.0, 19.0});
  CheckClose(a + b,
             MC2x2{{-4.0, 9.0}, {12.0, -7.0}, {-8.0, -9.0}, {24.0, 27.0}});
}

BOOST_AUTO_TEST_CASE(add_diagonal) {
  const MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const aocommon::MC2x2Diag b({-5.0, 7.0}, {17.0, 19.0});
  CheckClose(a + b, MC2x2{{-4.0, 9.0}, {3.0, 4.0}, {5.0, 6.0}, {24.0, 27.0}});
}

BOOST_AUTO_TEST_CASE(add_lhs_diagonal) {
  const aocommon::MC2x2Diag a({-5.0, 7.0}, {17.0, 19.0});
  const MC2x2 b({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  CheckClose(a + b, MC2x2{{-4.0, 9.0}, {3.0, 4.0}, {5.0, 6.0}, {24.0, 27.0}});
}

BOOST_AUTO_TEST_CASE(assign_add_diagonal) {
  MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const aocommon::MC2x2Diag b({5.0, 7.0}, {9.0, 11.0});
  const MC2x2 result({6.0, 9.0}, {3.0, 4.0}, {5.0, 6.0}, {16.0, 19.0});
  CheckClose(a += b, result);
  CheckClose(a, result);
}

BOOST_AUTO_TEST_CASE(subtraction) {
  const MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const MC2x2 b({5.0, 7.0}, {9.0, 11.0}, {13.0, 15.0}, {17.0, 19.0});
  CheckClose(a - b,
             MC2x2{{-4.0, -5.0}, {-6.0, -7.0}, {-8.0, -9.0}, {-10.0, -11.0}});
}

BOOST_AUTO_TEST_CASE(subtract_diagonal) {
  const MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const aocommon::MC2x2Diag b({5.0, 7.0}, {9.0, 11.0});
  CheckClose(a - b, MC2x2({-4.0, -5.0}, {3.0, 4.0}, {5.0, 6.0}, {-2.0, -3.0}));
}

BOOST_AUTO_TEST_CASE(subtract_lhs_diagonal) {
  const aocommon::MC2x2Diag a({5.0, 7.0}, {9.0, 11.0});
  const MC2x2 b({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  CheckClose(a - b, MC2x2({4.0, 5.0}, {-3.0, -4.0}, {-5.0, -6.0}, {2.0, 3.0}));
}

BOOST_AUTO_TEST_CASE(assign_subtract_diagonal) {
  MC2x2 a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const aocommon::MC2x2Diag b({5.0, 7.0}, {9.0, 11.0});
  const MC2x2 result({-4.0, -5.0}, {3.0, 4.0}, {5.0, 6.0}, {-2.0, -3.0});
  CheckClose(a -= b, result);
  CheckClose(a, result);
}

BOOST_AUTO_TEST_CASE(scalar_multiplication) {
  const MC2x2 m({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});

  const std::complex<double> complex_scalar(5.0, -1.0);
  const MC2x2 expected_1({7.0, 9.0}, {19.0, 17.0}, {31.0, 25.0}, {43.0, 33.0});
  CheckClose(m * complex_scalar, expected_1);
  CheckClose(complex_scalar * m, expected_1);

  constexpr double real_scalar = -2.0;
  const MC2x2 expected_2({-2.0, -4.0}, {-6.0, -8.0}, {-10.0, -12.0},
                         {-14.0, -16.0});
  CheckClose(m * real_scalar, expected_2);
  CheckClose(real_scalar * m, expected_2);
}

BOOST_AUTO_TEST_CASE(multiply_complex_times_real_array) {
  MC2x2 a({1.0, 2.0}, {0, 0}, {0, 0}, {3.0, 4.0});
  // Flattened real 2x2 array
  const double r[4] = {10, 20, 30, 40};

  // Multiply
  MC2x2 b = a * r;
  CheckClose(b, MC2x2{{10, 20}, {20, 40}, {90, 120}, {120, 160}});

  // Multiply-assign
  a *= r;
  CheckClose(a, MC2x2{{10, 20}, {20, 40}, {90, 120}, {120, 160}});
}

BOOST_AUTO_TEST_CASE(multiply_full_and_diagonal) {
  const MC2x2F a({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0});
  const MC2x2FDiag b({-4.0, -5.0}, {-2.0, -3.0});
  const MC2x2F a_times_b{
      {6.0, -13.0}, {6.0, -17.0}, {10.0, -49.0}, {10, -37.0}};
  const MC2x2F b_times_a{{6.0, -13.0}, {8.0, -31.0}, {8.0, -27.0}, {10, -37.0}};
  // Multiply
  CheckClose(a * b, a_times_b);
  CheckClose(b * a, b_times_a);
  CheckClose(MC2x2(a) * MC2x2Diag(b), MC2x2(a_times_b));
  CheckClose(MC2x2Diag(b) * MC2x2(a), MC2x2(b_times_a));
}

BOOST_AUTO_TEST_CASE(complex_division_with_real) {
  MC2x2F a({4.0, 2.0}, {0, 40}, {12, 16}, {8.0, 4.0});

  // Divide and assign
  a /= 4.0f;
  CheckClose(a, MC2x2F{{1, 0.5}, {0, 10}, {3, 4}, {2, 1}});
}

BOOST_AUTO_TEST_CASE(assign_to) {
  MC2x2 a({1.0, 2.0}, {0, 0}, {0, 0}, {3.0, 4.0});
  std::complex<double> r1[4];

  // Assign to complex double buffer
  a.AssignTo(r1);
  BOOST_CHECK_CLOSE(r1[0].real(), 1, 1e-6);
  BOOST_CHECK_CLOSE(r1[0].imag(), 2, 1e-6);
  BOOST_CHECK_CLOSE(r1[3].real(), 3, 1e-6);
  BOOST_CHECK_CLOSE(r1[3].imag(), 4, 1e-6);

  // Assign to complex float buffer.
  std::complex<float> r2[4];
  a.AssignTo(r2);
  BOOST_CHECK_CLOSE(r2[0].real(), 1, 1e-6);
  BOOST_CHECK_CLOSE(r2[0].imag(), 2, 1e-6);
  BOOST_CHECK_CLOSE(r2[3].real(), 3, 1e-6);
  BOOST_CHECK_CLOSE(r2[3].imag(), 4, 1e-6);
}

BOOST_AUTO_TEST_CASE(hermitian_square) {
  MC2x2 a({1.0, 2.0}, {10.0, 11.0}, {20, 21}, {30.0, 31.0});
  BOOST_CHECK(a.HermitianSquare() == a.HermTranspose() * a);
}

BOOST_AUTO_TEST_CASE(real_hermitian_square_float) {
  MC2x2F a1({1.0, 2.0}, {0.0, 0.0}, {0, 0}, {30.0, 31.0});
  MC2x2FDiag a2({1.0, 2.0}, {30.0, 31.0});

  MC2x2F a1_herm_square = a1.HermitianSquare();
  std::array<float, 2> a2_real_herm_square = a2.RealHermitianSquare();

  BOOST_CHECK_CLOSE(std::real(a1_herm_square.Get(0)), a2_real_herm_square[0],
                    1e-12);
  BOOST_CHECK(std::real(a1_herm_square.Get(1)) == 0.0);
  BOOST_CHECK(std::real(a1_herm_square.Get(2)) == 0.0);
  BOOST_CHECK_CLOSE(std::real(a1_herm_square.Get(3)), a2_real_herm_square[1],
                    1e-12);
}

BOOST_AUTO_TEST_CASE(real_hermitian_square_double) {
  MC2x2 a1({1.0, 2.0}, {0.0, 0.0}, {0, 0}, {30.0, 31.0});
  MC2x2Diag a2({1.0, 2.0}, {30.0, 31.0});

  MC2x2 a1_herm_square = a1.HermitianSquare();
  std::array<double, 2> a2_real_herm_square = a2.RealHermitianSquare();

  BOOST_CHECK_CLOSE(std::real(a1_herm_square.Get(0)), a2_real_herm_square[0],
                    1e-12);
  BOOST_CHECK(std::real(a1_herm_square.Get(1)) == 0.0);
  BOOST_CHECK(std::real(a1_herm_square.Get(2)) == 0.0);
  BOOST_CHECK_CLOSE(std::real(a1_herm_square.Get(3)), a2_real_herm_square[1],
                    1e-12);
}

BOOST_AUTO_TEST_CASE(eigenvalue1) {
  double unit[4] = {1.0, 0.0, 0.0, 1.0};
  double e1, e2;
  Matrix2x2::EigenValues(unit, e1, e2);
  BOOST_CHECK_CLOSE(e1, 1.0, 1e-6);
  BOOST_CHECK_CLOSE(e2, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(eigenvalue2) {
  double unit[4] = {0.0, 1.0, -2.0, -3.0};
  double e1, e2;
  Matrix2x2::EigenValues(unit, e1, e2);
  if (e1 < e2) std::swap(e1, e2);
  BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
  BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(eigenvalue3) {
  double unit[4] = {0.0, -2.0, 1.0, -3.0};
  double e1, e2;
  Matrix2x2::EigenValues(unit, e1, e2);
  if (e1 < e2) std::swap(e1, e2);
  BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
  BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(eigenvalue4) {
  double unit[4] = {0.0, 1.0, -1.0, 0.0};
  double e1, e2;
  Matrix2x2::EigenValues(unit, e1, e2);
  if (e1 < e2) std::swap(e1, e2);
  BOOST_CHECK(!std::isfinite(e1));
  BOOST_CHECK(!std::isfinite(e2));
}

BOOST_AUTO_TEST_CASE(eigenvector2) {
  double unit[4] = {0.0, 1.0, -2.0, -3.0};
  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
  if (e1 < e2) {
    std::swap(e1, e2);
    std::swap(vec1, vec2);
  }
  BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
  BOOST_CHECK_CLOSE(vec1[0] / vec1[1], -1.0, 1e-6);  // vec1 = c [-1, 1]
  BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
  BOOST_CHECK_CLOSE(vec2[0] / vec2[1], -0.5, 1e-6);  // vec2 = c [-1, 2]
}

BOOST_AUTO_TEST_CASE(eigenvector3) {
  double unit[4] = {0.0, -2.0, 1.0, -3.0};
  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
  if (e1 < e2) {
    std::swap(e1, e2);
    std::swap(vec1, vec2);
  }
  BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
  BOOST_CHECK_CLOSE(vec1[0] / vec1[1], 2.0, 1e-6);  // vec1 = c [2, 1]
  BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
  BOOST_CHECK_CLOSE(vec2[0] / vec2[1], 1.0, 1e-6);  // vec2 = c [1, 1]
}

BOOST_AUTO_TEST_CASE(eigenvector4) {
  double unit[4] = {1.0, 2.0, 3.0, -4.0};
  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
  if (e1 < e2) {
    std::swap(e1, e2);
    std::swap(vec1, vec2);
  }
  BOOST_CHECK_CLOSE(e1, 2.0, 1e-6);
  BOOST_CHECK_CLOSE(vec1[0] / vec1[1], 2.0, 1e-6);  // vec1 = c [2, 1]
  BOOST_CHECK_CLOSE(e2, -5.0, 1e-6);
  BOOST_CHECK_CLOSE(vec2[1] / vec2[0], -3.0, 1e-6);  // vec2 = c [-2, 6]
}

BOOST_AUTO_TEST_CASE(eigenvector5) {
  double m[4] = {1.0, 0.0, 0.0, 0.5};
  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);
  if (e1 < e2) {
    std::swap(e1, e2);
    std::swap(vec1, vec2);
  }
  BOOST_CHECK_CLOSE(e1, 1.0, 1e-6);
  BOOST_CHECK_CLOSE(vec1[1] / vec1[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(e2, 0.5, 1e-6);
  BOOST_CHECK_CLOSE(vec2[0] / vec2[1], 0.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(cholesky_real) {
  std::complex<double> matrixA[4] = {1., 2., 2., 13.};
  std::complex<double> matrixB[4] = {1., 2., 2., 13.};
  const std::complex<double> answer[4] = {1., 0., 2., 3.};

  BOOST_CHECK(Matrix2x2::Cholesky(matrixA));
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(matrixA[i].real(), answer[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(matrixA[i].imag(), answer[i].imag(), 1e-6);
  }

  Matrix2x2::UncheckedCholesky(matrixB);
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(matrixB[i].real(), answer[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(matrixB[i].imag(), answer[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(cholesky_complex) {
  std::complex<double> matrixA[4] = {{1., 0.}, {2., -5.}, {2., 5.}, {38., 0.}};
  std::complex<double> matrixB[4] = {{1., 0.}, {2., -5.}, {2., 5.}, {38., 0.}};
  std::complex<double> answer[4] = {{1., 0.}, {0., 0.}, {2., 5.}, {3., 0.}};
  BOOST_CHECK(Matrix2x2::CheckedCholesky(matrixA));
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(matrixA[i].real(), answer[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(matrixA[i].imag(), answer[i].imag(), 1e-6);
  }

  Matrix2x2::UncheckedCholesky(matrixB);
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(matrixB[i].real(), answer[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(matrixB[i].imag(), answer[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(cholesky_not_positive) {
  std::complex<double> diag_not_positive[4] = {
      {0., 0.}, {0., 0.}, {0., 0.}, {1., 0.}};  // diagonal not positive
  BOOST_CHECK(!Matrix2x2::CheckedCholesky(diag_not_positive));
  std::complex<double> diag_not_real[4] = {
      {1., 0.}, {0., 0.}, {0., 0.}, {1., 1.}};  // diagonal not real
  BOOST_CHECK(!Matrix2x2::CheckedCholesky(diag_not_real));
  std::complex<double> not_hermitian[4] = {
      {1., 0.}, {1., 0.}, {2., 0.}, {1., 0.}};  // not hermitian
  BOOST_CHECK(!Matrix2x2::CheckedCholesky(not_hermitian));
}

BOOST_AUTO_TEST_CASE(eigen_value_and_vectors_real) {
  double m[] = {4.0, 1.0, 0.0, 4.0};

  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);

  BOOST_CHECK_CLOSE(e1, 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e2, 4.0, 1e-5);

  BOOST_CHECK_CLOSE(vec1[0], -1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1], 0.0, 1e-5);

  BOOST_CHECK_CLOSE(vec2[0], -1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1], 0.0, 1e-5);

  // Of course this is no longer necessary when the above checks
  // are already done, but e.g. signs are actually ambiguous in
  // above equations, so this is the real equation that should hold:
  BOOST_CHECK_CLOSE(m[0] * vec1[0] + m[1] * vec1[1], e1 * vec1[0], 1e-5);
  BOOST_CHECK_CLOSE(m[2] * vec1[0] + m[3] * vec1[1], e1 * vec1[1], 1e-5);
}

BOOST_AUTO_TEST_CASE(eigen_value_and_vectors_complex) {
  std::complex<double> m[] = {
      std::complex<double>(4.0, 1.0), std::complex<double>(1.0, 0.0),
      std::complex<double>(0.0, 0.0), std::complex<double>(4.0, 1.0)};

  std::complex<double> e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);

  BOOST_CHECK_CLOSE(e1.real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e1.imag(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.imag(), 1.0, 1e-5);

  BOOST_CHECK_CLOSE(vec1[0].real(), -1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[0].imag(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1].real(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1].imag(), 0.0, 1e-5);

  BOOST_CHECK_CLOSE(vec2[0].real(), -1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[0].imag(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1].real(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1].imag(), 0.0, 1e-5);

  BOOST_CHECK_LT(std::abs(m[0] * vec1[0] + m[1] * vec1[1] - e1 * vec1[0]),
                 1e-5);
  BOOST_CHECK_LT(std::abs(m[2] * vec1[0] + m[3] * vec1[1] - e1 * vec1[1]),
                 1e-5);
}

BOOST_AUTO_TEST_CASE(eigen_value_order_real) {
  // Test a specific case for which the eigen vector order
  // is "ambiguous". vec1 should always be associated with
  // e1, and vec2 with e2.
  // vec1 = { 0 , 1 }
  // vec2 = { 1 , 0 }
  // e1 = 4, e2 = 3
  // m {0, 1}^T = {0, 4} and m {1, 0}^T = {3, 0}
  // m = [ 3 0 ; 0 4 ]
  double m[] = {3.0, 0.0, 0.0, 4.0};

  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);

  BOOST_CHECK_CLOSE(e1, 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e2, 3.0, 1e-5);

  BOOST_CHECK_CLOSE(vec1[0], 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1], 1.0, 1e-5);

  BOOST_CHECK_CLOSE(vec2[0], 1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1], 0.0, 1e-5);

  BOOST_CHECK_CLOSE(m[0] * vec1[0] + m[1] * vec1[1], e1 * vec1[0], 1e-5);
  BOOST_CHECK_CLOSE(m[2] * vec1[0] + m[3] * vec1[1], e1 * vec1[1], 1e-5);
}

BOOST_AUTO_TEST_CASE(eigen_value_order1_complex) {
  // Test a specific case for which the eigen vector order
  // is "ambiguous". vec1 should always be associated with
  // e1, and vec2 with e2.
  // vec1 = { 0 , 1 }
  // vec2 = { 1 , 0 }
  // e1 = 4 + i, e2 = 3 + i
  // m {0, 1}^T = {0, 4+i} and m {1, 0}^T = {3+i, 0}
  // m = [ 3+i 0 ; 0 4+i ]
  std::complex<double> m[] = {
      std::complex<double>(3.0, 1.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0), std::complex<double>(4.0, 1.0)};

  std::complex<double> e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);

  BOOST_CHECK_CLOSE(e1.real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e1.imag(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.real(), 3.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.imag(), 1.0, 1e-5);

  BOOST_CHECK_CLOSE(vec1[0].real(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1].real(), 1.0, 1e-5);

  BOOST_CHECK_CLOSE(vec2[0].real(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1].real(), 0.0, 1e-5);

  const std::complex<double> lhs1 = m[0] * vec1[0] + m[1] * vec1[1],
                             rhs1 = e1 * vec1[0],
                             lhs2 = m[2] * vec1[0] + m[3] * vec1[1],
                             rhs2 = e1 * vec1[1];
  BOOST_CHECK_LT(std::abs(lhs1 - rhs1), 1e-5);
  BOOST_CHECK_LT(std::abs(lhs2 - rhs2), 1e-5);
}

BOOST_AUTO_TEST_CASE(eigen_value_order2_complex) {
  // vec1 = { 1 , 0 }
  // vec2 = { 0 , 1 }
  // e1 = 4 + i, e2 = 3 + i
  // m {1, 0}^T = {4+i, 0} and m {0, 1}^T = {0, 3+i}
  // m = [ 4+i 0 ; 0 3+i ]
  std::complex<double> m[] = {
      std::complex<double>(4.0, 1.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0), std::complex<double>(3.0, 1.0)};

  std::complex<double> e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(m, e1, e2, vec1, vec2);

  BOOST_CHECK_CLOSE(e1.real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE(e1.imag(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.real(), 3.0, 1e-5);
  BOOST_CHECK_CLOSE(e2.imag(), 1.0, 1e-5);

  BOOST_CHECK_CLOSE(vec1[0].real(), 1.0, 1e-5);
  BOOST_CHECK_CLOSE(vec1[1].real(), 0.0, 1e-5);

  BOOST_CHECK_CLOSE(vec2[0].real(), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(vec2[1].real(), 1.0, 1e-5);

  const std::complex<double> lhs1 = m[0] * vec1[0] + m[1] * vec1[1],
                             rhs1 = e1 * vec1[0],
                             lhs2 = m[2] * vec1[0] + m[3] * vec1[1],
                             rhs2 = e1 * vec1[1];
  BOOST_CHECK_LT(std::abs(lhs1 - rhs1), 1e-5);
  BOOST_CHECK_LT(std::abs(lhs2 - rhs2), 1e-5);
}

BOOST_AUTO_TEST_CASE(evdecomposition) {
  MC2x2 a(1, 2, 3, 4), b(5, 6, 7, 8);
  MC2x2 jones = a.MultiplyHerm(b) + b.MultiplyHerm(a);
  MC2x2 r = jones;
  r *= r.HermTranspose();
  std::complex<double> e1, e2, vec1[2], vec2[2];
  std::complex<double> r_data[4];
  r.AssignTo(r_data);
  Matrix2x2::EigenValuesAndVectors(r_data, e1, e2, vec1, vec2);
  double v1norm = std::norm(vec1[0]) + std::norm(vec1[1]);
  vec1[0] /= sqrt(v1norm);
  vec1[1] /= sqrt(v1norm);
  double v2norm = std::norm(vec2[0]) + std::norm(vec2[1]);
  vec2[0] /= sqrt(v2norm);
  vec2[1] /= sqrt(v2norm);

  MC2x2 u(vec1[0], vec2[0], vec1[1], vec2[1]), e(e1, 0, 0, e2);
  MC2x2 res = u.Multiply(e).MultiplyHerm(u);
  for (size_t i = 0; i != 4; ++i)
    BOOST_CHECK_CLOSE(res.Get(i).real(), r.Get(i).real(), 1e-6);

  MC2x2 decomposed = r.DecomposeHermitianEigenvalue();
  decomposed *= decomposed.HermTranspose();
  for (size_t i = 0; i != 4; ++i)
    BOOST_CHECK_CLOSE(decomposed.Get(i).real(), r.Get(i).real(), 1e-6);
}

BOOST_AUTO_TEST_CASE(herm_transpose) {
  const std::complex<double> a(1, 2);
  const std::complex<double> b(3, 4);
  const std::complex<double> c(5, 6);
  const std::complex<double> d(7, 8);
  const MC2x2 m(a, b, c, d);
  MC2x2 result = m.HermTranspose();
  CheckClose(result, MC2x2{{a.real(), -a.imag()},
                           {c.real(), -c.imag()},
                           {b.real(), -b.imag()},
                           {d.real(), -d.imag()}});
  result -= HermTranspose(m);
  for (size_t i = 0; i != 4; ++i)
    BOOST_CHECK_LT(std::norm(result.Get(i)), 1e-6);
}

BOOST_AUTO_TEST_CASE(conjugate, *boost::unit_test::tolerance(1e8)) {
  const std::complex<double> a(1, 2);
  const std::complex<double> b(3, 4);
  const std::complex<double> c(5, 6);
  const std::complex<double> d(7, 8);

  const MC2x2 m(a, b, c, d);
  const MC2x2 m_conj = m.Conjugate();
  BOOST_TEST(m_conj.Get(0) == std::conj(a));
  BOOST_TEST(m_conj.Get(1) == std::conj(b));
  BOOST_TEST(m_conj.Get(2) == std::conj(c));
  BOOST_TEST(m_conj.Get(3) == std::conj(d));
}

template <typename Num, typename Matrix>
void TestDoubleDot() {
  const std::complex<Num> a(1, 2);
  const std::complex<Num> b(3, 4);
  const std::complex<Num> c(5, 6);
  const std::complex<Num> d(7, 8);
  const Matrix m(a, b, c, d);

  // Double contraction with conjugate of itself should equal the matrix norm
  const std::complex<Num> result0 = m.DoubleDot(m.Conjugate());
  BOOST_CHECK_CLOSE(result0.real(), Norm(m), 1e-8);
  BOOST_CHECK_CLOSE(result0.imag(), 0.0, 1e-8);

  const std::complex<Num> result1 = m.DoubleDot(m);
  const std::complex<Num> result_ref = a * a + b * b + c * c + d * d;
  BOOST_CHECK_CLOSE(result1.real(), result_ref.real(), 1e-8);
  BOOST_CHECK_CLOSE(result1.imag(), result_ref.imag(), 1e-8);
}

BOOST_AUTO_TEST_CASE(double_dot) {
  TestDoubleDot<float, MC2x2F>();
  TestDoubleDot<double, MC2x2>();
}

BOOST_AUTO_TEST_CASE(trace) {
  const std::complex<double> a(1, 2);
  const std::complex<double> b(3, 4);
  const std::complex<double> c(5, 6);
  const std::complex<double> d(7, 8);
  const MC2x2 m(a, b, c, d);
  BOOST_CHECK_CLOSE(Trace(m).real(), (a + d).real(), 1e-6);
  BOOST_CHECK_CLOSE(Trace(m).imag(), (a + d).imag(), 1e-6);
  BOOST_CHECK_CLOSE((Trace(m) * 0.0).real(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE((Trace(m) * 0.0).imag(), 0.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(norm) {
  const std::complex<double> a(1, 2);
  const std::complex<double> b(3, 4);
  const std::complex<double> c(5, 6);
  const std::complex<double> d(7, 8);
  const MC2x2 m(a, b, c, d);
  double norm_result =
      1 * 1 + 2 * 2 + 3 * 3 + 4 * 4 + 5 * 5 + 6 * 6 + 7 * 7 + 8 * 8;
  BOOST_CHECK_CLOSE(Norm(m), norm_result, 1e-6);
  BOOST_CHECK_CLOSE(Norm(m * std::complex<double>(0.0, 0.0)), 0.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(sum_of_absolute) {
  const std::complex<double> a(3, 4);
  const std::complex<double> b(-7, 0);
  const std::complex<double> c(0, -8);
  const std::complex<double> d(6, 8);
  const MC2x2 m(a, b, c, d);
  double norm_result = 5 + 7 + 8 + 10;
  BOOST_CHECK_CLOSE(SumOfAbsolute(m), norm_result, 1e-6);
  BOOST_CHECK_CLOSE(SumOfAbsolute(m * std::complex<double>(0.0, 0.0)), 0.0,
                    1e-6);
}

BOOST_AUTO_TEST_CASE(element_wise_product_double) {
  const MC2x2 a({0.0, 1.0}, {2.0, 3.0}, {4.0, 5.0}, {6.0, 7.0});
  const MC2x2 b({10.0, 20.0}, {30.0, 40.0}, {50.0, 60.0}, {70.0, 80.0});
  const MC2x2 result = ElementProduct(a, b);
  BOOST_CHECK_CLOSE(result.Get(0).real(), -20.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(0).imag(), 10.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(1).real(), -60.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(1).imag(), 170.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(2).real(), -100.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(2).imag(), 490.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(3).real(), -140.0, 1e-9);
  BOOST_CHECK_CLOSE(result.Get(3).imag(), 970.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(element_wise_product_float) {
  const MC2x2F a({0.0f, 1.0f}, {2.0f, 3.0f}, {4.0f, 5.0f}, {6.0f, 7.0f});
  const MC2x2F b({10.0f, 20.0f}, {30.0f, 40.0f}, {50.0f, 60.0f},
                 {70.0f, 80.0f});
  const MC2x2F result = ElementProduct(a, b);
  BOOST_CHECK_CLOSE(result.Get(0).real(), -20.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(0).imag(), 10.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(1).real(), -60.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(1).imag(), 170.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(2).real(), -100.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(2).imag(), 490.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(3).real(), -140.0f, 1e-6);
  BOOST_CHECK_CLOSE(result.Get(3).imag(), 970.0f, 1e-6);
}

BOOST_AUTO_TEST_CASE(indexing_single) {
  const MC2x2F a({0.0f, 1.0f}, {2.0f, 3.0f}, {4.0f, 5.0f}, {6.0f, 7.0f});
  MC2x2F b = MC2x2F::Zero();
  b.Set(0, a.Get(0));
  b.Set(1, a.Get(1));
  b.Set(2, a.Get(2));
  BOOST_CHECK_NE(a, b);
  b.Set(3, a.Get(3));
  BOOST_CHECK_EQUAL(a, b);
}

BOOST_AUTO_TEST_CASE(indexing_double) {
  const MC2x2 a({0.0, 1.0}, {2.0, 3.0}, {4.0, 5.0}, {6.0, 7.0});
  MC2x2 b = MC2x2::Zero();
  b.Set(0, a.Get(0));
  b.Set(1, a.Get(1));
  b.Set(2, a.Get(2));
  BOOST_CHECK_NE(a, b);
  b.Set(3, a.Get(3));
  BOOST_CHECK_EQUAL(a, b);
}

#ifndef USE_AVX_MATRIX
BOOST_AUTO_TEST_CASE(index_real) {
  const MC2x2F m({4.0f, 5.0f}, {6.0f, 7.0f}, {8.0f, 9.0f}, {10.0f, 11.0f});
  for (size_t i = 0; i != 8; ++i) {
    BOOST_CHECK_CLOSE(m.IndexReal(i), i + 4, 1e-6);
  }
}
#endif

BOOST_AUTO_TEST_CASE(dubious_float_cast) {
  const MC2x2F const_matrix({0.0f, 1.0f}, {2.0f, 3.0f}, {4.0f, 5.0f},
                            {6.0f, 7.0f});
  MC2x2F writable_matrix(const_matrix);
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_EQUAL(aocommon::DubiousComplexPointerCast(const_matrix)[i],
                      std::complex<float>(i * 2, i * 2 + 1));

    aocommon::DubiousComplexPointerCast(writable_matrix)[i] =
        std::complex<float>(8.0f, 9.0f);
    BOOST_CHECK_EQUAL(aocommon::DubiousComplexPointerCast(writable_matrix)[i],
                      std::complex<float>(8.0f, 9.0f));
  }
}

BOOST_AUTO_TEST_CASE(dubious_double_cast) {
  const MC2x2 const_matrix({0.0, 1.0}, {2.0, 3.0}, {4.0, 5.0}, {6.0, 7.0});
  MC2x2 writable_matrix(const_matrix);
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_EQUAL(aocommon::DubiousDComplexPointerCast(const_matrix)[i],
                      std::complex<double>(i * 2, i * 2 + 1));

    MC2x2 writable_matrix({0.0, 1.0}, {2.0, 3.0}, {4.0, 5.0}, {6.0, 7.0});
    aocommon::DubiousDComplexPointerCast(writable_matrix)[i] =
        std::complex<double>(8.0, 9.0);
    BOOST_CHECK_EQUAL(aocommon::DubiousDComplexPointerCast(writable_matrix)[i],
                      std::complex<double>(8.0, 9.0));
  }
}

BOOST_AUTO_TEST_SUITE_END()
