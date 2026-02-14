#ifndef TEST_MATRIX4X4_H
#define TEST_MATRIX4X4_H

#include <boost/test/unit_test.hpp>

#include <aocommon/matrix4x4.h>

using aocommon::Matrix4x4;
using aocommon::MC2x2;
using aocommon::MC4x4;
using aocommon::Vector4;

BOOST_AUTO_TEST_SUITE(matrix4x4)

namespace {
constexpr MC4x4 GetExampleMatrix() {
  return MC4x4{{1.0, -2.0},    {3.0, 4.0},     {-5.0, 6.0},   {-7.0, -8.0},
               {9.0, -10.0},   {11.0, 12.0},   {13.0, -14.0}, {-15.0, -16.0},
               {-17.0, -18.0}, {-19.0, 20.0},  {21.0, 22.0},  {23.0, -24.0},
               {25.0, -26.0},  {-27.0, -28.0}, {29.0, 30.0},  {-31.0, 32.0}};
}

constexpr MC4x4 GetExampleMatrixHermTransposed() {
  return MC4x4{{1.0, 2.0},   {9.0, 10.0},   {-17.0, 18.0},  {25.0, 26.0},
               {3.0, -4.0},  {11.0, -12.0}, {-19.0, -20.0}, {-27.0, 28.0},
               {-5.0, -6.0}, {13.0, 14.0},  {21.0, -22.0},  {29.0, -30.0},
               {-7.0, 8.0},  {-15.0, 16.0}, {23.0, 24.0},   {-31.0, -32.0}};
}
}  // namespace

static void CheckMatrix(const Matrix4x4& result, const Matrix4x4& groundtruth) {
  for (size_t i = 0; i != 16; ++i) {
    BOOST_CHECK_CLOSE(result[i].real(), groundtruth[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(result[i].imag(), groundtruth[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(construction) {
  CheckMatrix(MC4x4(), MC4x4{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
}

BOOST_AUTO_TEST_CASE(unit) {
  MC4x4 unit = MC4x4::Unit();
  MC4x4 ref{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  CheckMatrix(unit, ref);
}

BOOST_AUTO_TEST_CASE(multiplication) {
  const MC4x4 unit = MC4x4::Unit();
  CheckMatrix(unit * unit, unit);

  const MC4x4 a{1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,
                9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
  const MC4x4 b = a * 3.0;
  const MC4x4 result{270, 300,  330,  360,  606,  684,  762,  840,
                     942, 1068, 1194, 1320, 1278, 1452, 1626, 1800};
  CheckMatrix(a * b, result);

  const MC4x4 result2{810,  900,  990,  1080, 1818, 2052, 2286, 2520,
                      2826, 3204, 3582, 3960, 3834, 4356, 4878, 5400};
  CheckMatrix(b.Square(), result2);
}

BOOST_AUTO_TEST_CASE(herm_transpose) {
  CheckMatrix(MC4x4::Zero(), MC4x4::Zero().HermTranspose());
  CheckMatrix(MC4x4::Unit() * 3, (MC4x4::Unit() * 3).HermTranspose());
  CheckMatrix(GetExampleMatrix().HermTranspose(),
              GetExampleMatrixHermTransposed());
  CheckMatrix(GetExampleMatrix(),
              GetExampleMatrixHermTransposed().HermTranspose());
  CheckMatrix(GetExampleMatrix().HermTranspose().HermTranspose(),
              GetExampleMatrix());
}

BOOST_AUTO_TEST_CASE(inversion) {
  MC4x4 m1(MC4x4::Unit());
  BOOST_CHECK(m1.Invert());
  CheckMatrix(m1, MC4x4::Unit());

  MC4x4 m2(MC4x4::Unit() * 2);
  BOOST_CHECK(m2.Invert());
  CheckMatrix(m2, MC4x4::Unit() * 0.5);
  BOOST_CHECK(m2.Invert());
  CheckMatrix(m2, MC4x4::Unit() * 2.0);

  MC4x4 m3;
  BOOST_CHECK(!m3.Invert());
}

BOOST_AUTO_TEST_CASE(hermitian_square) {
  CheckMatrix(MC4x4::Zero().HermitianSquare().ToMatrix(), MC4x4::Zero());
  CheckMatrix(MC4x4::Unit().HermitianSquare().ToMatrix(), MC4x4::Unit());

  constexpr MC4x4 m = GetExampleMatrix();
  CheckMatrix(m.HermitianSquare().ToMatrix(), m.HermTranspose() * m);
}

static void checkKroneckerProduct(const MC2x2& a, const MC2x2& x,
                                  const MC2x2& b) {
  Vector4 ref = a.Multiply(x).MultiplyHerm(b).Vec();
  MC4x4 product = MC4x4::KroneckerProduct(b.HermTranspose().Transpose(), a);
  Vector4 v = product * x.Vec();
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(v[i].real(), ref[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(v[i].imag(), ref[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(kronecker_product) {
  checkKroneckerProduct(MC2x2::Unity(), MC2x2::Unity(), MC2x2::Unity());

  MC2x2 a1{1.0, 2.0, 2.0, 4.0}, x1(MC2x2::Unity()), b1{1.0, 2.0, 2.0, 4.0};
  checkKroneckerProduct(a1, x1, b1);

  MC2x2 a2{0.0, 1.0, 2.0, 3.0}, x2(MC2x2::Unity()), b2{0.0, 1.0, 2.0, 3.0};
  checkKroneckerProduct(a2, x2, b2);

  MC2x2 a3{0.0, 1.0, 2.0, 3.0}, x3{0.0, 1.0, 2.0, 3.0}, b3{0.0, 1.0, 2.0, 3.0};
  checkKroneckerProduct(a3, x3, b3);

  std::complex<double> x(8, 2), y(6, 3);
  MC2x2 a4{0.0, 1.0 * y, 2.0 * x, 3.0 * y},
      x4{1.0 * y, 2.0 * x, 3.0 * x, 4.0 * y},
      b4{1.0 * x, 2.0 * x, 3.0 * x, 4.0 * y};
  checkKroneckerProduct(a4, x4, b4);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
