#include <boost/test/unit_test.hpp>

#include <aocommon/hmatrix4x4.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix4x4.h>

using aocommon::HMC4x4;
using aocommon::MC2x2;
using aocommon::MC4x4;
using aocommon::Vector4;

#define CHECK_CLOSE_MESSAGE(VAL, REF, MSG)         \
  BOOST_CHECK_MESSAGE(std::fabs(VAL - REF) < 1e-6, \
                      MSG << " is " << VAL << ", should be " << REF);

namespace {
HMC4x4 GetExampleMatrix() {
  const std::complex<double> j(0, 1);
  return HMC4x4{
      1.0,  2.0 + 3.0 * j,   4.0 - 5.0 * j,   6.0 + 7.0 * j,   2.0 - 3.0 * j,
      8.0,  9.0 + 10.0 * j,  11.0 - 12.0 * j, 4.0 + 5.0 * j,   9.0 - 10.0 * j,
      13.0, 14.0 + 15.0 * j, 6.0 - 7.0 * j,   11.0 + 12.0 * j, 14.0 - 15.0 * j,
      16.0};
}
}  // namespace

BOOST_AUTO_TEST_SUITE(hmatrix4x4)

template <typename Matrix>
static void CheckMatrix(const Matrix& result, const Matrix& groundtruth) {
  for (size_t i = 0; i != 16; ++i) {
    BOOST_CHECK_CLOSE(result[i].real(), groundtruth[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(result[i].imag(), groundtruth[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(zero) {
  constexpr HMC4x4 zero = HMC4x4::Zero();
  constexpr HMC4x4 ref{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  CheckMatrix(zero, ref);
  CheckMatrix(zero.ToMatrix(), MC4x4::Zero());
}

BOOST_AUTO_TEST_CASE(unit) {
  constexpr HMC4x4 unit = HMC4x4::Unit();
  constexpr HMC4x4 ref{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  CheckMatrix(unit, ref);
  CheckMatrix(unit.ToMatrix(), MC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(equals) {
  static_assert(HMC4x4::Unit() == HMC4x4::Unit());
  static_assert(!(HMC4x4::Zero() == HMC4x4::Unit()));
  static_assert(HMC4x4::Zero() == HMC4x4::Zero());
  static_assert(!(HMC4x4::Unit() == HMC4x4::Zero()));
  for (size_t i = 0; i != 16; ++i) {
    HMC4x4 m;
    m.Data(i) = 1e-8;
    BOOST_CHECK(m == m);
    BOOST_CHECK(!(m == HMC4x4::Zero()));
  }
}

BOOST_AUTO_TEST_CASE(unequal) {
  static_assert(!(HMC4x4::Unit() != HMC4x4::Unit()));
  static_assert(HMC4x4::Zero() != HMC4x4::Unit());
  static_assert(!(HMC4x4::Zero() != HMC4x4::Zero()));
  static_assert(HMC4x4::Unit() != HMC4x4::Zero());
  for (size_t i = 0; i != 16; ++i) {
    HMC4x4 m;
    m.Data(i) = 1e-8;
    BOOST_CHECK(!(m != m));
    BOOST_CHECK(m != HMC4x4::Zero());
  }
}

BOOST_AUTO_TEST_CASE(buffer) {
  HMC4x4 ref = HMC4x4::Unit();
  std::array<double, 16> buffer;
  for (size_t i = 0; i != buffer.size(); ++i) {
    buffer[i] = ref.Data(i);
  }
  HMC4x4 unit(buffer.data());
  CheckMatrix(unit, ref);
  CheckMatrix(unit.ToMatrix(), MC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(diagonal_values) {
  HMC4x4 matrix = HMC4x4::Unit();
  std::array<double, 4> diagonal = matrix.DiagonalValues();

  const std::array<double, 4> unit{1.0, 1.0, 1.0, 1.0};
  BOOST_CHECK_EQUAL_COLLECTIONS(diagonal.begin(), diagonal.end(), unit.begin(),
                                unit.end());

  // 0,3,8,15 come from the ordering of the data inside the HMC4x4 matrix
  // (see help for HCM4x4::Data()).
  matrix.Data(0) = 12.0;
  matrix.Data(3) = 13.0;
  matrix.Data(8) = 14.0;
  matrix.Data(15) = 15.0;
  diagonal = matrix.DiagonalValues();
  const std::array<double, 4> ref{12.0, 13.0, 14.0, 15.0};
  BOOST_CHECK_EQUAL_COLLECTIONS(diagonal.begin(), diagonal.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(addition) {
  const HMC4x4 unit = HMC4x4::Unit();
  CheckMatrix(unit + HMC4x4::Zero(), unit);
  CheckMatrix(unit + unit, unit * 2);
  HMC4x4 a = GetExampleMatrix();
  MC4x4 b;
  MC4x4 reference;
  for (size_t col = 0; col != 4; ++col) {
    for (size_t row = col; row != 4; ++row) {
      if (col == row)
        b[col + row * 4] = std::complex<double>(col + row, 0.0);
      else
        b[col + row * 4] = std::complex<double>(col + row, col);
      reference[col + row * 4] = b[col + row * 4] + a[col + row * 4];
    }
  }
  CheckMatrix(a + HMC4x4(b), HMC4x4(reference));
}

BOOST_AUTO_TEST_CASE(inversion) {
  HMC4x4 m1(HMC4x4::Unit());
  BOOST_CHECK(m1.Invert());
  CheckMatrix(m1, HMC4x4::Unit());

  HMC4x4 m2(HMC4x4::Unit() * 2);
  BOOST_CHECK(m2.Invert());
  CheckMatrix(m2, HMC4x4::Unit() * 0.5);
  BOOST_CHECK(m2.Invert());
  CheckMatrix(m2, HMC4x4::Unit() * 2.0);

  HMC4x4 m3;
  BOOST_CHECK(!m3.Invert());
}

BOOST_AUTO_TEST_CASE(from_data) {
  // Identity matrix
  HMC4x4 m = HMC4x4::FromData({1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
  CheckMatrix(m, HMC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(indexing1) {
  HMC4x4 m{1.0, 2.0, 4.0, 7.0, 2.0, 3.0, 5.0, 8.0,
           4.0, 5.0, 6.0, 9.0, 7.0, 8.0, 9.0, 10.0};
  const double vals[16] = {1.0, 2.0, 4.0, 7.0, 2.0, 3.0, 5.0, 8.0,
                           4.0, 5.0, 6.0, 9.0, 7.0, 8.0, 9.0, 10.0};
  for (size_t i = 0; i != 16; ++i) {
    BOOST_CHECK_CLOSE(m[i].real(), vals[i], 1e-6);
    BOOST_CHECK_CLOSE(m[i].imag(), 0.0, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(indexing2) {
  std::complex<double> j(0.0, 1.0);
  HMC4x4 m{1.0, 2.0 - j, 4.0, 7.0, 2.0 + j, 3.0,     5.0, 8.0 + j,
           4.0, 5.0,     6.0, 9.0, 7.0,     8.0 - j, 9.0, 10.0};
  const std::complex<double> vals[16] = {1.0, 2.0 - j, 4.0, 7.0, 2.0 + j, 3.0,
                                         5.0, 8.0 + j, 4.0, 5.0, 6.0,     9.0,
                                         7.0, 8.0 - j, 9.0, 10.0};
  for (size_t i = 0; i != 16; ++i) {
    CHECK_CLOSE_MESSAGE(m[i].real(), vals[i].real(), "Real element " << i);
    CHECK_CLOSE_MESSAGE(m[i].imag(), vals[i].imag(), "Imag element " << i);
  }
}

BOOST_AUTO_TEST_CASE(scalar_product) {
  HMC4x4 ref = HMC4x4{2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
                      0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0};
  // scalar multiplication
  CheckMatrix(HMC4x4::Unit() * 2.0, ref);
  // scalar multiplication-assignment
  HMC4x4 m = HMC4x4::Unit();
  m *= 2.0;
  CheckMatrix(m, ref);
}

BOOST_AUTO_TEST_CASE(scalar_division) {
  HMC4x4 m = HMC4x4{2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
                    0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0};
  // scalar division-assignment
  m /= 2.0;
  CheckMatrix(m, HMC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(product_with_vector4) {
  Vector4 v1(2.0, 2.0, 2.0, 2.0);
  Vector4 res = HMC4x4::Unit() * v1;
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(res[i].real(), 2.0, 1e-6);
    BOOST_CHECK_CLOSE(res[i].imag(), 0.0, 1e-6);
  }
  Vector4 v2(std::complex<double>(2.0, 3.0), std::complex<double>(4.0, 5.0),
             std::complex<double>(5.0, 6.0), std::complex<double>(7.0, 8.0));
  res = HMC4x4::Unit() * 0.5 * v2;
  Vector4 ref = MC4x4::Unit() * 0.5 * v2;
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(res[i].real(), ref[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(res[i].imag(), ref[i].imag(), 1e-6);
  }
  std::complex<double> j(0.0, 1.0);
  MC4x4 m{1.0,           2.0 + 1.0 * j, 3.0 + 2.0 * j, 4.0 + 3.0 * j,
          2.0 - 1.0 * j, 2.0,           3.0 + 2.0 * j, 4.0 + 2.0 * j,
          3.0 - 2.0 * j, 3.0 - 2.0 * j, 3.0,           4.0 - 3.0 * j,
          4.0 - 3.0 * j, 4.0 - 2.0 * j, 4.0 + 3.0 * j, 4.0};

  res = HMC4x4(m) * v2;
  ref = m * v2;
  for (size_t i = 0; i != 4; ++i) {
    CHECK_CLOSE_MESSAGE(res[i].real(), ref[i].real(), "Element " << i);
    CHECK_CLOSE_MESSAGE(res[i].imag(), ref[i].imag(), "Element " << i);
  }
}

static void checkKroneckerProduct(const MC2x2& a, const MC2x2& x,
                                  const MC2x2& b) {
  Vector4 ref = a.Multiply(x).MultiplyHerm(b).Vec();
  HMC4x4 product = HMC4x4::KroneckerProduct(b.HermTranspose().Transpose(), a);
  Vector4 v = product * x.Vec();
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_CLOSE(v[i].real(), ref[i].real(), 1e-6);
    BOOST_CHECK_CLOSE(v[i].imag(), ref[i].imag(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(kronecker_product_a) {
  checkKroneckerProduct(MC2x2::Unity(), MC2x2::Unity(), MC2x2::Unity());
}

BOOST_AUTO_TEST_CASE(kronecker_product_b) {
  MC2x2 a1{1.0, 2.0, 2.0, 4.0}, x1(MC2x2::Unity()), b1{1.0, 2.0, 2.0, 4.0};
  checkKroneckerProduct(a1, x1, b1);
}

BOOST_AUTO_TEST_CASE(kronecker_product_c) {
  MC2x2 a3{0.0, 1.0, 1.0, 3.0}, x3{0.0, 1.0, 2.0, 3.0}, b3{0.0, 1.0, 1.0, 3.0};
  checkKroneckerProduct(a3, x3, b3);
}

BOOST_AUTO_TEST_CASE(kronecker_product_d) {
  std::complex<double> x(8, 2), y(6, 3), xc = std::conj(x), yc = std::conj(y);
  MC2x2 a4{0.0, 2.0 * y, 2.0 * yc, 3.0}, x4{1.0, 2.0 * xc, 2.0 * x, 4.0},
      b4{1.0, 3.0 * x, 3.0 * xc, 4.0};
  checkKroneckerProduct(a4, x4, b4);
}

BOOST_AUTO_TEST_CASE(norm) {
  BOOST_CHECK_CLOSE((HMC4x4::Unit() * 2.0).Norm(), (MC4x4::Unit() * 2.0).Norm(),
                    1e-6);
  const HMC4x4 m = GetExampleMatrix();
  BOOST_CHECK_CLOSE(m.Norm(), m.ToMatrix().Norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(square) {
  CheckMatrix(HMC4x4::Unit().Square(), HMC4x4::Unit());
  CheckMatrix((HMC4x4::Unit() * 2.0).Square(), HMC4x4::Unit() * 4.0);
  const HMC4x4 m = GetExampleMatrix();
  const MC4x4 square = m.Square().ToMatrix();
  const std::complex<double> m00 = m[0];
  const std::complex<double> m01 = m[1];
  const std::complex<double> m02 = m[2];
  const std::complex<double> m03 = m[3];
  const std::complex<double> m10 = m[4];
  const std::complex<double> m11 = m[5];
  const std::complex<double> m12 = m[6];
  const std::complex<double> m13 = m[7];
  const std::complex<double> m20 = m[8];
  const std::complex<double> m21 = m[9];
  const std::complex<double> m22 = m[10];
  const std::complex<double> m23 = m[11];
  const std::complex<double> m30 = m[12];
  const std::complex<double> m31 = m[13];
  const std::complex<double> m32 = m[14];
  const std::complex<double> m33 = m[15];
  // This is written out quite verbosely because when using functions to make
  // things shorter, it is much harder to trace errors.
  std::complex r00 = (m00 * m00) + (m01 * m10) + (m02 * m20) + (m03 * m30);
  BOOST_CHECK_CLOSE(square[0].real(), r00.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[0].imag(), r00.imag(), 1e-6);
  std::complex r01 = (m00 * m01) + (m01 * m11) + (m02 * m21) + (m03 * m31);
  BOOST_CHECK_CLOSE(square[1].real(), r01.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[1].imag(), r01.imag(), 1e-6);
  std::complex r02 = (m00 * m02) + (m01 * m12) + (m02 * m22) + (m03 * m32);
  BOOST_CHECK_CLOSE(square[2].real(), r02.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[2].imag(), r02.imag(), 1e-6);
  std::complex r03 = (m00 * m03) + (m01 * m13) + (m02 * m23) + (m03 * m33);
  BOOST_CHECK_CLOSE(square[3].real(), r03.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[3].imag(), r03.imag(), 1e-6);

  std::complex r11 = (m10 * m01) + (m11 * m11) + (m12 * m21) + (m13 * m31);
  BOOST_CHECK_CLOSE(square[5].real(), r11.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[5].imag(), r11.imag(), 1e-6);
  std::complex r12 = (m10 * m02) + (m11 * m12) + (m12 * m22) + (m13 * m32);
  BOOST_CHECK_CLOSE(square[6].real(), r12.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[6].imag(), r12.imag(), 1e-6);
  std::complex r13 = (m10 * m03) + (m11 * m13) + (m12 * m23) + (m13 * m33);
  BOOST_CHECK_CLOSE(square[7].real(), r13.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[7].imag(), r13.imag(), 1e-6);

  std::complex r22 = (m20 * m02) + (m21 * m12) + (m22 * m22) + (m23 * m32);
  BOOST_CHECK_CLOSE(square[10].real(), r22.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[10].imag(), r22.imag(), 1e-6);
  std::complex r23 = (m20 * m03) + (m21 * m13) + (m22 * m23) + (m23 * m33);
  BOOST_CHECK_CLOSE(square[11].real(), r23.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[11].imag(), r23.imag(), 1e-6);

  std::complex r33 = (m30 * m03) + (m31 * m13) + (m32 * m23) + (m33 * m33);
  BOOST_CHECK_CLOSE(square[15].real(), r33.real(), 1e-6);
  BOOST_CHECK_CLOSE(square[15].imag(), r33.imag(), 1e-6);
}

BOOST_AUTO_TEST_CASE(serialize) {
  aocommon::SerialOStream o_stream;
  const HMC4x4 m1 = HMC4x4::Zero();
  const HMC4x4 m2 = GetExampleMatrix();
  m1.Serialize(o_stream);
  m2.Serialize(o_stream);
  aocommon::SerialIStream i_stream(std::move(o_stream));
  HMC4x4 m1result;
  HMC4x4 m2result;
  m1result.Unserialize(i_stream);
  CheckMatrix(m1result, m1);
  m2result.Unserialize(i_stream);
  CheckMatrix(m2result, m2);
}

BOOST_AUTO_TEST_SUITE_END()
