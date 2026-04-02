#ifndef AOCOMMON_HMATRIX_4X4_H_
#define AOCOMMON_HMATRIX_4X4_H_

#include <array>
#include <complex>
#include <string>
#include <sstream>
#include <stdexcept>

#include "aocommon/io/serialostream.h"
#include "aocommon/io/serialistream.h"
#include "aocommon/matrix2x2.h"

namespace aocommon {

class Matrix4x4;

/**
 * Class implements a Hermitian 4x4 matrix. Internally, the data is
 * stored as 16 doubles, rather than storing 16 complex doubles.
 */
class HMatrix4x4 {
 public:
  /**
   * A Hermitian matrix with zeros.
   */
  constexpr HMatrix4x4()
      : _data{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} {}

  /**
   * Construct a HMatrix4x4 object from a buffer with size 16
   */
  constexpr HMatrix4x4(const double buffer[16])
      : _data{buffer[0],  buffer[1],  buffer[2],  buffer[3],
              buffer[4],  buffer[5],  buffer[6],  buffer[7],
              buffer[8],  buffer[9],  buffer[10], buffer[11],
              buffer[12], buffer[13], buffer[14], buffer[15]} {}

  /**
   * Construct a new HMatrix4x4 object from the lower triangular entries in a
   * Matrix4x4 object. No error is thrown if the input is not Hermitian.
   */
  constexpr explicit HMatrix4x4(const aocommon::Matrix4x4& src);

  /**
   * Construct a new HMatrix4x4 object from an initializer list of
   * std::complex<double> having length 16. No error is thrown if
   * the input matrix is not Hermitian. Initialization is done from
   * the diagonal + lower triangle of the given matrix.
   *
   * @param list Initializer list of std::complex<double>, should have length 16
   */
  constexpr HMatrix4x4(std::initializer_list<std::complex<double>> list)
      :  // row 0
        _data{list.begin()[0].real(),
              // row 1
              list.begin()[4].real(), list.begin()[4].imag(),
              list.begin()[5].real(),
              // row 2
              list.begin()[8].real(), list.begin()[8].imag(),
              list.begin()[9].real(), list.begin()[9].imag(),
              list.begin()[10].real(),
              // row 3
              list.begin()[12].real(), list.begin()[12].imag(),
              list.begin()[13].real(), list.begin()[13].imag(),
              list.begin()[14].real(), list.begin()[14].imag(),
              list.begin()[15].real()} {}

  /**
   * Make a HMatrix4x4 object from an initializer list of (16) doubles
   *
   * @param list Initializer list of doubles with length 16
   * @return HMatrix4x4
   */
  static constexpr HMatrix4x4 FromData(std::initializer_list<double> list) {
    assert(list.size() == 16);
    HMatrix4x4 m;
    double* ptr = m._data.data();
    for (double e : list) {
      *ptr = e;
      ++ptr;
    }
    return m;
  }

  /**
   * Return HMatrix4x4 filled with zeros.
   */
  static constexpr HMatrix4x4 Zero() { return HMatrix4x4(); }

  /**
   * Return Hermitian 4x4 identity matrix.
   */
  static constexpr HMatrix4x4 Unit() {
    HMatrix4x4 unit;
    unit._data[0] = 1.0;
    unit._data[3] = 1.0;
    unit._data[8] = 1.0;
    unit._data[15] = 1.0;
    return unit;
  }

  constexpr bool operator==(const HMatrix4x4& rhs) const {
    // operator== of std::array is constexpr only since C++20
    for (size_t i = 0; i != 16; ++i)
      if (_data[i] != rhs._data[i]) return false;
    return true;
  }

  constexpr bool operator!=(const HMatrix4x4& rhs) const {
    return !(*this == rhs);
  }

  /**
   * Return the values that are on the diagonal. These are by definition
   * real, because the matrix is Hermitian.
   */
  constexpr std::array<double, 4> DiagonalValues() const {
    return {_data[0], _data[3], _data[8], _data[15]};
  }

  constexpr HMatrix4x4 operator+(const HMatrix4x4& rhs) const {
    HMatrix4x4 result;
    for (size_t i = 0; i != 16; ++i) result._data[i] = _data[i] + rhs._data[i];
    return result;
  }

  /**
   * Addition assignment operator
   */
  constexpr HMatrix4x4& operator+=(const HMatrix4x4& rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] += rhs._data[i];
    return *this;
  }

  /**
   * Scalar multiplication operator
   */
  constexpr HMatrix4x4 operator*(double rhs) const {
    HMatrix4x4 m;
    for (size_t i = 0; i != 16; ++i) m._data[i] = _data[i] * rhs;
    return m;
  }

  /**
   * Matrix-vector dot product
   */
  aocommon::Vector4 operator*(const aocommon::Vector4& rhs) const {
    aocommon::Vector4 v(_data[0] * rhs[0], (*this)[4] * rhs[0],
                        (*this)[8] * rhs[0], (*this)[12] * rhs[0]);
    v[0] += (*this)[1] * rhs[1];
    v[1] += _data[3] * rhs[1];
    v[2] += (*this)[1 + 8] * rhs[1];
    v[3] += (*this)[1 + 12] * rhs[1];

    v[0] += (*this)[2] * rhs[2];
    v[1] += (*this)[2 + 4] * rhs[2];
    v[2] += _data[8] * rhs[2];
    v[3] += (*this)[2 + 12] * rhs[2];

    v[0] += (*this)[3] * rhs[3];
    v[1] += (*this)[3 + 4] * rhs[3];
    v[2] += (*this)[3 + 8] * rhs[3];
    v[3] += _data[15] * rhs[3];

    return v;
  }

  /**
   * Scalar multiplication-assignment
   */
  constexpr HMatrix4x4& operator*=(double rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] *= rhs;
    return *this;
  }

  /**
   * Scalar division assignment
   */
  constexpr HMatrix4x4& operator/=(double rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] /= rhs;
    return *this;
  }

  /**
   * Invert matrix. Returns false if Hermitian not invertible.
   */
  bool Invert();

  /**
   * Indexing operator
   */
  constexpr std::complex<double> operator[](size_t i) const {
    constexpr size_t lookup[16] = {32, 17, 20, 25, 1, 35, 22, 27,
                                   4,  6,  40, 29, 9, 11, 13, 47};
    const size_t l = lookup[i];
    return ((l & 32) == 0)
               ? (((l & 16) == 0)
                      ? std::complex<double>(_data[l], _data[l + 1])
                      : std::complex<double>(_data[l & (~16)],
                                             -_data[(l & (~16)) + 1]))
               : (_data[l & (~32)]);
  }

  /**
   * Convert Hermitian to regular (complex-valued) 4x4 matrix
   */
  constexpr aocommon::Matrix4x4 ToMatrix() const;
  /**
   * "Entrywise" square of the L2 norm of the Hermitian
   */
  double Norm() const {
    return
        // diagonal
        _data[0] * _data[0] + _data[3] * _data[3] + _data[8] * _data[8] +
        _data[15] * _data[15] +
        // lower half x 2
        2.0 * (std::norm(ToComplex(1)) + std::norm(ToComplex(4)) +
               std::norm(ToComplex(6)) + std::norm(ToComplex(9)) +
               std::norm(ToComplex(11)) + std::norm(ToComplex(13)));
  }

  /**
   * Convert matrix to pretty string
   */
  std::string String() const {
    std::ostringstream str;
    for (size_t y = 0; y != 4; ++y) {
      for (size_t x = 0; x != 3; ++x) {
        str << (*this)[x + y * 4] << '\t';
      }
      str << (*this)[3 + y * 4] << '\n';
    }
    return str.str();
  }

  /**
   * Compute Hermitian matrix as the product of two 2x2 complex valued matrices.
   * Typical use case is to convert the product of two Jones matrices
   * into a Mueller matrix.
   */
  static HMatrix4x4 KroneckerProduct(const aocommon::MC2x2& hma,
                                     const aocommon::MC2x2& hmb) {
    using T = const std::complex<double>;
    using RT = const double;
    using std::conj;
    using std::norm;

    RT p00 = hma.Get(0).real();
    T p10 = hma.Get(2);
    RT p11 = hma.Get(3).real();

    RT q00 = hmb.Get(0).real();
    T q10 = hmb.Get(2);  // = conj(q01)
    RT q11 = hmb.Get(3).real();

    RT m00 = p00 * q00;
    T m10 = p00 * q10;
    RT m11 = p00 * q11;
    T m20 = p10 * q00;
    T m21 = p10 * conj(q10);
    RT m22 = p11 * q00;
    T m30 = p10 * q10;
    T m31 = p10 * q11;
    T m32 = p11 * q10;
    RT m33 = p11 * q11;

    return aocommon::HMatrix4x4::FromData(
        {m00, m10.real(), m10.imag(), m11, m20.real(), m20.imag(), m21.real(),
         m21.imag(), m22, m30.real(), m30.imag(), m31.real(), m31.imag(),
         m32.real(), m32.imag(), m33});
  }

  /**
   * Optimized version of @ref KroneckerProduct for diagonal matrix.
   * In this case most computation are optimized away and only the real values
   * from each matrix are required. So we take as input only the real portion
   * instead of the full complex numbers. Input can optimally be computed using
   * @ref RealHermitianSquare() instead of @ref HermitianSquare()
   */
  static HMatrix4x4 KroneckerProduct(const std::array<double, 2>& hma,
                                     const std::array<double, 2>& hmb) {
    return aocommon::HMatrix4x4::FromData({
        hma[0] * hmb[0],
        0.0f,
        0.0f,
        hma[0] * hmb[1],
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        hma[1] * hmb[0],
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        hma[1] * hmb[1],
    });
  }

  /**
   * Get underlying data by index, where 0 <= index <= 15. This indexing
   * is used since the data is internally stored as 16 doubles. The diagonal
   * is real, and only the lower (complex) half is stored, in
   * column-first order. The elements can therefore be indexed in
   * the following way:
   *  0
   *  1  3
   *  4  6 8
   *  9 11 13 15
   *
   * Note that "skipped indices" are the imaginary entries of the
   * Hermitian matrix
   */
  constexpr const double& Data(size_t index) const { return _data[index]; }
  constexpr double& Data(size_t index) { return _data[index]; }

  /**
   * Returns 0, 3, 8 and 15; the element indices of the diagonal
   * entries.
   */
  constexpr static std::array<size_t, 4> kDiagonalIndices{0, 3, 8, 15};

  /**
   * Returns A times A. Because the matrix is Hermitian, this is the
   * same as doing A^H times A. The result of this is also Hermitian,
   * proof: (AB)^H = B^H A^H, and given A and B are Hermitian,
   * (AB)^H = BA. Replacing B by A: (AA)^H = AA, and thus AA must be
   * Hermitian as well.
   */
  HMatrix4x4 Square() const {
    HMatrix4x4 square;
    square._data[0] = _data[0] * _data[0] + _data[1] * _data[1] +
                      _data[2] * _data[2] + _data[4] * _data[4] +
                      _data[5] * _data[5] + _data[9] * _data[9] +
                      _data[10] * _data[10];
    square.SetComplex(1, _data[0] * ToComplex(1) + ToComplex(1) * _data[3] +
                             ToComplex(4) * ToConjugate(6) +
                             ToComplex(9) * ToConjugate(11));
    square._data[3] = _data[1] * _data[1] + _data[2] * _data[2] +
                      _data[3] * _data[3] + _data[6] * _data[6] +
                      _data[7] * _data[7] + _data[11] * _data[11] +
                      _data[12] * _data[12];
    square.SetComplex(4, _data[0] * ToComplex(4) + ToComplex(1) * ToComplex(6) +
                             ToComplex(4) * _data[8] +
                             ToComplex(9) * ToConjugate(13));
    square.SetComplex(6, ToConjugate(1) * ToComplex(4) +
                             _data[3] * ToComplex(6) + ToComplex(6) * _data[8] +
                             ToComplex(11) * ToConjugate(13));
    square._data[8] = _data[4] * _data[4] + _data[5] * _data[5] +
                      _data[6] * _data[6] + _data[7] * _data[7] +
                      _data[8] * _data[8] + _data[13] * _data[13] +
                      _data[14] * _data[14];
    square.SetComplex(
        9, _data[0] * ToComplex(9) + ToComplex(1) * ToComplex(11) +
               ToComplex(4) * ToComplex(13) + ToComplex(9) * _data[15]);
    square.SetComplex(
        11, ToConjugate(1) * ToComplex(9) + _data[3] * ToComplex(11) +
                ToComplex(6) * ToComplex(13) + ToComplex(11) * _data[15]);
    square.SetComplex(
        13, ToConjugate(4) * ToComplex(9) + ToConjugate(6) * ToComplex(11) +
                _data[8] * ToComplex(13) + ToComplex(13) * _data[15]);
    square._data[15] = _data[9] * _data[9] + _data[10] * _data[10] +
                       _data[11] * _data[11] + _data[12] * _data[12] +
                       _data[13] * _data[13] + _data[14] * _data[14] +
                       _data[15] * _data[15];
    return square;
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    for (double d : _data) stream.Double(d);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    for (double& d : _data) stream.Double(d);
  }

 private:
  /**
   * Combines a real value and its next into a complex value.
   * @param single_index should be one of 1, 4, 6, 9, 11 or 13.
   * @see @ref Data() for info about the indices.
   */
  constexpr std::complex<double> ToComplex(size_t single_index) const {
    return std::complex<double>(_data[single_index], _data[single_index + 1]);
  }
  /**
   * Like @ref ToComplex(), but returns its conjugate.
   */
  constexpr std::complex<double> ToConjugate(size_t single_index) const {
    return std::complex<double>(_data[single_index], -_data[single_index + 1]);
  }
  /**
   * Set two consecutive real values to the specific complex value.
   * Counterpart of @ref ToComplex().
   */
  constexpr void SetComplex(size_t single_index,
                            std::complex<double> new_value) {
    _data[single_index] = new_value.real();
    _data[single_index + 1] = new_value.imag();
  }

  // See documentation for Data method on the internal data storage
  std::array<double, 16> _data;
};

typedef HMatrix4x4 HMC4x4;
}  // namespace aocommon

// Functions below require Matrix4x4 to be available. They're separated from the
// class because of a circular dependency.

#include "aocommon/matrix4x4.h"

namespace aocommon {

inline constexpr HMatrix4x4::HMatrix4x4(const aocommon::Matrix4x4& src)
    :  // row 0
      _data{src[0].real(),
            // row 1
            src[4].real(), src[4].imag(), src[5].real(),
            // row 2
            src[8].real(), src[8].imag(), src[9].real(), src[9].imag(),
            src[10].real(),
            // row 3
            src[12].real(), src[12].imag(), src[13].real(), src[13].imag(),
            src[14].real(), src[14].imag(), src[15].real()} {}

inline bool HMatrix4x4::Invert() {
  aocommon::Matrix4x4 inv = ToMatrix();
  if (!inv.Invert())
    return false;
  else {
    *this = HMatrix4x4(inv);
    return true;
  }
}

inline constexpr aocommon::Matrix4x4 HMatrix4x4::ToMatrix() const {
  aocommon::Matrix4x4 m;
  for (size_t i = 0; i != 16; ++i) {
    m[i] = (*this)[i];
  }
  return m;
}

}  // namespace aocommon

#endif
