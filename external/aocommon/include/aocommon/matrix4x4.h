#ifndef AOCOMMON_MATRIX_4X4_H_
#define AOCOMMON_MATRIX_4X4_H_

#include <cassert>
#include <cmath>
#include <complex>
#include <string>
#include <sstream>
#include <stdexcept>

#include "matrix2x2.h"

#include "scalar/vector4.h"

namespace aocommon {

class HMatrix4x4;

class Matrix4x4 {
 public:
  constexpr Matrix4x4() {}

  constexpr Matrix4x4(std::initializer_list<std::complex<double>> list) {
    assert(list.size() == 16);
    size_t index = 0;
    for (const std::complex<double>& el : list) {
      _data[index] = el;
      ++index;
    }
  }

  static constexpr Matrix4x4 Zero() { return Matrix4x4(); }

  static constexpr Matrix4x4 Unit() {
    return Matrix4x4{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                     0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  }

  Matrix4x4 operator+(const Matrix4x4& rhs) const {
    Matrix4x4 result;
    for (size_t i = 0; i != 16; ++i) result[i] = this->_data[i] + rhs._data[i];
    return result;
  }

  Matrix4x4& operator+=(const Matrix4x4& rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] += rhs._data[i];
    return *this;
  }

  Matrix4x4 operator*(const std::complex<double>& rhs) const {
    Matrix4x4 m;
    for (size_t i = 0; i != 16; ++i) m[i] = _data[i] * rhs;
    return m;
  }

  Vector4 operator*(const Vector4& rhs) const {
    Vector4 v(_data[0] * rhs[0], _data[4] * rhs[0], _data[8] * rhs[0],
              _data[12] * rhs[0]);
    for (size_t i = 1; i != 4; ++i) {
      v[0] += _data[i] * rhs[i];
      v[1] += _data[i + 4] * rhs[i];
      v[2] += _data[i + 8] * rhs[i];
      v[3] += _data[i + 12] * rhs[i];
    }
    return v;
  }

  Matrix4x4 operator*(const Matrix4x4& rhs) const {
    return {_data[0] * rhs[0] + _data[1] * rhs[4] + _data[2] * rhs[8] +
                _data[3] * rhs[12],
            _data[0] * rhs[1] + _data[1] * rhs[5] + _data[2] * rhs[9] +
                _data[3] * rhs[13],
            _data[0] * rhs[2] + _data[1] * rhs[6] + _data[2] * rhs[10] +
                _data[3] * rhs[14],
            _data[0] * rhs[3] + _data[1] * rhs[7] + _data[2] * rhs[11] +
                _data[3] * rhs[15],
            _data[4] * rhs[0] + _data[5] * rhs[4] + _data[6] * rhs[8] +
                _data[7] * rhs[12],
            _data[4] * rhs[1] + _data[5] * rhs[5] + _data[6] * rhs[9] +
                _data[7] * rhs[13],
            _data[4] * rhs[2] + _data[5] * rhs[6] + _data[6] * rhs[10] +
                _data[7] * rhs[14],
            _data[4] * rhs[3] + _data[5] * rhs[7] + _data[6] * rhs[11] +
                _data[7] * rhs[15],
            _data[8] * rhs[0] + _data[9] * rhs[4] + _data[10] * rhs[8] +
                _data[11] * rhs[12],
            _data[8] * rhs[1] + _data[9] * rhs[5] + _data[10] * rhs[9] +
                _data[11] * rhs[13],
            _data[8] * rhs[2] + _data[9] * rhs[6] + _data[10] * rhs[10] +
                _data[11] * rhs[14],
            _data[8] * rhs[3] + _data[9] * rhs[7] + _data[10] * rhs[11] +
                _data[11] * rhs[15],
            _data[12] * rhs[0] + _data[13] * rhs[4] + _data[14] * rhs[8] +
                _data[15] * rhs[12],
            _data[12] * rhs[1] + _data[13] * rhs[5] + _data[14] * rhs[9] +
                _data[15] * rhs[13],
            _data[12] * rhs[2] + _data[13] * rhs[6] + _data[14] * rhs[10] +
                _data[15] * rhs[14],
            _data[12] * rhs[3] + _data[13] * rhs[7] + _data[14] * rhs[11] +
                _data[15] * rhs[15]};
  }

  Matrix4x4 HermTranspose() const {
    auto C = [](std::complex<double> z) -> std::complex<double> {
      return std::conj(z);
    };
    return Matrix4x4{C(_data[0]), C(_data[4]), C(_data[8]),  C(_data[12]),
                     C(_data[1]), C(_data[5]), C(_data[9]),  C(_data[13]),
                     C(_data[2]), C(_data[6]), C(_data[10]), C(_data[14]),
                     C(_data[3]), C(_data[7]), C(_data[11]), C(_data[15])};
  }

  bool Invert() {
    std::complex<double> inv[16];
    const std::complex<double>* m = _data;

    inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] +
             m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

    inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] +
             m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] +
             m[12] * m[7] * m[10];

    inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] +
             m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

    inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] +
              m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] +
              m[12] * m[6] * m[9];

    inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] +
             m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] +
             m[13] * m[3] * m[10];

    inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] +
             m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

    inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] -
             m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

    inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] +
              m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

    inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] +
             m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

    inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] -
             m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

    inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] +
              m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

    inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] -
              m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] -
             m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] -
              m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

    std::complex<double> det =
        m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0.0) return false;

    det = 1.0 / det;

    for (size_t i = 0; i < 16; i++) _data[i] = inv[i] * det;

    return true;
  }

  constexpr std::complex<double>& operator[](size_t i) { return _data[i]; }

  constexpr const std::complex<double>& operator[](size_t i) const {
    return _data[i];
  }

  double Norm() const {
    double n = 0.0;
    for (size_t i = 0; i != 16; ++i) {
      n += std::norm(_data[i]);
    }
    return n;
  }

  /// @returns A^2.
  Matrix4x4 Square() const { return (*this) * (*this); }

  /**
   * Calculates A^H A, i.e. the Hermitian transpose of itself multiplied by
   * itself. The result is Hermitian (also for non-Hermitian inputs) and
   * therefore a HMatrix4x4 is returned. This avoids calculating one half of the
   * result and is therefore faster than A.HermTranspose() * A.
   */
  HMatrix4x4 HermitianSquare() const;

  std::string String() const {
    std::ostringstream str;
    for (size_t y = 0; y != 4; ++y) {
      for (size_t x = 0; x != 3; ++x) {
        str << _data[x + y * 4] << '\t';
      }
      str << _data[3 + y * 4] << '\n';
    }
    return str.str();
  }

  static Matrix4x4 KroneckerProduct(const MC2x2& veca, const MC2x2& vecb) {
    Matrix4x4 result;
    const size_t posa[4] = {0, 2, 8, 10};
    for (size_t i = 0; i != 4; ++i) {
      result[posa[i]] = veca.Get(i) * vecb.Get(0);
      result[posa[i] + 1] = veca.Get(i) * vecb.Get(1);
      result[posa[i] + 4] = veca.Get(i) * vecb.Get(2);
      result[posa[i] + 5] = veca.Get(i) * vecb.Get(3);
    }
    return result;
  }

 private:
  std::complex<double> _data[16];
};

typedef Matrix4x4 MC4x4;

}  // namespace aocommon

// Functions below require HMatrix4x4 to be available. They're separated from
// the class because of a circular dependency.

#include "hmatrix4x4.h"

namespace aocommon {

inline HMatrix4x4 Matrix4x4::HermitianSquare() const {
  auto C = [](std::complex<double> z) -> std::complex<double> {
    return std::conj(z);
  };
  auto N = [](std::complex<double> z) -> double { return std::norm(z); };
  return {N(_data[0]) + N(_data[4]) + N(_data[8]) + N(_data[12]),
          0.0,
          0.0,
          0.0,
          C(_data[1]) * _data[0] + C(_data[5]) * _data[4] +
              C(_data[9]) * _data[8] + C(_data[13]) * _data[12],
          N(_data[1]) + N(_data[5]) + N(_data[9]) + N(_data[13]),
          0.0,
          0.0,
          C(_data[2]) * _data[0] + C(_data[6]) * _data[4] +
              C(_data[10]) * _data[8] + C(_data[14]) * _data[12],
          C(_data[2]) * _data[1] + C(_data[6]) * _data[5] +
              C(_data[10]) * _data[9] + C(_data[14]) * _data[13],
          N(_data[2]) + N(_data[6]) + N(_data[10]) + N(_data[14]),
          0.0,
          C(_data[3]) * _data[0] + C(_data[7]) * _data[4] +
              C(_data[11]) * _data[8] + C(_data[15]) * _data[12],
          C(_data[3]) * _data[1] + C(_data[7]) * _data[5] +
              C(_data[11]) * _data[9] + C(_data[15]) * _data[13],
          C(_data[3]) * _data[2] + C(_data[7]) * _data[6] +
              C(_data[11]) * _data[10] + C(_data[15]) * _data[14],
          N(_data[3]) + N(_data[7]) + N(_data[11]) + N(_data[15])};
}

}  // namespace aocommon

#endif
