#ifndef AOCOMMON_SCALAR_VECTOR4_H_
#define AOCOMMON_SCALAR_VECTOR4_H_

namespace aocommon {

/**
 * Class implements a vector with length 4 of complex-valued doubles.
 *
 */
class Vector4 {
 public:
  /**
   * Default constructor
   *
   */
  Vector4(){};

  /**
   * Construct a Vector4 object from four complex-valued numbers.
   *
   */
  Vector4(std::complex<double> a, std::complex<double> b,
          std::complex<double> c, std::complex<double> d) {
    _data[0] = a;
    _data[1] = b;
    _data[2] = c;
    _data[3] = d;
  };

  /**
   * Index Vector4 object
   */
  std::complex<double>& operator[](size_t i) { return _data[i]; }
  const std::complex<double>& operator[](size_t i) const { return _data[i]; }

  /**
   * Get pointer to underlying data
   */
  std::complex<double>* data() { return _data; }
  const std::complex<double>* data() const { return _data; }

 private:
  std::complex<double> _data[4];
};

}  // namespace aocommon

#endif
