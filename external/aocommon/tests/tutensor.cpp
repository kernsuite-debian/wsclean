#include <algorithm>
#include <cstddef>
#include <memory>
#include <complex>

#include <aocommon/xt/utensor.h>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(utensor)

BOOST_AUTO_TEST_CASE(create_complex) {
  // This allocator wraps an std::allocator. After allocating memory, it fills
  // the allocated memory with fixed non-zero values.
  class InitializingAllocator {
   public:
    using value_type = std::complex<float>;

    value_type* allocate(std::size_t n) {
      value_type* result = allocator_.allocate(n);
      std::fill_n(result, n, value_type{42.0f, 43.0f});
      return result;
    }

    void deallocate(value_type* p, std::size_t n) {
      allocator_.deallocate(p, n);
    }

    bool operator==(const InitializingAllocator& other) const {
      return allocator_ == other.allocator_;
    }

    bool operator!=(const InitializingAllocator& other) const {
      return allocator_ != other.allocator_;
    }

   private:
    std::allocator<value_type> allocator_;
  };

  const std::size_t kTensorSize = 3;
  const std::array<std::size_t, 1> kShape{kTensorSize};
  xt::xtensor<std::complex<float>, 1, XTENSOR_DEFAULT_LAYOUT,
              InitializingAllocator>
      initialized(kShape);
  aocommon::xt::UTensor<std::complex<float>, 1, XTENSOR_DEFAULT_LAYOUT,
                        InitializingAllocator>
      uninitialized(kShape);

  for (std::size_t i = 0; i < kTensorSize; ++i) {
    // The allocator filled the allocated space with {42.0f, 43.0f}.

    // The regular 'initialized' xtensor overwrites the contents.
    BOOST_CHECK_EQUAL(initialized(i).real(), 0.0f);
    BOOST_CHECK_EQUAL(initialized(i).imag(), 0.0f);

    // The 'uninitialized' UTensor should not have overwritten the contents.
    BOOST_CHECK_EQUAL(uninitialized(i).real(), 42.0f);
    BOOST_CHECK_EQUAL(uninitialized(i).imag(), 43.0f);
  }
}

BOOST_AUTO_TEST_SUITE_END()
