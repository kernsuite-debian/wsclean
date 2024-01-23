#include "aocommon/transform_if.h"

#include <boost/test/unit_test.hpp>

// Since boost test can't do compile-time validation the tests use asserts.
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <array>
#include <cassert>

#if __cplusplus > 201703L
// std::equal is constexpr since C++20
using std::equal;
#else
template <class InputIt1, class InputIt2>
#if __cplusplus > 201402L
constexpr
#endif
    bool
    equal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
  while (first1 != last1) {
    if (!(*first1 == *first2)) {
      return false;
    }
    ++first1;
    ++first2;
  }
  return true;
}
#endif

static constexpr std::array<int, 10> kInput{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

BOOST_AUTO_TEST_SUITE(transform_if)

#if __cplusplus > 201402L
constexpr
#endif
    bool
    TestAllFalse() {
  std::array<int, 10> output{};
  std::array<int, 10>::iterator out = aocommon::transform_if(
      kInput.begin(), kInput.end(), output.begin(), [](int) { return false; },
      [](int i) { return i; });

  assert(out == output.begin());

  return true;
}

BOOST_AUTO_TEST_CASE(transform_if_all_false) {
  TestAllFalse();
#if __cplusplus > 201402L
  static_assert(TestAllFalse());
#endif
}

#if __cplusplus > 201402L
constexpr
#endif
    bool
    TestAllTrue() {
  std::array<int, 10> output{};
  std::array<int, 10>::iterator out = aocommon::transform_if(
      kInput.begin(), kInput.end(), output.begin(), [](int) { return true; },
      [](int i) { return i; });

  assert(out == output.end());
  assert(equal(kInput.begin(), kInput.end(), output.begin()));

  return true;
}

BOOST_AUTO_TEST_CASE(transform_if_all_true) {
  TestAllTrue();
#if __cplusplus > 201402L
  static_assert(TestAllTrue());
#endif
}

#if __cplusplus > 201402L
constexpr
#endif
    bool
    TestLessThanFive() {
  std::array<int, 10> output{};
  std::array<int, 10>::iterator out = aocommon::transform_if(
      kInput.begin(), kInput.end(), output.begin(), [](int i) { return i < 5; },
      [](int i) { return i; });

  assert(out == output.begin() + 5);
  assert(equal(kInput.begin(), kInput.begin() + 5, output.begin()));

  return true;
}

BOOST_AUTO_TEST_CASE(transform_if_less_than_five) {
  TestLessThanFive();
#if __cplusplus > 201402L
  static_assert(TestLessThanFive());
#endif
}

#if __cplusplus > 201402L
constexpr
#endif
    bool
    TestLessThanThreeToDouble() {
  std::array<double, 10> output{};
  std::array<double, 10>::iterator out = aocommon::transform_if(
      kInput.begin(), kInput.end(), output.begin(), [](int i) { return i < 3; },
      [](int i) { return i + 3; });

  assert(out == output.begin() + 3);
  assert(equal(kInput.begin() + 3, kInput.begin() + 6, output.begin()));

  return true;
}

BOOST_AUTO_TEST_CASE(transform_if_less_than_three_to_double) {
  TestLessThanThreeToDouble();
#if __cplusplus > 201402L
  static_assert(TestLessThanThreeToDouble());
#endif
}

#if __cplusplus > 201402L
constexpr
#endif
    bool
    TestToUnrelatedType() {
  constexpr std::array<const char*, 3> kData{"abc", "def", "ghi"};
  std::array<const char*, 3> output{};
  std::array<const char*, 3>::iterator out = aocommon::transform_if(
      kInput.begin(), kInput.end(), output.begin(),
      [](int i) { return i == 5; }, [&](int i) { return kData[1]; });

  assert(out == output.begin() + 1);
  assert(output.front() == kData[1]);

  return true;
}

BOOST_AUTO_TEST_CASE(transform_if_to_unrelated_type) {
  TestToUnrelatedType();
#if __cplusplus > 201402L
  static_assert(TestToUnrelatedType());
#endif
}

BOOST_AUTO_TEST_SUITE_END()
