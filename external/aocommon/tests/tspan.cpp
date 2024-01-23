#include <aocommon/xt/span.h>
#include <aocommon/xt/utensor.h>

#include <algorithm>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include <xtensor/xbuilder.hpp>
#include <xtensor/xtensor.hpp>

using aocommon::xt::CreateSpan;
using aocommon::xt::Span;

namespace {
const std::array<size_t, 2> kShape{2, 3};
const size_t kSize = kShape[0] * kShape[1];
}  // namespace

BOOST_AUTO_TEST_SUITE(span)

BOOST_AUTO_TEST_CASE(non_const_span) {
  float values[kSize] = {42.0, 42.0, 42.0, 42.0, 42.0, 42.0};
  Span<float, 2> span = CreateSpan(&values[0], kShape);
  static_assert(std::is_reference_v<decltype(span(0, 0))>);
  static_assert(
      !std::is_const_v<std::remove_reference_t<decltype(span(0, 0))>>);
}

BOOST_AUTO_TEST_CASE(const_span) {
  const float values[kSize] = {42.0, 42.0, 42.0, 42.0, 42.0, 42.0};
  Span<const float, 2> span = CreateSpan(&values[0], kShape);
  static_assert(std::is_reference_v<decltype(span(0, 0))>);
  static_assert(std::is_const_v<std::remove_reference_t<decltype(span(0, 0))>>);
}

BOOST_AUTO_TEST_CASE(raw_pointer) {
  {
    double values[kSize] = {42.0, 42.0, 42.0, 42.0, 42.0, 42.0};
    Span<double, 2> span = CreateSpan(&values[0], kShape);
    BOOST_TEST(span == xt::full_like(xt::ones<double>(kShape), 42.0));
    BOOST_TEST(span.data() == &values[0]);
  }

  {
    const char values[kSize] = {1, 2, 3, 4, 5, 6};
    Span<const char, 2> span = CreateSpan(&values[0], kShape);
    xt::xtensor<char, 2> tensor = {{1, 2, 3}, {4, 5, 6}};
    BOOST_TEST(span == tensor);
    BOOST_TEST(span.data() == &values[0]);
  }
}

BOOST_AUTO_TEST_CASE(vector) {
  {
    std::vector<int> vector{42, 42, 42, 42, 42, 42};
    Span<int, 2> span = CreateSpan(vector, kShape);
    BOOST_TEST(span == xt::full_like(xt::ones<int>(kShape), 42));
    BOOST_TEST(span.data() == vector.data());
  }
  {
    const std::vector<std::size_t> vector{1, 2, 3, 4, 5, 6};
    Span<const std::size_t, 2> span = CreateSpan(vector, kShape);
    xt::xtensor<std::size_t, 2> tensor = {{1, 2, 3}, {4, 5, 6}};
    BOOST_TEST(span == tensor);
    BOOST_TEST(span.data() == vector.data());
  }
}

namespace {
template <typename TensorType>
void TestTensor1D() {
  TensorType tensor1d{41.0, 42.0, 43.0, 44.0, 45.0};
  Span<typename TensorType::value_type, 1> span = CreateSpan(tensor1d);
  BOOST_TEST(span == tensor1d);
  BOOST_TEST(span.data() == tensor1d.data());
}

template <typename TensorType>
void TestTensor2D() {
  TensorType tensor2d_square{
      {true, false, true}, {false, true, true}, {true, false, false}};
  Span<bool, 2> span = CreateSpan(tensor2d_square);
  BOOST_TEST(span == tensor2d_square);
  BOOST_TEST(span.data() == tensor2d_square.data());
}

template <typename TensorType>
void TestTensor3D() {
  const TensorType tensor3d_const{{{1, 2}, {3, 4}, {5, 6}},
                                  {{7, 8}, {9, 10}, {11, 12}},
                                  {{13, 14}, {15, 16}, {17, 18}},
                                  {{19, 20}, {21, 22}, {23, 24}}};
  Span<const typename TensorType::value_type, 3> span =
      CreateSpan(tensor3d_const);
  BOOST_TEST(span == tensor3d_const);
  BOOST_TEST(span.data() == tensor3d_const.data());
}
}  // namespace

BOOST_AUTO_TEST_CASE(tensor) {
  TestTensor1D<::xt::xtensor<float, 1>>();
  TestTensor2D<::xt::xtensor<bool, 2>>();
  TestTensor3D<::xt::xtensor<int, 3>>();
}

BOOST_AUTO_TEST_CASE(utensor) {
  TestTensor1D<aocommon::xt::UTensor<long double, 1>>();
  TestTensor2D<aocommon::xt::UTensor<bool, 2>>();
  TestTensor3D<aocommon::xt::UTensor<short, 3>>();
}

BOOST_AUTO_TEST_SUITE_END()
