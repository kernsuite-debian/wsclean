#include <aocommon/py/uniqueptr.h>

#include <complex>
#include <type_traits>

#include <boost/test/unit_test.hpp>

using aocommon::py::PyUniquePointer;

BOOST_AUTO_TEST_SUITE(pyuniquepointer)

BOOST_AUTO_TEST_CASE(construct_from_unique_pointer) {
  static_assert(!std::is_copy_assignable_v<PyUniquePointer<int>>);
  static_assert(!std::is_copy_constructible_v<PyUniquePointer<int>>);
  static_assert(std::is_move_assignable_v<PyUniquePointer<int>>);
  static_assert(std::is_move_constructible_v<PyUniquePointer<int>>);

  auto i = std::make_unique<int>(42);
  auto j = std::make_unique<int>(43);
  int* i_raw = i.get();
  int* j_raw = j.get();
  PyUniquePointer<int> pypointer(std::move(i));
  const PyUniquePointer<int> pypointer_const(std::move(j));

  // Test operator*
  BOOST_TEST(&(*pypointer) == i_raw);
  BOOST_TEST(&(*pypointer_const) == j_raw);
  BOOST_TEST(*pypointer == 42);
  BOOST_TEST(*pypointer_const == 43);

  // Test operator->
  BOOST_TEST(pypointer.operator->() == i_raw);
  BOOST_TEST(pypointer_const.operator->() == j_raw);
}

BOOST_AUTO_TEST_CASE(construct_from_value) {
  const PyUniquePointer<std::complex<float>> pypointer_default;
  BOOST_TEST(*pypointer_default == std::complex<float>(0.0f, 0.0f));

  const PyUniquePointer<std::complex<float>> pypointer_one_value{42.0f};
  BOOST_TEST(*pypointer_one_value == std::complex<float>(42.0f, 0.0f));

  const PyUniquePointer<std::complex<float>> pypointer_two_values{42.0f,
                                                                  -43.0f};
  BOOST_TEST(*pypointer_two_values == std::complex<float>(42.0f, -43.0f));
}

BOOST_AUTO_TEST_CASE(take) {
  auto i = std::make_unique<int>(42);
  int* i_raw = i.get();
  PyUniquePointer<int> pypointer(std::move(i));

  i = pypointer.take();
  BOOST_TEST(i.get() == i_raw);
  BOOST_CHECK_THROW(pypointer.take(), std::runtime_error);
  BOOST_CHECK_THROW((void)*pypointer, std::runtime_error);
  BOOST_CHECK_THROW((void)pypointer.operator->(), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
