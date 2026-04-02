#include "aocommon/optionalnumber.h"

#include <cstdint>

#include <boost/test/unit_test.hpp>

namespace {

using aocommon::OptionalNumber;

template <typename T>
class Test {
  static_assert(!OptionalNumber<T>());
  static_assert(!OptionalNumber<T>().HasValue());
  static_assert(!OptionalNumber<T>(std::nullopt));
  static_assert(!OptionalNumber<T>{});
  static_assert(OptionalNumber<T>(T(0)));
  static_assert(OptionalNumber<T>(0).HasValue());
  static_assert(OptionalNumber<T>(T(13)).Value() == T(13));

  static constexpr OptionalNumber<T> GetAssigned() {
    OptionalNumber<T> number;
    return number = T(0.0);
  }

  static_assert(GetAssigned() == 0);
  static_assert(*GetAssigned() == 0);

  static constexpr OptionalNumber<T> GetReset() {
    OptionalNumber<T> number = GetAssigned();
    number.Reset();
    return number;
  }

  static_assert(!GetReset());

  static constexpr OptionalNumber<T> GetUnassigned() {
    OptionalNumber<T> number(T(1));
    number = {};
    return number;
  }

  static constexpr OptionalNumber<T> GetSwappedLhs(OptionalNumber<T> a,
                                                   OptionalNumber<T> b) {
    swap(a, b);
    return a;
  }

  static constexpr OptionalNumber<T> GetSwappedRhs(OptionalNumber<T> a,
                                                   OptionalNumber<T> b) {
    swap(a, b);
    return b;
  }

  static_assert(!GetUnassigned());

  static_assert(*OptionalNumber(OptionalNumber<T>(5)) == 5);

  static_assert(OptionalNumber<T>().ValueOr(3) == 3);
  static_assert(OptionalNumber<T>(2).ValueOr(3) == 2);

  // Comparisons with optional numbers on both sides
  static_assert(OptionalNumber<T>() == OptionalNumber<T>());
  static_assert(OptionalNumber<T>(3) == OptionalNumber<T>(3));
  static_assert(!(OptionalNumber<T>() == OptionalNumber<T>(0)));
  static_assert(!(OptionalNumber<T>(0) == OptionalNumber<T>()));

  static_assert(OptionalNumber<T>() != OptionalNumber<T>(3));
  static_assert(OptionalNumber<T>(3) != OptionalNumber<T>());
  static_assert(OptionalNumber<T>(0) != OptionalNumber<T>(3));
  static_assert(!(OptionalNumber<T>() != OptionalNumber<T>()));
  static_assert(!(OptionalNumber<T>(3) != OptionalNumber<T>(3)));

  static_assert(OptionalNumber<T>() <= OptionalNumber<T>());
  static_assert(OptionalNumber<T>() <= OptionalNumber<T>(3));
  static_assert(OptionalNumber<T>(3) <= OptionalNumber<T>(3));
  static_assert(!(OptionalNumber<T>(3) <= OptionalNumber<T>()));
  static_assert(!(OptionalNumber<T>(3) <= OptionalNumber<T>(0)));

  static_assert(OptionalNumber<T>() < OptionalNumber<T>(3));
  static_assert(OptionalNumber<T>(0) < OptionalNumber<T>(3));
  static_assert(!(OptionalNumber<T>(0) < OptionalNumber<T>()));
  static_assert(!(OptionalNumber<T>(3) < OptionalNumber<T>(0)));

  static_assert(OptionalNumber<T>(3) >= OptionalNumber<T>());
  static_assert(OptionalNumber<T>() >= OptionalNumber<T>());
  static_assert(!(OptionalNumber<T>() >= OptionalNumber<T>(3)));
  static_assert(!(OptionalNumber<T>(0) >= OptionalNumber<T>(3)));

  static_assert(OptionalNumber<T>(3) > OptionalNumber<T>());
  static_assert(OptionalNumber<T>(3) > OptionalNumber<T>(0));
  static_assert(!(OptionalNumber<T>() > OptionalNumber<T>(0)));
  static_assert(!(OptionalNumber<T>(0) > OptionalNumber<T>(0)));

  // Comparisons with optional number lhs and literal rhs
  static_assert(OptionalNumber<T>(3) == 3);
  static_assert(!(OptionalNumber<T>() == 3));
  static_assert(!(OptionalNumber<T>(3) == 0));

  static_assert(OptionalNumber<T>() != 3);
  static_assert(OptionalNumber<T>(0) != 3);
  static_assert(!(OptionalNumber<T>(3) != 3));

  static_assert(OptionalNumber<T>() <= 3);
  static_assert(OptionalNumber<T>(0) <= 3);
  static_assert(!(OptionalNumber<T>(3) <= 0));

  static_assert(OptionalNumber<T>() < 0);
  static_assert(OptionalNumber<T>(0) < 3);
  static_assert(!(OptionalNumber<T>(3) < 0));

  static_assert(OptionalNumber<T>(0) >= 0);
  static_assert(OptionalNumber<T>(3) >= 3);
  static_assert(!(OptionalNumber<T>() >= 0));
  static_assert(!(OptionalNumber<T>(0) >= 3));

  static_assert(OptionalNumber<T>(3) > 0);
  static_assert(!(OptionalNumber<T>(0) > 3));
  static_assert(!(OptionalNumber<T>() > 0));

  static_assert(GetSwappedLhs(OptionalNumber<T>(5), OptionalNumber<T>(5)) == 5);
  static_assert(GetSwappedLhs(OptionalNumber<T>(5), OptionalNumber<T>(1)) == 1);
  static_assert(GetSwappedRhs(OptionalNumber<T>(5), OptionalNumber<T>(1)) == 5);
  static_assert(GetSwappedLhs(OptionalNumber<T>(), OptionalNumber<T>(5)) == 5);
  static_assert(!GetSwappedRhs(OptionalNumber<T>(), OptionalNumber<T>(5)));
  static_assert(!GetSwappedLhs(OptionalNumber<T>(5), OptionalNumber<T>()));
  static_assert(GetSwappedRhs(OptionalNumber<T>(5), OptionalNumber<T>()) == 5);
  static_assert(!GetSwappedLhs(OptionalNumber<T>(), OptionalNumber<T>()));
  static_assert(!GetSwappedRhs(OptionalNumber<T>(), OptionalNumber<T>()));
};

template class Test<unsigned char>;
template class Test<int>;
template class Test<std::size_t>;
template class Test<float>;
template class Test<double>;

static_assert(OptionalNumber<std::uint8_t>::UnsetValue == 255);
static_assert(OptionalNumber<std::int8_t>::UnsetValue == -128);
static_assert(OptionalNumber<std::uint16_t>::UnsetValue == 65535);
static_assert(OptionalNumber<std::int16_t>::UnsetValue == -32768);
static_assert(OptionalNumber<float>::UnsetValue ==
              std::numeric_limits<float>::lowest());
static_assert(OptionalNumber<double>::UnsetValue ==
              std::numeric_limits<double>::lowest());

BOOST_AUTO_TEST_SUITE(optional_number)

BOOST_AUTO_TEST_CASE(from_std_optional) {
  const std::optional<size_t> std_unset;
  const std::optional<size_t> std_five(5);
  aocommon::OptionalNumber<size_t> value_a(std_five);
  BOOST_CHECK(value_a.HasValue());
  BOOST_CHECK_EQUAL(*std_five, *value_a);

  aocommon::OptionalNumber<size_t> value_b(std_unset);
  BOOST_CHECK(!value_b.HasValue());

  const std::optional<size_t> std_six(6);
  value_a = std_six;
  BOOST_CHECK(value_a.HasValue());
  BOOST_CHECK_EQUAL(*std_six, *value_a);

  // from unset to set
  value_b = std_six;
  BOOST_CHECK(value_b.HasValue());
  BOOST_CHECK_EQUAL(*std_six, *value_b);

  // from set to unset
  value_a = std_unset;
  BOOST_CHECK(!value_a.HasValue());
  // from unset to unset
  value_a = std_unset;
  BOOST_CHECK(!value_a.HasValue());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
