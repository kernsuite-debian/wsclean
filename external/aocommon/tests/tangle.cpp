#include <boost/test/unit_test.hpp>

#include <aocommon/units/angle.h>

#include <cmath>

using aocommon::units::Angle;

BOOST_AUTO_TEST_SUITE(angle)

BOOST_AUTO_TEST_CASE(parse_values) {
  BOOST_CHECK_EQUAL(Angle::Parse("0", ""), 0.0);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("90", "", Angle::kDegrees), M_PI_2,
                             1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("2700ArCmIn", ""), M_PI_4, 1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("5400", "", Angle::kArcminutes),
                             M_PI_2, 1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("1296000ASEC", ""), 2.0 * M_PI, 1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("1296000mas", ""),
                             2.0 * M_PI / 1000.0, 1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("90 deg", ""), M_PI_2, 1e-7);
  BOOST_CHECK_CLOSE_FRACTION(Angle::Parse("-180 deg", ""), -M_PI, 1e-7);
  BOOST_CHECK_THROW(Angle::Parse("10 aocommon", ""), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
