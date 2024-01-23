
#include <aocommon/polarization.h>

#include <boost/test/unit_test.hpp>

using aocommon::Polarization;
using aocommon::PolarizationEnum;

BOOST_AUTO_TEST_SUITE(polarization)

BOOST_AUTO_TEST_CASE(parse_list) {
  BOOST_CHECK(Polarization::ParseList("").empty());

  std::set<PolarizationEnum> result = Polarization::ParseList("xx");
  BOOST_REQUIRE_EQUAL(result.size(), 1);
  BOOST_CHECK_EQUAL(*result.begin(), PolarizationEnum::XX);

  result = Polarization::ParseList("iquv");
  BOOST_CHECK_EQUAL(result.size(), 4);
  BOOST_CHECK(result.count(PolarizationEnum::StokesI) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesQ) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesU) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesV) == 1);

  result = Polarization::ParseList("xxxyyy");
  BOOST_CHECK_EQUAL(result.size(), 3);
  BOOST_CHECK(result.count(PolarizationEnum::XX) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::XY) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::YY) == 1);

  result = Polarization::ParseList("yy,rr,i,ll,v");
  BOOST_CHECK_EQUAL(result.size(), 5);
  BOOST_CHECK(result.count(PolarizationEnum::YY) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::RR) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesI) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::LL) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesV) == 1);

  result = Polarization::ParseList("I/RR");
  BOOST_CHECK_EQUAL(result.size(), 2);
  BOOST_CHECK(result.count(PolarizationEnum::StokesI) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::RR) == 1);

  result = Polarization::ParseList("Xy I Yx");
  BOOST_CHECK_EQUAL(result.size(), 3);
  BOOST_CHECK(result.count(PolarizationEnum::XY) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::StokesI) == 1);
  BOOST_CHECK(result.count(PolarizationEnum::YX) == 1);

  BOOST_CHECK_THROW(Polarization::ParseList("3"), std::runtime_error);
  BOOST_CHECK_THROW(Polarization::ParseList("iq3v"), std::runtime_error);
  BOOST_CHECK_THROW(Polarization::ParseList("  "), std::runtime_error);
  BOOST_CHECK_THROW(Polarization::ParseList("xx  yy"), std::runtime_error);
  BOOST_CHECK_THROW(Polarization::ParseList("x"), std::runtime_error);
  BOOST_CHECK_THROW(Polarization::ParseList("yyr"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
