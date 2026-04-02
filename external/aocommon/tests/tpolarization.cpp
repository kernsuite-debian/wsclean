
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

BOOST_AUTO_TEST_CASE(convert_linear) {
  constexpr std::complex<double> linear_values[] = {
      {10.0, 1.0}, {-11.0, 2.0}, {12.0, -3.0}, {-13.0, -4.0}};
  BOOST_CHECK_EQUAL(
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::XX),
      linear_values[0]);
  BOOST_CHECK_EQUAL(
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::XY),
      linear_values[1]);
  BOOST_CHECK_EQUAL(
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::YX),
      linear_values[2]);
  BOOST_CHECK_EQUAL(
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::YY),
      linear_values[3]);
  double real_stokes_values[4];
  Polarization::LinearToStokes(linear_values, real_stokes_values);
  const std::complex<double> i =
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::StokesI);
  const std::complex<double> q =
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::StokesQ);
  const std::complex<double> u =
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::StokesU);
  const std::complex<double> v =
      Polarization::ConvertFromLinear(linear_values, PolarizationEnum::StokesV);
  BOOST_CHECK_EQUAL(i.real(), real_stokes_values[0]);
  BOOST_CHECK_EQUAL(q.real(), real_stokes_values[1]);
  BOOST_CHECK_EQUAL(u.real(), real_stokes_values[2]);
  BOOST_CHECK_EQUAL(v.real(), real_stokes_values[3]);
  std::complex<double> linear_i[4];
  Polarization::ConvertToLinear(i, PolarizationEnum::StokesI, linear_i);
  std::complex<double> linear_q[4];
  Polarization::ConvertToLinear(q, PolarizationEnum::StokesQ, linear_q);
  std::complex<double> linear_u[4];
  Polarization::ConvertToLinear(u, PolarizationEnum::StokesU, linear_u);
  std::complex<double> linear_v[4];
  Polarization::ConvertToLinear(v, PolarizationEnum::StokesV, linear_v);
  std::complex<double> linear_reconverted_values[4];
  for (size_t p = 0; p != 4; ++p) {
    linear_reconverted_values[p] =
        linear_i[p] + linear_q[p] + linear_u[p] + linear_v[p];
    BOOST_CHECK_CLOSE_FRACTION(linear_values[p].real(),
                               linear_reconverted_values[p].real(), 1e-8);
    BOOST_CHECK_CLOSE_FRACTION(linear_values[p].imag(),
                               linear_reconverted_values[p].imag(), 1e-8);
  }
  Polarization::StokesToLinear(real_stokes_values, linear_reconverted_values);
  Polarization::LinearToStokes(linear_reconverted_values, real_stokes_values);
  BOOST_CHECK_CLOSE_FRACTION(i.real(), real_stokes_values[0], 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(q.real(), real_stokes_values[1], 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(u.real(), real_stokes_values[2], 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(v.real(), real_stokes_values[3], 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
