#include <boost/test/unit_test.hpp>

#include <aocommon/imagecoordinates.h>

#include <cmath>

using aocommon::ImageCoordinates;

BOOST_AUTO_TEST_SUITE(image_coordinates)

BOOST_AUTO_TEST_CASE(mean_position) {
  std::pair<double, double> result =
      ImageCoordinates::MeanPosition(std::vector<std::pair<double, double>>{});
  BOOST_CHECK_EQUAL(result.first, 0.0);
  BOOST_CHECK_EQUAL(result.second, 0.0);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, 0.0}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{1.0, 0.0}});
  BOOST_CHECK_CLOSE_FRACTION(result.first, 1.0, 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{1.0, 1.0}});
  BOOST_CHECK_CLOSE_FRACTION(result.first, 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.second, 1.0, 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{-1.0, 0.0}});
  BOOST_CHECK_CLOSE_FRACTION(result.first, -1.0, 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{-1.0, -1.0}});
  BOOST_CHECK_CLOSE_FRACTION(result.first, -1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.second, -1.0, 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, M_PI * 0.5}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.second, M_PI * 0.5, 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, 0.0}, {0.0, 0.0}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{1.0, 1.0}, {1.0, 1.0}});
  BOOST_CHECK_CLOSE_FRACTION(result.first, 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.second, 1.0, 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, 0.0}, {0.0, -0.5 * M_PI}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.second, -0.25 * M_PI, 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, 0.0},
                                             {0.0, -0.25 * M_PI},
                                             {0.0, 0.25 * M_PI},
                                             {0.0, 0.0},
                                             {0.0, -0.25 * M_PI},
                                             {0.0, 0.25 * M_PI}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  result = ImageCoordinates::MeanPosition(
      std::vector<std::pair<double, double>>{{0.0, 0.0}});
  BOOST_CHECK_LT(std::abs(result.first), 1e-6);
  BOOST_CHECK_LT(std::abs(result.second), 1e-6);

  // Positions around the pole
  std::vector<std::pair<float, float>> around_pole;
  for (size_t i = 0; i != 20; ++i) {
    around_pole.emplace_back(2.0 * M_PI * i / 20.0 + 0.03, 0.5 * M_PI - 0.07);
  }
  // Here we also test the span overload
  result = ImageCoordinates::MeanPosition(std::span(around_pole));
  BOOST_CHECK(std::isfinite(result.first));
  BOOST_CHECK_CLOSE_FRACTION(result.second, 0.5 * M_PI, 1e-6);

  const std::vector<std::pair<double, double>> origin_test{{1.0, -0.5 * M_PI},
                                                           {0.0, 0.5 * M_PI}};
  BOOST_CHECK_NO_THROW(ImageCoordinates::MeanPosition(std::span(origin_test)));

  std::array<float, 2> array_result = ImageCoordinates::MeanPosition(
      std::vector<std::array<float, 2>>{{1.0, 1.0}, {1.0, 1.0}});
  BOOST_CHECK_CLOSE_FRACTION(array_result[0], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(array_result[1], 1.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
