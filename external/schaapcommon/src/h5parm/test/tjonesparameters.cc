// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "jonesparameters.h"
#include "h5parm.h"

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::GainType;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::JonesParameters;
using schaapcommon::h5parm::SolTab;

const std::vector<double> kFreqs{130e6, 131e6};
const std::vector<double> kTimes{0., 1.};
const unsigned int kNAnts = 4;
const JonesParameters::InterpolationType kInterpolationType =
    JonesParameters::InterpolationType::LINEAR;
const hsize_t kDirection = 42;

class SolTabMock : public SolTab {
 public:
  SolTabMock() : called(0) {}
  std::vector<double> GetValuesOrWeights(const std::string& val_or_weight,
                                         const std::string& ant_name,
                                         const std::vector<double>& times,
                                         const std::vector<double>& freqs,
                                         unsigned int pol, unsigned int dir,
                                         bool nearest) const override {
    ++called;
    auto res = std::vector<double>(kNAnts, 200.);

    if (ant_name.back() - '0' >= int(kNAnts)) {
      // Number represented by last character of ant_name is >= kNants
      // E.g. an antenna 'Antenna5' is requested which is not in the soltab
      throw(std::runtime_error("SolTab has no element Antenna" +
                               std::string(1, ant_name.back()) +
                               " in antenna"));
    }

    if (val_or_weight == "val") {
      res.back() = 100.;
    } else {
      res.back() = 0.;
    }
    return res;
  }

  mutable int called;
};

JonesParameters PrepareJonesParameters(GainType ct, bool invert = false,
                                       float sigma_mmse = 0.,
                                       unsigned int parm_size = 0) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();
  SolTabMock mock2 = SolTabMock();

  JonesParameters jones(kFreqs, kTimes, antNames, ct, kInterpolationType,
                        kDirection, &mock, &mock2, invert, sigma_mmse,
                        parm_size);
  return jones;
}

BOOST_AUTO_TEST_SUITE(jonesparameters)

BOOST_AUTO_TEST_CASE(make_complex_gain) {
  JonesParameters jones = PrepareJonesParameters(GainType::kDiagonalComplex);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_EQUAL(parms.shape()[1], kNAnts);
  BOOST_CHECK_EQUAL(parms.shape()[2], kTimes.size() * kFreqs.size());
  // Amplitude and phase are 200
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 97.437, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -174.659, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_scalar_gain) {
  JonesParameters jones = PrepareJonesParameters(GainType::kScalarComplex);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_EQUAL(parms.shape()[1], kNAnts);
  BOOST_CHECK_EQUAL(parms.shape()[2], kTimes.size() * kFreqs.size());
  // Amplitude and phase are 200
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 97.437, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -174.659, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_fulljones) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kFullJones, true, 1.);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 4);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 0.5006091, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0.00108991, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_tec) {
  JonesParameters jones = PrepareJonesParameters(GainType::kTec, false, 0., 1U);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), -0.993747, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0.1116511, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_tec2) {
  JonesParameters jones = PrepareJonesParameters(GainType::kTec, true, 0., 2U);
  const auto parms = jones.GetParms();

  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), -0.993747, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -0.1116511, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_r_angle) {
  JonesParameters jones = PrepareJonesParameters(GainType::kRotationAngle);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 4);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 0.487187, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0., 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_r_angle_inverted) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kRotationAngle, true);
  const auto parms = jones.GetParms();

  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 0.487187, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0.0, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_r_measure) {
  JonesParameters jones = PrepareJonesParameters(GainType::kRotationMeasure);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 4);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), -0.185403, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0., 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_r_measure_inverted) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kRotationMeasure, true, 0.5, 4);
  const auto parms = jones.GetParms();

  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), -0.185403, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0.0, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_diagonal_phase) {
  JonesParameters jones = PrepareJonesParameters(GainType::kDiagonalPhase);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 0.487187, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -0.87329, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_scalar_phase) {
  JonesParameters jones = PrepareJonesParameters(GainType::kScalarPhase);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 0.487187, 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -0.87329, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_amplitude) {
  JonesParameters jones = PrepareJonesParameters(GainType::kDiagonalAmplitude);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 200., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0., 1e-3);
}

BOOST_AUTO_TEST_CASE(make_scalar_amplitude) {
  JonesParameters jones = PrepareJonesParameters(GainType::kScalarAmplitude);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 200., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 0., 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_clock) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kClock, false, 0., 1U);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 1., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 2.08822e-06, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_clock2) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kClock, true, 50, 2U);
  const auto parms = jones.GetParms();

  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 1., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), -2.08822e-06, 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_gain_re_im) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kDiagonalRealImaginary);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 200., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 200., 1e-3);
}

BOOST_AUTO_TEST_CASE(make_complex_fulljones_re_im) {
  JonesParameters jones =
      PrepareJonesParameters(GainType::kFullJonesRealImaginary);
  const auto parms = jones.GetParms();

  BOOST_CHECK_EQUAL(parms.shape()[0], 4);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).real(), 200., 1e-3);
  BOOST_CHECK_CLOSE(parms(0, 0, 0).imag(), 200., 1e-3);
}

BOOST_AUTO_TEST_CASE(fulljones_with_nullptr) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();

  BOOST_CHECK_THROW(
      JonesParameters jones(kFreqs, kTimes, antNames, GainType::kFullJones,
                            kInterpolationType, kDirection, &mock, nullptr),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(missing_antenna_error) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts + 1; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();

  BOOST_CHECK_THROW(JonesParameters jones(
                        kFreqs, kTimes, antNames, GainType::kDiagonalAmplitude,
                        kInterpolationType, kDirection, &mock, nullptr),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(missing_antenna_flag) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts + 1; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();

  JonesParameters jones(kFreqs, kTimes, antNames, GainType::kDiagonalAmplitude,
                        kInterpolationType, kDirection, &mock, nullptr, false,
                        0, 0, JonesParameters::MissingAntennaBehavior::kFlag);

  const auto parms = jones.GetParms();
  BOOST_CHECK_EQUAL(parms.shape()[1], kNAnts + 1);
  BOOST_CHECK(std::isfinite(parms(0, kNAnts - 1, 0).real()));
  BOOST_CHECK(!std::isfinite(parms(0, kNAnts, 0).real()));
}

BOOST_AUTO_TEST_CASE(missing_antenna_unit_diag) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts + 1; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();

  JonesParameters jones(kFreqs, kTimes, antNames, GainType::kDiagonalAmplitude,
                        kInterpolationType, kDirection, &mock, nullptr, false,
                        0, 0, JonesParameters::MissingAntennaBehavior::kUnit);

  const auto parms = jones.GetParms();
  BOOST_CHECK_EQUAL(parms.shape()[0], 2);
  BOOST_CHECK_EQUAL(parms.shape()[1], kNAnts + 1);
  BOOST_CHECK_EQUAL(parms(0, kNAnts, 0).real(), 1.);
  BOOST_CHECK_EQUAL(parms(1, kNAnts, 0).real(), 1.);
  BOOST_CHECK_EQUAL(parms(0, kNAnts, 0).imag(), 0.);
}

BOOST_AUTO_TEST_CASE(missing_antenna_unit_full) {
  std::vector<std::string> antNames;
  for (unsigned int i = 0; i < kNAnts + 1; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
  }

  SolTabMock mock = SolTabMock();

  JonesParameters jones(kFreqs, kTimes, antNames, GainType::kRotationAngle,
                        kInterpolationType, kDirection, &mock, nullptr, false,
                        0, 0, JonesParameters::MissingAntennaBehavior::kUnit);

  const auto parms = jones.GetParms();
  BOOST_CHECK_EQUAL(parms.shape()[0], 4);
  BOOST_CHECK_EQUAL(parms.shape()[1], kNAnts + 1);
  BOOST_CHECK_EQUAL(parms(0, kNAnts, 0).real(), 1.);
  BOOST_CHECK_EQUAL(parms(0, kNAnts, 0).imag(), 0.);
  BOOST_CHECK_EQUAL(parms(1, kNAnts, 0).real(), 0.);
  BOOST_CHECK_EQUAL(parms(1, kNAnts, 0).imag(), 0.);
  BOOST_CHECK_EQUAL(parms(2, kNAnts, 0).real(), 0.);
  BOOST_CHECK_EQUAL(parms(2, kNAnts, 0).imag(), 0.);
  BOOST_CHECK_EQUAL(parms(3, kNAnts, 0).real(), 1.);
  BOOST_CHECK_EQUAL(parms(3, kNAnts, 0).imag(), 0.);
}

BOOST_AUTO_TEST_CASE(gain_type_to_human_readable_str) {
  const std::vector<std::string> kGainTypeNames = {
      "diagonal complex",
      "full-Jones",
      "scalar complex",
      "TEC",
      "clock",
      "rotation angle",
      "scalar phase",
      "diagonal phase",
      "rotation measure",
      "scalar amplitude",
      "diagonal amplitude",
      "diagonal complex (real/imaginary)",
      "full-Jones (real/imaginary)"};
  const std::vector<GainType> kGainTypes = {GainType::kDiagonalComplex,
                                            GainType::kFullJones,
                                            GainType::kScalarComplex,
                                            GainType::kTec,
                                            GainType::kClock,
                                            GainType::kRotationAngle,
                                            GainType::kScalarPhase,
                                            GainType::kDiagonalPhase,
                                            GainType::kRotationMeasure,
                                            GainType::kScalarAmplitude,
                                            GainType::kDiagonalAmplitude,
                                            GainType::kDiagonalRealImaginary,
                                            GainType::kFullJonesRealImaginary};
  BOOST_REQUIRE_EQUAL(kGainTypeNames.size(), kGainTypes.size());
  for (size_t i = 0; i != kGainTypeNames.size(); ++i) {
    const std::string human_string =
        JonesParameters::GainTypeToHumanReadableString(kGainTypes[i]);
    BOOST_CHECK_EQUAL(kGainTypeNames[i], human_string);
  }
}

BOOST_AUTO_TEST_CASE(h5_parm_type_string_to_gain_type) {
  const std::vector<std::string> kH5TypeStrings = {
      "gain",        "fulljones",       "tec",           "clock",
      "scalarphase", "scalaramplitude", "rotationangle", "rotationmeasure",
      "phase",       "amplitude"};
  const std::vector<GainType> kGainTypes = {GainType::kDiagonalComplex,
                                            GainType::kFullJones,
                                            GainType::kTec,
                                            GainType::kClock,
                                            GainType::kScalarPhase,
                                            GainType::kScalarAmplitude,
                                            GainType::kRotationAngle,
                                            GainType::kRotationMeasure,
                                            GainType::kDiagonalPhase,
                                            GainType::kDiagonalAmplitude};
  BOOST_REQUIRE_EQUAL(kH5TypeStrings.size(), kGainTypes.size());
  for (size_t i = 0; i != kH5TypeStrings.size(); ++i) {
    const GainType result_gain_type =
        JonesParameters::H5ParmTypeStringToGainType(kH5TypeStrings[i]);
    BOOST_CHECK(result_gain_type == kGainTypes[i]);
  }
  BOOST_CHECK_THROW(
      JonesParameters::H5ParmTypeStringToGainType("nosuchcorrection"),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
