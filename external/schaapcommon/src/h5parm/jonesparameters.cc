// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "jonesparameters.h"

#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/BasicSL/Constants.h>

#include <cassert>
#include <utility>

#include <aocommon/matrix2x2.h>

#include <boost/algorithm/string.hpp>

namespace schaapcommon {
namespace h5parm {
JonesParameters::JonesParameters(
    const std::vector<double>& freqs, const std::vector<double>& times,
    const std::vector<std::string>& antenna_names, GainType gain_type,
    InterpolationType interpolation_type, hsize_t direction,
    std::vector<std::vector<std::vector<double>>>&& parm_values, bool invert,
    float sigma_mmse)
    : parm_values_(std::move(parm_values)) {
  const unsigned int num_parms = GetNParms(gain_type);
  parms_.resize(num_parms, antenna_names.size(), freqs.size() * times.size());
  for (size_t ant = 0; ant < antenna_names.size(); ++ant) {
    MakeComplex(ant, freqs, gain_type, invert);
  }
  if (invert) {
    Invert(parms_, sigma_mmse, gain_type);
  }
}

JonesParameters::JonesParameters(
    const std::vector<double>& freqs, const std::vector<double>& times,
    const std::vector<std::string>& antenna_names, GainType gain_type,
    const std::vector<std::vector<std::complex<double>>>& solution, bool invert,
    float sigma_mmse) {
  /*
   * Solution is in format [n_ants * n_pols, n_chans]
   * Expected format of the parm_values is [n_chans, n_ants, n_pols]
   * So we reorder the matrices below
   */
  const unsigned int num_parms = GetNParms(gain_type);
  parms_.resize(num_parms, antenna_names.size(), freqs.size() * times.size());

  for (uint i = 0; i < solution.size(); i++) {
    uint ant = i / antenna_names.size();
    uint pol = i % antenna_names.size();

    for (uint chan = 0; chan < solution[i].size(); chan++) {
      parms_(chan, ant, pol) = solution[i][chan];
    }
  }

  if (invert) {
    Invert(parms_, sigma_mmse, gain_type);
  }
}

JonesParameters::JonesParameters(
    const std::vector<double>& freqs, const std::vector<double>& times,
    const std::vector<std::string>& antenna_names, GainType gain_type,
    InterpolationType interpolation_type, hsize_t direction,
    schaapcommon::h5parm::SolTab* sol_tab,
    schaapcommon::h5parm::SolTab* sol_tab2, bool invert, float sigma_mmse,
    unsigned int parm_size, MissingAntennaBehavior missing_antenna_behavior) {
  const unsigned int num_parms = GetNParms(gain_type);

  parms_.resize(num_parms, antenna_names.size(), freqs.size() * times.size());
  if (parm_size == 0U) {
    parm_size = GetNParmValues(gain_type);
  }
  parm_values_.resize(parm_size);
  for (auto& parm_values : parm_values_) {
    parm_values.resize(antenna_names.size());
  }
  for (size_t ant = 0; ant < antenna_names.size(); ++ant) {
    try {
      FillParmValues(sol_tab, sol_tab2, freqs, times, antenna_names, ant,
                     gain_type, interpolation_type, direction);
      MakeComplex(ant, freqs, gain_type, invert);
    } catch (const std::exception& e) {
      if (std::string(e.what()).rfind("SolTab has no element", 0) != 0) {
        throw;
      }
      if (missing_antenna_behavior == MissingAntennaBehavior::kError) {
        throw;
      } else if (missing_antenna_behavior == MissingAntennaBehavior::kFlag) {
        // Insert flagged solution
        const size_t tfsize = parms_.shape()[2];
        for (size_t parm = 0; parm < num_parms; ++parm) {
          for (size_t tf = 0; tf < tfsize; ++tf) {
            parms_(0, ant, tf) = std::numeric_limits<float>::quiet_NaN();
          }
        }
      } else if (missing_antenna_behavior == MissingAntennaBehavior::kUnit) {
        // Insert unit solution
        const size_t tfsize = parms_.shape()[2];
        if (num_parms == 2) {
          for (size_t tf = 0; tf < tfsize; ++tf) {
            parms_(0, ant, tf) = 1.;
            parms_(1, ant, tf) = 1.;
          }
        } else {
          for (size_t tf = 0; tf < tfsize; ++tf) {
            parms_(0, ant, tf) = 1.;
            parms_(1, ant, tf) = 0.;
            parms_(2, ant, tf) = 0.;
            parms_(3, ant, tf) = 1.;
          }
        }
      }
    }
  }
  if (invert) {
    Invert(parms_, sigma_mmse, gain_type);
  }
}

void JonesParameters::MakeComplex(size_t ant, const std::vector<double>& freqs,
                                  GainType gain_type, bool invert) {
  for (unsigned int tf = 0; tf < parms_.shape()[2]; ++tf) {
    const double freq = freqs[tf % freqs.size()];

    switch (gain_type) {
      case GainType::kDiagonalComplex:
        parms_(0, ant, tf) =
            casacore::polar(parm_values_[0][ant][tf], parm_values_[1][ant][tf]);
        parms_(1, ant, tf) =
            casacore::polar(parm_values_[2][ant][tf], parm_values_[3][ant][tf]);
        break;
      case GainType::kFullJones:
        parms_(0, ant, tf) =
            casacore::polar(parm_values_[0][ant][tf], parm_values_[1][ant][tf]);
        parms_(1, ant, tf) =
            casacore::polar(parm_values_[2][ant][tf], parm_values_[3][ant][tf]);
        parms_(2, ant, tf) =
            casacore::polar(parm_values_[4][ant][tf], parm_values_[5][ant][tf]);
        parms_(3, ant, tf) =
            casacore::polar(parm_values_[6][ant][tf], parm_values_[7][ant][tf]);
        break;
      case GainType::kScalarComplex:
        parms_(0, ant, tf) =
            casacore::polar(parm_values_[0][ant][tf], parm_values_[1][ant][tf]);
        parms_(1, ant, tf) = parms_(0, ant, tf);
        break;
      case GainType::kTec:
        parms_(0, ant, tf) = casacore::polar(
            1., parm_values_[0][ant][tf] * -8.44797245e9 / freq);
        if (parm_values_.size() == 1) {  // No TEC:0, only TEC:
          parms_(1, ant, tf) = casacore::polar(
              1., parm_values_[0][ant][tf] * -8.44797245e9 / freq);
        } else {  // TEC:0 and TEC:1
          parms_(1, ant, tf) = casacore::polar(
              1., parm_values_[1][ant][tf] * -8.44797245e9 / freq);
        }
        break;
      case GainType::kClock:
        parms_(0, ant, tf) = casacore::polar(
            1., parm_values_[0][ant][tf] * freq * casacore::C::_2pi);
        if (parm_values_.size() == 1) {  // No Clock:0, only Clock:
          parms_(1, ant, tf) = casacore::polar(
              1., parm_values_[0][ant][tf] * freq * casacore::C::_2pi);
        } else {  // Clock:0 and Clock:1
          parms_(1, ant, tf) = casacore::polar(
              1., parm_values_[1][ant][tf] * freq * casacore::C::_2pi);
        }
        break;
      case GainType::kRotationAngle: {
        double phi = parm_values_[0][ant][tf];
        if (invert) {
          phi = -phi;
        }
        const float sinv = std::sin(phi);
        const float cosv = std::cos(phi);
        parms_(0, ant, tf) = cosv;
        parms_(1, ant, tf) = -sinv;
        parms_(2, ant, tf) = sinv;
        parms_(3, ant, tf) = cosv;
      } break;
      case GainType::kRotationMeasure: {
        const double lambda2 =
            (casacore::C::c / freq) * (casacore::C::c / freq);
        double chi = parm_values_[0][ant][tf] * lambda2;
        if (invert) {
          chi = -chi;
        }
        const float sinv = std::sin(chi);
        const float cosv = std::cos(chi);
        parms_(0, ant, tf) = cosv;
        parms_(1, ant, tf) = -sinv;
        parms_(2, ant, tf) = sinv;
        parms_(3, ant, tf) = cosv;
      } break;
      case GainType::kDiagonalPhase:
      case GainType::kScalarPhase:
        parms_(0, ant, tf) = casacore::polar(1., parm_values_[0][ant][tf]);
        if (gain_type == GainType::kScalarPhase) {  // Same value for x and y
          parms_(1, ant, tf) = casacore::polar(1., parm_values_[0][ant][tf]);
        } else {  // Different value for x and y
          parms_(1, ant, tf) = casacore::polar(1., parm_values_[1][ant][tf]);
        }
        break;
      case GainType::kDiagonalAmplitude:
      case GainType::kScalarAmplitude:
        parms_(0, ant, tf) = parm_values_[0][ant][tf];
        if (gain_type ==
            GainType::kScalarAmplitude) {  // Same value for x and y
          parms_(1, ant, tf) = parm_values_[0][ant][tf];
        } else {  // Different value for x and y
          parms_(1, ant, tf) = parm_values_[1][ant][tf];
        }
        break;
      case GainType::kDiagonalRealImaginary:
        parms_(0, ant, tf) = std::complex<float>(parm_values_[0][ant][tf],
                                                 parm_values_[1][ant][tf]);
        parms_(1, ant, tf) = std::complex<float>(parm_values_[2][ant][tf],
                                                 parm_values_[3][ant][tf]);
        break;
      case GainType::kFullJonesRealImaginary:
        parms_(0, ant, tf) = std::complex<float>(parm_values_[0][ant][tf],
                                                 parm_values_[1][ant][tf]);
        parms_(1, ant, tf) = std::complex<float>(parm_values_[2][ant][tf],
                                                 parm_values_[3][ant][tf]);
        parms_(2, ant, tf) = std::complex<float>(parm_values_[4][ant][tf],
                                                 parm_values_[5][ant][tf]);
        parms_(3, ant, tf) = std::complex<float>(parm_values_[6][ant][tf],
                                                 parm_values_[7][ant][tf]);
    }
  }
}

unsigned int JonesParameters::GetNParms(GainType gain_type) {
  switch (gain_type) {
    case GainType::kFullJones:
    case GainType::kRotationAngle:
    case GainType::kRotationMeasure:
    case GainType::kFullJonesRealImaginary:
      return 4;
    default:
      return 2;
  }
}

unsigned int JonesParameters::GetNParmValues(GainType gain_type) {
  switch (gain_type) {
    case GainType::kFullJones:
    case GainType::kFullJonesRealImaginary:
      return 8;
    case GainType::kDiagonalComplex:
    case GainType::kDiagonalRealImaginary:
      return 4;
    case GainType::kDiagonalPhase:
    case GainType::kDiagonalAmplitude:
    case GainType::kScalarComplex:
      return 2;
    case GainType::kRotationAngle:
    case GainType::kScalarPhase:
    case GainType::kRotationMeasure:
    case GainType::kScalarAmplitude:
      return 1;
    case GainType::kTec:
    case GainType::kClock:
      throw std::runtime_error(
          "Gain type is variable. Use parameter parm_size instead.");
  }
  throw std::runtime_error("Unknown correction type");
}

void JonesParameters::FillParmValues(
    schaapcommon::h5parm::SolTab* sol_tab,
    schaapcommon::h5parm::SolTab* sol_tab2, const std::vector<double>& freqs,
    const std::vector<double>& times,
    const std::vector<std::string>& antenna_names, size_t ant,
    GainType gain_type, InterpolationType interpolation_type,
    hsize_t direction) {
  auto get_flagged_values = [times, freqs, direction, interpolation_type](
                                SolTab* sol_tab,
                                const std::string& antenna_name,
                                unsigned int pol) -> std::vector<double> {
    std::vector<double> parmvalues = sol_tab->GetValuesOrWeights(
        "val", antenna_name, times, freqs, pol, direction,
        interpolation_type == InterpolationType::NEAREST);
    std::vector<double> weights = sol_tab->GetValuesOrWeights(
        "weight", antenna_name, times, freqs, pol, direction,
        interpolation_type == InterpolationType::NEAREST);
    ApplyFlags(parmvalues, weights);
    return parmvalues;
  };

  const std::string& ant_name = antenna_names[ant];
  if (gain_type == GainType::kFullJones ||
      gain_type == GainType::kDiagonalComplex ||
      gain_type == GainType::kScalarComplex) {
    for (size_t pol = 0; pol < parm_values_.size() / 2; ++pol) {
      // Place amplitude in even and phase in odd elements
      parm_values_[pol * 2][ant] = get_flagged_values(sol_tab, ant_name, pol);
      if (sol_tab2 == nullptr) {
        throw std::runtime_error(
            "soltab2 cannot be a nullpointer for correct_type=FULLJONES and "
            "correct_type=GAIN");
      }
      parm_values_[pol * 2 + 1][ant] =
          get_flagged_values(sol_tab2, ant_name, pol);
    }
  } else {
    for (size_t pol = 0; pol < parm_values_.size(); ++pol) {
      parm_values_[pol][ant] = get_flagged_values(sol_tab, ant_name, pol);
    }
  }
}

void JonesParameters::ApplyFlags(std::vector<double>& values,
                                 const std::vector<double>& weights) {
  assert(values.size() == weights.size());
  auto weight_it = weights.cbegin();

  for (double& value : values) {
    if (*weight_it == 0.) {
      value = std::numeric_limits<double>::quiet_NaN();
    }
    ++weight_it;
  }
}

void JonesParameters::Invert(casacore::Cube<casacore::Complex>& parms,
                             float sigma_mmse, GainType gain_type) {
  for (unsigned int tf = 0; tf < parms.shape()[2]; ++tf) {
    for (unsigned int ant = 0; ant < parms.shape()[1]; ++ant) {
      if (parms.shape()[0] == 2) {
        parms(0, ant, tf) = 1.f / parms(0, ant, tf);
        parms(1, ant, tf) = 1.f / parms(1, ant, tf);
      } else if (gain_type == GainType::kFullJones ||
                 gain_type == GainType::kFullJonesRealImaginary) {
        aocommon::MC2x2F v(&parms(0, ant, tf));

        // Add the variance of the nuisance term to the elements on the
        // diagonal.
        const float variance = sigma_mmse * sigma_mmse;
        v[0] += variance;
        v[3] += variance;
        v.Invert();
        v.AssignTo(&parms(0, ant, tf));
      } else if (gain_type == GainType::kRotationMeasure ||
                 gain_type == GainType::kRotationAngle) {
        // rotationmeasure and commonrotationangle have been inverted already in
        // MakeComplex
      } else {
        throw std::runtime_error("Invalid correction type");
      }
    }
  }
}

std::string JonesParameters::GainTypeToHumanReadableString(GainType ct) {
  switch (ct) {
    case GainType::kDiagonalComplex:
      return "diagonal complex";
    case GainType::kDiagonalRealImaginary:
      return "diagonal complex (real/imaginary)";
    case GainType::kFullJones:
      return "full-Jones";
    case GainType::kFullJonesRealImaginary:
      return "full-Jones (real/imaginary)";
    case GainType::kTec:
      return "TEC";
    case GainType::kClock:
      return "clock";
    case GainType::kScalarComplex:
      return "scalar complex";
    case GainType::kScalarPhase:
      return "scalar phase";
    case GainType::kScalarAmplitude:
      return "scalar amplitude";
    case GainType::kRotationAngle:
      return "rotation angle";
    case GainType::kRotationMeasure:
      return "rotation measure";
    case GainType::kDiagonalPhase:
      return "diagonal phase";
    case GainType::kDiagonalAmplitude:
      return "diagonal amplitude";
  }
  throw std::runtime_error("Unknown correction type: " +
                           std::to_string(static_cast<int>(ct)));
}

GainType JonesParameters::H5ParmTypeStringToGainType(const std::string& str) {
  if (str == "gain") {
    return GainType::kDiagonalComplex;
  } else if (str == "fulljones") {
    return GainType::kFullJones;
  } else if (str == "tec") {
    return GainType::kTec;
  } else if (str == "clock") {
    return GainType::kClock;
  } else if (str == "scalargain" || str == "commonscalargain") {
    return GainType::kScalarComplex;
  } else if (str == "scalarphase" || str == "commonscalarphase") {
    return GainType::kScalarPhase;
  } else if (str == "scalaramplitude" || str == "commonscalaramplitude") {
    return GainType::kScalarAmplitude;
  } else if (str == "phase") {
    return GainType::kDiagonalPhase;
  } else if (str == "amplitude") {
    return GainType::kDiagonalAmplitude;
  } else if (str == "rotationangle" || str == "commonrotationangle" ||
             str == "rotation") {
    return GainType::kRotationAngle;
  } else if (str == "rotationmeasure") {
    return GainType::kRotationMeasure;
  } else {
    // TODO:
    // - FULLJONES_RE_IM
    // - GAIN_RE_IM
    // are not listed, ask TJ
    throw std::runtime_error("Unknown correction type: " + str);
  }
}

JonesParameters::MissingAntennaBehavior
JonesParameters::StringToMissingAntennaBehavior(
    const std::string& behavior_str) {
  if (behavior_str == "error") {
    return JonesParameters::MissingAntennaBehavior::kError;
  } else if (behavior_str == "flag") {
    return JonesParameters::MissingAntennaBehavior::kFlag;
  } else if (behavior_str == "unit") {
    return JonesParameters::MissingAntennaBehavior::kUnit;
  } else {
    throw std::runtime_error(
        "missingantennabehavior should be one of 'error', 'flag' or "
        "'unit', not '" +
        behavior_str + "'");
  }
}

std::string JonesParameters::MissingAntennaBehaviorToString(
    JonesParameters::MissingAntennaBehavior behavior) {
  if (behavior == JonesParameters::MissingAntennaBehavior::kError) {
    return "error";
  } else if (behavior == JonesParameters::MissingAntennaBehavior::kFlag) {
    return "flag";
  } else if (behavior == JonesParameters::MissingAntennaBehavior::kUnit) {
    return "unit";
  }
  throw std::logic_error("Unexpected missingantennabehavior");
}

}  // namespace h5parm
}  // namespace schaapcommon
