// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_H5PARM_JONESPARAMETERS_H_
#define SCHAAPCOMMON_H5PARM_JONESPARAMETERS_H_

#include <complex>
#include <vector>

#include "h5parm.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/String.h>

namespace schaapcommon {
namespace h5parm {

/**
 * Type of Jones matrix.
 * NOTE: kScalarPhase, kScalarAmplitude, kRealImaginary and
 * kFullJonesRealImaginary are added to be compatible with ParmDB.
 */
enum class GainType {
  kDiagonalComplex,
  kFullJones,
  kScalarComplex,
  kTec,
  kClock,
  kRotationAngle,
  kScalarPhase,
  kDiagonalPhase,
  kRotationMeasure,
  kScalarAmplitude,
  kDiagonalAmplitude,
  kDiagonalRealImaginary,
  kFullJonesRealImaginary
};

/// @brief Class to extract Jones matrices from an h5parm.
/// Provides some compatibility with ParmDB.
class JonesParameters {
 public:
  /**
   * What to do with missing antennas
   */
  enum class MissingAntennaBehavior {
    kError,  ///< Raise an error on missing antennas
    kFlag,   ///< Insert flagged parameter values for missing antennas
    kUnit    ///< Insert a unit Jones matrix for missing antennas
  };

  /**
   * Interpolation in time and frequency.
   */
  enum class InterpolationType { NEAREST, LINEAR };

  /**
   * Constructor for JonesParameters with given parm_values. To be used if
   * parameter values are read externally (e.g. from a ParmDB)
   * \param freqs Output frequency for sampled values
   * \param times Output times for sampled values
   * \param antenna_names Names of the antennas
   * \param correct_type Correction type of the Jones matrix
   * \param interpolation_type Interpolation type of the Jones matrix
   * \param direction Direction number in the H5parm
   * \param parm_values Parameter values (e.g. TEC values)
   * \param invert (optional default=false) Invert the parameters
   * \param sigma_mmse (optional default=0.) Minimum mean square error parameter
   * (remnant from BBS times, leave at 0 unless you know what you're doing).
   */
  JonesParameters(const std::vector<double>& freqs,
                  const std::vector<double>& times,
                  const std::vector<std::string>& antenna_names,
                  GainType gain_type, InterpolationType interpolation_type,
                  hsize_t direction,
                  std::vector<std::vector<std::vector<double>>>&& parm_values,
                  bool invert = false, float sigma_mmse = 0.);

  /**
   * Constructor for JonesParameters with given parm_values. To be used if
   * solutions from prior steps are already in buffer. Allows the immediate
   * application of solutions by passing through the buffer.
   * \param freqs Output frequency for sampled values
   * \param times Output times for sampled values
   * \param antenna_names Names of the antennas
   * \param correct_type Correction type of the Jones matrix
   * \param solution Solution in format [n_ants * n_pols, n_chans]
   * \param invert (optional default=false) Invert the parameters
   * \param sigma_mmse (optional default=0.) Minimum mean square error parameter
   * (remnant from BBS times, leave at 0 unless you know what you're doing).
   */
  JonesParameters(
      const std::vector<double>& freqs, const std::vector<double>& times,
      const std::vector<std::string>& antenna_names, GainType gain_type,
      const std::vector<std::vector<std::complex<double>>>& solution,
      bool invert = false, float sigma_mmse = 0.);

  /**
   * Contructor for JonesParameters. JonesParameters will extract parameter
   * values itself from an H5Parm.
   * \param freqs Output frequency for sampled values
   * \param times Output times for sampled values
   * \param antenna_names Names of the antennas
   * \param correct_type Correction type of the Jones matrix
   * \param interpolation_type Interpolation type of the Jones matrix
   * \param direction Direction number in the H5parm
   * \param sol_tab soltab with parameters
   * \param sol_tab2 (optional default=nullptr) soltab with parameters for
   * complex values. Shapes of sol_tab and sol_tab2 can differ
   * \param invert (optional default=false) Invert the parameters
   * \param sigma_mmse (optional default=0.) Minimum mean square error parameter
   * (remnant from BBS times, leave at 0 unless you know what you're doing)
   * \param parm_size (optional default=0) allows to override the vector size
   * for parm_values
   * \param missing_antenna_behavior (optional default=kError) what to do with
   * missing antennas
   */
  JonesParameters(const std::vector<double>& freqs,
                  const std::vector<double>& times,
                  const std::vector<std::string>& antenna_names,
                  GainType gain_type, InterpolationType interpolation_type,
                  hsize_t direction, schaapcommon::h5parm::SolTab* sol_tab,
                  schaapcommon::h5parm::SolTab* sol_tab2 = nullptr,
                  bool invert = false, float sigma_mmse = 0.,
                  unsigned int parm_size = 0,
                  MissingAntennaBehavior missing_antenna_behavior =
                      MissingAntennaBehavior::kError);

  /**
   * Return the Jones matrices as a casacore cube with dimensions (nparms,
   * nantenna, ntime*nfreq), frequency varies fastest. nparms is 2 for diagonal,
   * 4 for full jones parameters.
   */
  const casacore::Cube<std::complex<float>>& GetParms() const { return parms_; }

  /**
   * Parse a H5Parm type string into an GainType enum value.
   * The H5 strings are not equal to the corresponding strings
   * that are returned by @ref GainTypeToHumanReadableString().
   * @throws std::runtime_error when string is unrecognized
   */
  static GainType H5ParmTypeStringToGainType(
      const std::string& h5parm_type_string);

  /**
   * Convert GainType to a string suitable for outputting to the user.
   */
  static std::string GainTypeToHumanReadableString(GainType);

  /**
   * Parse a missing antennabehavior string into an enum value
   */
  static MissingAntennaBehavior StringToMissingAntennaBehavior(
      const std::string&);

  /**
   * Convert MissingAntennaBehavior enum to string
   */
  static std::string MissingAntennaBehaviorToString(MissingAntennaBehavior);

 private:
  /**
   * Fill parms_ with the Jones matrices that correspond to parameters in
   * parmvalues_. Inverts the Jones matrix if invert is true.
   * \param ant Antenna number
   * \param invert (optional default=false) Invert the parameters. This will
   * ONLY have an effect on RotationMeasure and RotationAngle. Other effects
   * have to be inverted explicitly by calling Invert()
   */
  void MakeComplex(size_t ant, const std::vector<double>& freqs,
                   GainType correct_type, bool invert = false);

  /**
   * Get the number of parameters for a given @param correct_type.
   */
  static unsigned int GetNParms(GainType correct_type);

  /**
   * Get the dimension for parm_values, i.e. the number of parameter names in
   * the H5Parm.
   */
  static unsigned int GetNParmValues(GainType correct_type);

  /**
   * Fill the JonesParameters parameter values from the solution tables
   * \param sol_tab soltab with parameters
   * \param sol_tab2 (optional) soltab with parameters for complex values.
   * Shapes of sol_tab and sol_tab2 can differ
   * \param freqs Output frequency for sampled values
   * \param times Output times for sampled values
   * \param antenna_names Names of the antennas
   * \param ant Antenna number
   * \param correct_type Correction type of the Jones matrix
   * \param interpolation_type Interpolation type of the Jones matrix
   */
  void FillParmValues(schaapcommon::h5parm::SolTab* sol_tab,
                      schaapcommon::h5parm::SolTab* sol_tab2,
                      const std::vector<double>& freqs,
                      const std::vector<double>& times,
                      const std::vector<std::string>& antenna_names, size_t ant,
                      GainType gain_type, InterpolationType interpolation_type,
                      hsize_t direction);

  /**
   * Replace values by NaN on places where weight is zero
   */
  static void ApplyFlags(std::vector<double>& values,
                         const std::vector<double>& weights);

  /**
   *  Static function to invert the complex parameters. Is automatically called
   * in MakeComplex.
   * \param parms Reference to the complex parameters that will be inverted
   * obtained via MakeComplex
   * \param sigma_mmse (optional default=0.) Minimum mean square error parameter
   * (remnant from BBS times, leave at 0 unless you know what you're doing).
   * \param correct_type Correction type of the Jones matrix
   */
  static void Invert(casacore::Cube<std::complex<float>>& parms,
                     float sigma_mmse, GainType gain_type);

  /// Parameter values, inner vector has dimension num_times * num_frequencies,
  /// the middle vector has dimension number_antennas, outer vector has
  /// dimension num_parameters (e.g. phase and amplitude)
  /// These are the paramters as they are stored (e.g. TEC values)
  /// TODO: This probably does not need to be a member
  std::vector<std::vector<std::vector<double>>> parm_values_;
  /// Stored Jones matrices, dimensions (nparms, nantenna, ntime*nfreq),
  /// frequency varies fastest. nparms is 2 for diagonal, 4 for full jones
  /// parameters.
  casacore::Cube<std::complex<float>> parms_;
};

}  // namespace h5parm
}  // namespace schaapcommon

#endif
