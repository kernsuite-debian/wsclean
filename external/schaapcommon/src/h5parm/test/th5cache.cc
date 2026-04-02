// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "h5cache.h"

#include <boost/test/unit_test.hpp>

#include <filesystem>

#include <H5Cpp.h>

using schaapcommon::h5parm::H5Cache;

namespace {
const std::string kFilename = "test_h5cache.h5";
const std::string kGroupName = "sol000";
const std::string kTypeName = "type_name";
const std::string kAxisNames = "time,freq,ant,dir";

// This fixture creates a minimally valid h5 file for H5Cache and removes it
// upon destruction.
class H5Fixture {
 public:
  H5Fixture()
      : file_(kFilename, H5F_ACC_TRUNC), group_(file_.createGroup(kGroupName)) {
    // Add 'TITLE'
    const H5::DataType type_datatype(H5T_STRING, kTypeName.size());
    H5::Attribute type_attribute =
        group_.createAttribute("TITLE", type_datatype, H5::DataSpace());
    type_attribute.write(type_datatype, kTypeName);

    // Dimension of each axis
    constexpr size_t kNAntennas = 5;
    constexpr size_t kNDirections = 8;
    constexpr size_t kNFrequencies = 12;
    constexpr size_t kNTimes = 50;

    // Generating 'values' and 'weights' for testing.
    size_t total = kNAntennas * kNDirections * kNFrequencies * kNTimes;
    std::vector<double> values(total);
    std::vector<double> weights(total);
    for (size_t i = 0; i < total; ++i) {
      // Making each value have a correlative intger.
      values[i] = i;
      // Making the weights have zero for every seventh value.
      // Testing can verify that the value becomes NaN for these indices.
      weights[i] = 1;
      if ((i % 7) == 0) {
        weights[i] = 0;
      }
    }

    // Create the 'ant' dataset.
    size_t str_max_length = 9;
    H5::DataType ant_datatype = H5::StrType(H5::PredType::C_S1, str_max_length);
    H5::DataSet ant_dataset = CreateDataset("ant", ant_datatype, kNAntennas);

    // Write some data in 'ant' dataset.
    const std::vector<std::string>& sol_antennas = {
        "CS001HBA0", "CS001HBA1", "CS002HBA0", "CS002HBA1", "CS003HBA0"};
    WriteData(ant_dataset, ant_datatype, sol_antennas, str_max_length);

    // Create the 'dir' dataset.
    H5::DataType dir_datatype = H5::StrType(H5::PredType::C_S1, 11);
    CreateDataset("dir", dir_datatype, kNDirections);

    // Create the 'freq' dataset.
    H5::DataType freq_datatype = H5::PredType::NATIVE_DOUBLE;
    H5::DataSet freq_dataset =
        CreateDataset("freq", freq_datatype, kNFrequencies);

    // These time and frequency values come from actual observation
    // data. They are a tiny sample of the file used to reproduce
    // the problem that triggered the need to use cached values.

    // Write data in 'freq' dataset.
    std::vector<double> freq_data = {1.207e+8, 1.217e+8, 1.227e+8, 1.236e+8,
                                     1.246e+8, 1.256e+8, 1.266e+8, 1.275e+8,
                                     1.285e+8, 1.295e+8, 1.305e+8, 1.314e+8};
    freq_dataset.write(freq_data.data(), H5::PredType::IEEE_F64LE);

    // Create the 'time' dataset.
    H5::DataType time_datatype = H5::PredType::NATIVE_DOUBLE;
    H5::DataSet time_dataset = CreateDataset("time", time_datatype, kNTimes);

    // Write data in 'time' dataset.
    std::vector<double> time_data = {
        5034516928.039, 5034516936.050, 5034516944.061, 5034516952.072,
        5034516960.083, 5034516968.095, 5034516976.106, 5034516984.117,
        5034516992.128, 5034517000.139, 5034517008.150, 5034517016.161,
        5034517024.172, 5034517032.183, 5034517040.195, 5034517048.206,
        5034517056.217, 5034517064.228, 5034517072.239, 5034517080.250,
        5034517088.261, 5034517096.272, 5034517104.284, 5034517112.295,
        5034517120.306, 5034517128.317, 5034517136.328, 5034517144.339,
        5034517152.350, 5034517160.361, 5034517168.373, 5034517176.384,
        5034517184.395, 5034517192.406, 5034517200.417, 5034517208.428,
        5034517216.439, 5034517224.450, 5034517232.462, 5034517240.473,
        5034517248.484, 5034517256.495, 5034517264.506, 5034517272.517,
        5034517280.528, 5034517288.539, 5034517296.550, 5034517304.562,
        5034517312.573, 5034517320.584};
    time_dataset.write(time_data.data(), H5::PredType::IEEE_F64LE);

    // Data dimensions for 4-dimensional arrays: values and weights.
    const std::array<hsize_t, 4> set_dimensions{
        {kNTimes, kNFrequencies, kNAntennas, kNDirections}};

    // Create 'val' dataset.
    H5::DataSpace val_space(set_dimensions.size(), set_dimensions.data());
    H5::DataSet val =
        group_.createDataSet("val", H5::PredType::NATIVE_DOUBLE, val_space);

    // Write data in 'val' dataset.
    val.write(values.data(), H5::PredType::IEEE_F64LE);

    // Add 'AXES' attribute to 'val' dataset.
    H5::DataType axes_datatype(H5T_STRING, kAxisNames.size());
    H5::Attribute axes_attribute =
        val.createAttribute("AXES", axes_datatype, H5::DataSpace());
    axes_attribute.write(axes_datatype, kAxisNames);

    // Create 'weight' dataset.
    H5::DataSpace weight_space(set_dimensions.size(), set_dimensions.data());
    H5::DataSet weight =
        group_.createDataSet("weight", H5::PredType::NATIVE_DOUBLE, val_space);

    // Write data in 'weights' dataset.
    weight.write(weights.data(), H5::PredType::IEEE_F64LE);

    // Add 'AXES' attribute to 'weight' dataset.
    axes_attribute =
        weight.createAttribute("AXES", axes_datatype, H5::DataSpace());
    axes_attribute.write(axes_datatype, kAxisNames);
  }

  ~H5Fixture() { std::filesystem::remove(kFilename); }

  H5::Group& GetGroup() { return group_; };

 private:
  H5::DataSet CreateDataset(const std::string& name, H5::DataType datatype,
                            size_t size) {
    const std::array<hsize_t, 1> dimensions{{size}};
    H5::DataSpace space(dimensions.size(), dimensions.data());
    return group_.createDataSet(name, datatype, space);
  }

  void WriteData(const H5::DataSet& dataset, H5::DataType datatype,
                 const std::vector<std::string>& data, size_t max_length) {
    std::vector<char> output(data.size() * max_length, '\0');
    std::vector<char>::iterator it_out = output.begin();
    for (const std::string& text : data) {
      std::copy_n(text.begin(), text.size(), it_out);
      it_out += max_length;
    }
    dataset.write(output.data(), datatype);
  }

  H5::H5File file_;
  H5::Group group_;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(h5cache)

BOOST_FIXTURE_TEST_CASE(times_frequency_5x2, H5Fixture) {
  H5Cache h5cache(GetGroup());

  // These values of times and frequencies are chosen to fall
  // between the values of the fixture datasets. This approach
  // ensures comprehensive testing of all the methods involved,
  // including GridNearestNeighbor(..).
  const std::vector<double> times = {5034516940.0, 5034517000.0, 5034517100.0,
                                     5034517200.0, 5034517300.0};
  const std::vector<double> freqs = {1.220e+8, 1.305e+8};

  std::vector<double> flagged_values =
      h5cache.GetValues("CS002HBA1", times, freqs, 2, 0, true);
  BOOST_CHECK_EQUAL(flagged_values.size(), 10);
  BOOST_CHECK_CLOSE(flagged_values[0], 544, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[1], 904, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[2], 4384, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[3], 4744, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[4], 10144, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[5], 10504, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[6], 16384, 1e-8);
  BOOST_CHECK(std::isnan(flagged_values[7]));
  BOOST_CHECK_CLOSE(flagged_values[8], 22144, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[9], 22504, 1e-8);

  flagged_values = h5cache.GetValues("CS003HBA0", times, freqs, 2, 0, true);
  BOOST_CHECK_EQUAL(flagged_values.size(), 10);
  BOOST_CHECK_CLOSE(flagged_values[0], 552, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[1], 912, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[2], 4392, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[3], 4752, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[4], 10152, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[5], 10512, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[6], 16392, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[7], 16752, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[8], 22152, 1e-8);
  BOOST_CHECK(std::isnan(flagged_values[9]));
}

BOOST_FIXTURE_TEST_CASE(times_frequency_4x3, H5Fixture) {
  H5Cache h5cache(GetGroup());

  // These values of times and frequencies are chosen to fall
  // between the values of the fixture datasets. This approach
  // ensures comprehensive testing of all the methods involved,
  // including GridNearestNeighbor(..).
  const std::vector<double> times = {5034516900.0, 5034517150.0, 5034517250.0,
                                     5034517350.0};
  const std::vector<double> freqs = {1.250e+8, 1.280e+8, 1.300e+8};

  std::vector<double> flagged_values =
      h5cache.GetValues("CS001HBA0", times, freqs, 2, 0, true);
  BOOST_CHECK_EQUAL(flagged_values.size(), 12);
  BOOST_CHECK_CLOSE(flagged_values[0], 160, 1e-8);
  BOOST_CHECK(std::isnan(flagged_values[1]));
  BOOST_CHECK_CLOSE(flagged_values[2], 360, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[3], 13600, 1e-8);
  BOOST_CHECK(std::isnan(flagged_values[4]));
  BOOST_CHECK_CLOSE(flagged_values[5], 13800, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[6], 19360, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[7], 19480, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[8], 19560, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[9], 23680, 1e-8);
  BOOST_CHECK(std::isnan(flagged_values[10]));
  BOOST_CHECK_CLOSE(flagged_values[11], 23880, 1e-8);

  flagged_values = h5cache.GetValues("CS002HBA1", times, freqs, 2, 0, true);
  BOOST_CHECK_EQUAL(flagged_values.size(), 12);
  BOOST_CHECK_CLOSE(flagged_values[0], 184, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[1], 304, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[2], 384, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[3], 13624, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[4], 13744, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[5], 13824, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[6], 19384, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[7], 19504, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[8], 19584, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[9], 23704, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[10], 23824, 1e-8);
  BOOST_CHECK_CLOSE(flagged_values[11], 23904, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
