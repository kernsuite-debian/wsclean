// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "soltab.h"

#include <boost/test/unit_test.hpp>

#include <filesystem>

#include <H5Cpp.h>

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::SolTab;

namespace {
const std::string kFilename = "test_soltab.h5";
const std::string kGroupName = "group_name";
const std::string kTypeName = "type_name";
const std::string kAxisNames = "axis_name_0,axis_name_1";

// This fixture creates a minimally valid h5 file for SolTab and removes it upon
// destruction.
struct H5Fixture {
  H5Fixture()
      : file(kFilename, H5F_ACC_TRUNC), group(file.createGroup(kGroupName)) {
    // Add type ("TITLE").
    const H5::DataType type_datatype(H5T_STRING, kTypeName.size());
    H5::Attribute type_attribute =
        group.createAttribute("TITLE", type_datatype, H5::DataSpace());
    type_attribute.write(type_datatype, kTypeName);

    // Add "val" dataset.
    const std::array<hsize_t, 2> val_dimensions{{42, 43}};
    const H5::DataSpace val_space(val_dimensions.size(), val_dimensions.data());
    H5::DataSet val =
        group.createDataSet("val", H5::PredType::NATIVE_DOUBLE, val_space);

    // Add "AXES" to "val" dataset.
    const H5::DataType axes_datatype(H5T_STRING, kAxisNames.size());
    H5::Attribute axes_attribute =
        val.createAttribute("AXES", axes_datatype, H5::DataSpace());
    axes_attribute.write(axes_datatype, kAxisNames);
  }

  ~H5Fixture() { std::filesystem::remove(kFilename); }

  H5::H5File file;
  H5::Group group;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(soltab)

BOOST_AUTO_TEST_CASE(constructor) {
  const SolTab soltab;

  BOOST_TEST(soltab.GetAxes().size() == 0);
  BOOST_TEST(soltab.NumAxes() == 0);
  BOOST_TEST(!soltab.HasAxis("foo"));
  BOOST_CHECK_THROW(soltab.GetAxis("foo"), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetAxisIndex("foo"), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetRealAxis("foo"), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetStringAxis("foo"), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetFreqIndex(0.0), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetTimeIndex(0.0), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetDirIndex("foo"), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetTimeInterval(), std::runtime_error);
  BOOST_CHECK_THROW(soltab.GetFreqInterval(), std::runtime_error);

  BOOST_TEST(soltab.GetName() == "<invalid>");
  BOOST_TEST(soltab.GetType() == "");

  BOOST_CHECK_THROW(soltab.GetValues("ant", 0, 0, 0, 0, 0, 0, 0, 0),
                    H5::Exception);

  BOOST_CHECK_THROW(soltab.GetWeights("ant", 0, 0, 0, 0, 0, 0, 0, 0),
                    H5::Exception);

  const std::vector<double> times{0.0};
  const std::vector<double> frequencies{42.0};
  BOOST_CHECK_THROW(
      soltab.GetValuesOrWeights("foo", "ant", times, frequencies, 0, 0, false),
      H5::Exception);
}

BOOST_FIXTURE_TEST_CASE(construct_from_group, H5Fixture) {
  SolTab soltab(group);

  BOOST_TEST(soltab.GetName() == kGroupName);
  BOOST_TEST(soltab.GetType() == kTypeName);
}

BOOST_FIXTURE_TEST_CASE(construct_no_type, H5Fixture) {
  group.removeAttr("TITLE");
  BOOST_CHECK_THROW(SolTab soltab(group), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(construct_no_values, H5Fixture) {
  group.unlink("val");
  BOOST_CHECK_THROW(SolTab soltab(group), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(construct_no_axes, H5Fixture) {
  H5::DataSet val = group.openDataSet("val");
  val.removeAttr("AXES");
  BOOST_CHECK_THROW(SolTab soltab(group), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(construct_single_axis, H5Fixture) {
  // Test that SolTab throws an exception if the number of axes does not
  // match the number of dimensions in the value space.

  const std::string kSingleAxisName = "single_axis";

  H5::DataSet val = group.openDataSet("val");
  H5::Attribute axes_attribute = val.openAttribute("AXES");
  const H5::DataType axes_datatype(H5T_STRING, kSingleAxisName.size());
  axes_attribute.write(axes_datatype, kSingleAxisName);

  BOOST_CHECK_THROW(SolTab soltab(group), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(dir_ant_order, H5Fixture) {
  // Tests that GetStringAxis() returns the values in their original order.
  // Also tests that GetAntIndex behave correctly.
  // The values for "dir" are ordered, the values for "ant" are not ordered.
  const std::vector<std::string> kDirValues{"d0", "d1", "d2"};
  const std::vector<std::string> kAntValues{"a2", "a3", "a1", "a0"};
  const std::size_t kMaxStringLength = 2;
  const H5::DataType kDataType(H5T_STRING, kMaxStringLength);

  const std::array<hsize_t, 1> dir_dimensions{{kDirValues.size()}};
  const std::array<hsize_t, 1> ant_dimensions{{kAntValues.size()}};
  const H5::DataSpace dir_space(dir_dimensions.size(), dir_dimensions.data());
  const H5::DataSpace ant_space(ant_dimensions.size(), ant_dimensions.data());
  H5::DataSet dir_set = group.createDataSet("dir", kDataType, dir_space);
  H5::DataSet ant_set = group.createDataSet("ant", kDataType, ant_space);
  dir_set.write("d0d1d2", kDataType);
  ant_set.write("a2a3a1a0", kDataType);

  const SolTab soltab(group);

  const std::vector<std::string>& dir_values = soltab.GetStringAxis("dir");
  BOOST_CHECK_THROW(soltab.GetStringAxis("foo"), std::runtime_error);
  const std::vector<std::string>& ant_values = soltab.GetStringAxis("ant");

  BOOST_CHECK_EQUAL_COLLECTIONS(dir_values.begin(), dir_values.end(),
                                kDirValues.begin(), kDirValues.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(ant_values.begin(), ant_values.end(),
                                kAntValues.begin(), kAntValues.end());

  for (std::size_t i = 0; i < kDirValues.size(); ++i) {
    BOOST_TEST(soltab.GetDirIndex(kDirValues[i]) == i);
  }
  BOOST_CHECK_THROW(soltab.GetDirIndex("invalid_dir"), std::runtime_error);

  for (std::size_t i = 0; i < kAntValues.size(); ++i) {
    BOOST_TEST(soltab.GetAntIndex(kAntValues[i]) == i);
  }
  BOOST_CHECK_THROW(soltab.GetAntIndex("invalid_ant"), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(set_antennas, H5Fixture) {
  const std::vector<std::string> kNames{"foo", "bar", "long-antenna-name"};
  const std::vector<std::string> kNewNames{"new_name"};

  SolTab soltab(group);
  // Table "ant" does not exist yet
  BOOST_CHECK_THROW(soltab.GetStringAxis("ant"), std::runtime_error);

  soltab.SetAntennas(kNames);
  const std::vector<std::string>& names = soltab.GetStringAxis("ant");
  BOOST_CHECK_EQUAL_COLLECTIONS(kNames.begin(), kNames.end(), names.begin(),
                                names.end());
  for (std::size_t i = 0; i < kNames.size(); ++i) {
    BOOST_TEST(soltab.GetAntIndex(kNames[i]) == i);
  }

  soltab.SetAntennas(kNewNames);
  const std::vector<std::string>& new_names = soltab.GetStringAxis("ant");
  BOOST_CHECK_EQUAL_COLLECTIONS(kNewNames.begin(), kNewNames.end(),
                                new_names.begin(), new_names.end());
  BOOST_TEST(soltab.GetAntIndex(kNewNames[0]) == 0);
}

BOOST_AUTO_TEST_SUITE_END()
