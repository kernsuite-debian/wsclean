// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "filewriter.h"
#include "reordering.h"
#include "aocommon/polarization.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/filesystem.hpp>

#include <stdexcept>
#include <vector>
#include <set>
#include <complex>
#include <cstddef>

using aocommon::Polarization;
using aocommon::PolarizationEnum;

namespace schaapcommon::reordering {

class FixtureDirectory {
 public:
  /// Create the temporary directory and set it as working directory
  FixtureDirectory() {
    boost::filesystem::create_directories(kPath);
    boost::filesystem::current_path(kPath);
  }

  FixtureDirectory(const FixtureDirectory&) = delete;
  FixtureDirectory& operator=(const FixtureDirectory&) = delete;

  /// Remove the temporary diectory
  /// Will always run
  ~FixtureDirectory() {
    boost::filesystem::current_path(kWorkDir);
    boost::filesystem::remove_all(kPath);
  }

 private:
  const boost::filesystem::path kPath = boost::filesystem::unique_path();
  const boost::filesystem::path kWorkDir = boost::filesystem::current_path();
};

using ComplexVector = std::vector<std::complex<float>>;

const ComplexVector kTestData{
    10.0f, 11.0f, 12.0f, 13.0f, 20.0f, 21.0f, 22.0f, 23.0f,
    30.0f, 31.0f, 32.0f, 33.0f, 40.0f, 41.0f, 42.0f, 43.0f,
};
const std::vector<float> kTestWeights{0.1f, 0.1f, 0.1f, 0.1f, 0.2f, 0.2f,
                                      0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.3f,
                                      0.4f, 0.4f, 0.4f, 0.4f};
// std::vector<bool> does not have a .data() method
const bool kTestFlags[] = {false, false, false, false, false, false,
                           false, false, true,  true,  true,  true,
                           false, false, false, false};

const std::vector<ChannelRange> kChannelRanges{{0, 0, 2}, {0, 6, 8}, {1, 0, 4}};
const schaapcommon::reordering::MSSelection kSelection;
const std::set<PolarizationEnum> kPolsOut{Polarization::StokesI,
                                          Polarization::StokesQ};
const aocommon::MultiBandData kBands;
const std::string kTemporaryDirectory = "tmp";
const double kStartTime = 1.0;

const HandleData kData("test.ms", "DATA", "MODEL_DATA",
                       StorageManagerType::Default, kTemporaryDirectory,
                       MakeRegularChannelMap(kChannelRanges), true, false,
                       kPolsOut, kSelection, {kBands}, 6, true,
                       [](HandleData) {});

const std::map<size_t, std::set<aocommon::PolarizationEnum>>
    kMsPolarizationsPerDataDescId{{0, {Polarization::XX, Polarization::YY}},
                                  {1, {Polarization::XX, Polarization::YY}}};

const std::vector<aocommon::OptionalNumber<size_t>> kDataDescIdPerPart{
    aocommon::OptionalNumber<size_t>(0), aocommon::OptionalNumber<size_t>(1)};

BOOST_AUTO_TEST_SUITE(reordered_filewriter)

BOOST_FIXTURE_TEST_CASE(file_creation, FixtureDirectory) {
  {
    boost::filesystem::create_directory(kTemporaryDirectory);
    HandleData data(kData);
    data.metadata_indices_ = MakeMetaFilesMap(data.channels_).first;
    FileWriter reordering_writer(data, kMsPolarizationsPerDataDescId,
                                 kDataDescIdPerPart, kStartTime);
  }

  const std::vector<std::string> expected_reorder_files{
      "test.ms-meta-0000.tmp",    "test.ms-meta-0001.tmp",
      "test.ms-part0000-I.tmp",   "test.ms-part0000-I-w.tmp",
      "test.ms-part0000-I-m.tmp", "test.ms-part0001-I.tmp",
      "test.ms-part0001-I-w.tmp", "test.ms-part0001-I-m.tmp",
      "test.ms-part0002-I.tmp",   "test.ms-part0002-I-w.tmp",
      "test.ms-part0002-I-m.tmp", "test.ms-part0000-Q-m.tmp",
      "test.ms-part0000-Q.tmp",   "test.ms-part0000-Q-w.tmp",
      "test.ms-part0001-Q-m.tmp", "test.ms-part0001-Q.tmp",
      "test.ms-part0001-Q-w.tmp", "test.ms-part0002-Q-m.tmp",
      "test.ms-part0002-Q.tmp",   "test.ms-part0002-Q-w.tmp"};

  for (const std::string& reorder_file : expected_reorder_files) {
    std::string file = kTemporaryDirectory + "/" + reorder_file;
    BOOST_CHECK(boost::filesystem::exists(file));
  }
}

BOOST_FIXTURE_TEST_CASE(file_content, FixtureDirectory) {
  {
    // Write content for spw0 and assert the contents of each file
    boost::filesystem::create_directory(kTemporaryDirectory);
    HandleData data(kData);
    data.metadata_indices_ = MakeMetaFilesMap(data.channels_).first;
    FileWriter reordering_writer(data, kMsPolarizationsPerDataDescId,
                                 kDataDescIdPerPart, kStartTime);

    reordering_writer.WriteMetaRow(0.1, 0.2, 0.3, 1.1, 0, 0, 1, 1);

    reordering_writer.WriteDataRow(kTestData.data(), kTestData.data(),
                                   kTestWeights.data(), kTestFlags, 0);

    reordering_writer.UpdateMetaHeaders();
    reordering_writer.UpdatePartHeaders(true);
  }

  // Assert the header of the meta file before updating the headers
  std::ifstream meta_file("tmp/test.ms-meta-0000.tmp");

  MetaHeader meta_header;
  meta_header.Read(meta_file);
  BOOST_CHECK(meta_file.good());
  BOOST_CHECK_CLOSE_FRACTION(meta_header.start_time, 1.0, 1e-5);
  BOOST_CHECK_EQUAL(meta_header.selected_row_count, 1);
  BOOST_CHECK_EQUAL(meta_header.filename_length, 7);

  std::vector<char> ms_path(meta_header.filename_length);
  meta_file.read(ms_path.data(), meta_header.filename_length);
  std::string ms_path_str(ms_path.begin(), ms_path.end());

  BOOST_CHECK_EQUAL(ms_path_str, "test.ms");

  // Assert meta row
  MetaRecord meta_record;
  meta_record.Read(meta_file);

  BOOST_CHECK_CLOSE_FRACTION(meta_record.u, 0.1, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.v, 0.2, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.w, 0.3, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.time, 1.1, 1e-5);
  BOOST_CHECK_EQUAL(meta_record.antenna1, 0);
  BOOST_CHECK_EQUAL(meta_record.antenna2, 1);
  BOOST_CHECK_EQUAL(meta_record.field_id, 1);

  // Assert the contents of data
  std::ifstream data_file("tmp/test.ms-part0000-I.tmp");
  PartHeader part_header;
  part_header.Read(data_file);
  BOOST_CHECK(data_file.good());
  BOOST_CHECK_EQUAL(part_header.max_channel_count, 2);
  BOOST_CHECK(part_header.has_model);

  constexpr size_t kMaxChannels = 2;
  std::vector<std::complex<float>> data_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_data_buffer{{10.5f, 0.0f},
                                                        {12.5f, 0.0f}};

  data_file.read(reinterpret_cast<char*>(data_buffer.data()),
                 2 * sizeof(std::complex<float>));
  BOOST_CHECK(data_file.good());
  BOOST_CHECK_EQUAL_COLLECTIONS(data_buffer.begin(), data_buffer.end(),
                                expected_data_buffer.begin(),
                                expected_data_buffer.end());

  std::ifstream model_file("tmp/test.ms-part0000-I-m.tmp");
  std::vector<std::complex<float>> model_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_model_buffer{{10.5f, 0.0f},
                                                         {12.5f, 0.0f}};

  model_file.read(reinterpret_cast<char*>(model_buffer.data()),
                  2 * sizeof(std::complex<float>));
  BOOST_CHECK(model_file.good());
  BOOST_CHECK_EQUAL_COLLECTIONS(model_buffer.begin(), model_buffer.end(),
                                expected_model_buffer.begin(),
                                expected_model_buffer.end());

  std::ifstream weight_file("tmp/test.ms-part0000-I-w.tmp");
  std::vector<float> weight_buffer(kMaxChannels);
  std::vector<float> expected_weight_buffer{0.4f, 0.4f};

  weight_file.read(reinterpret_cast<char*>(weight_buffer.data()),
                   2 * sizeof(float));
  BOOST_CHECK(weight_file.good());
  BOOST_CHECK_EQUAL_COLLECTIONS(weight_buffer.begin(), weight_buffer.end(),
                                expected_weight_buffer.begin(),
                                expected_weight_buffer.end());
}

BOOST_FIXTURE_TEST_CASE(write_with_bda, FixtureDirectory) {
  boost::filesystem::create_directory(kTemporaryDirectory);

  constexpr size_t kDataDesc1 = 0;
  constexpr size_t kDataDesc2 = 2;
  constexpr size_t kFieldId = 5;

  std::map<size_t, std::set<aocommon::PolarizationEnum>> polarizations{
      {kDataDesc1, {Polarization::XX, Polarization::YY}},
      {kDataDesc2, {Polarization::XX, Polarization::YY}}};

  // data desc id 2: 3 channels, 2 polarizations (xx, yy)
  constexpr std::complex<float> kRow1Data[] = {100.0f, 102.0f, 101.0f,
                                               103.0f, 102.0f, 104.0f};
  // data desc id 1: 6 channels, 2 polarizations (xx, yy)
  constexpr std::complex<float> kRow2Data[] = {200.0f, 202.0f, 201.0f, 203.0f,
                                               202.0f, 204.0f, 203.0f, 205.0f,
                                               204.0f, 206.0f, 205.0f, 207.0f};
  constexpr std::complex<float> kRow3Data[] = {300.0f, 302.0f, 301.0f, 303.0f,
                                               302.0f, 304.0f, 303.0f, 305.0f,
                                               304.0f, 306.0f, 305.0f, 307.0f};
  constexpr std::complex<float> kRow4Data[] = {400.0f, 402.0f, 401.0f,
                                               403.0f, 402.0f, 404.0f};
  constexpr float kRowWeights[]{0.1f, 0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.4f,
                                0.4f, 0.5f, 0.5f, 0.6f, 0.6f, 0.7f, 0.7f};
  constexpr bool kRowFlagsDDI1[] = {true, true, false, false, false, false};
  constexpr bool kRowFlagsDDI2[] = {false, false, true,  true,  false, false,
                                    true,  true,  false, false, false, false};

  {
    // Write content for spw0 and assert the contents of each file
    const std::vector<aocommon::ChannelInfo> channels1{
        {100e6, 1e6}, {110e6, 1e6}, {120e6, 1e6},
        {130e6, 1e6}, {140e6, 1e6}, {150e6, 1e6}};
    aocommon::BandData band1(channels1, 100e6);
    const std::vector<aocommon::ChannelInfo> channels2{
        {105e6, 1e6}, {120e6, 1e6}, {145e6, 1e6}};
    aocommon::BandData band2(channels2, 60e6);

    aocommon::MultiBandData part_a_bands;
    part_a_bands.SetBand(kDataDesc1, aocommon::BandData(band1, 0, 3));
    part_a_bands.SetBand(kDataDesc2, aocommon::BandData(band2, 0, 2));
    aocommon::MultiBandData part_b_bands;
    part_b_bands.SetBand(kDataDesc1, aocommon::BandData(band1, 3, 6));
    part_b_bands.SetBand(kDataDesc2, aocommon::BandData(band2, 2, 3));

    const aocommon::VectorMap<ChannelRange> kChannelRangesPartA{
        {kDataDesc1, 0, 3}, {}, {kDataDesc2, 0, 2}};
    const aocommon::VectorMap<ChannelRange> kChannelRangesPartB{
        {kDataDesc1, 3, 6}, {}, {kDataDesc2, 2, 3}};
    const std::vector ranges_per_part{kChannelRangesPartA, kChannelRangesPartB};

    std::map<size_t, size_t> map;
    std::vector<aocommon::OptionalNumber<size_t>> data_desc_id_per_part;
    std::tie(map, data_desc_id_per_part) = MakeMetaFilesMap(ranges_per_part);
    BOOST_REQUIRE_EQUAL(data_desc_id_per_part.size(), 2);
    BOOST_REQUIRE_EQUAL(map.size(), 2);
    BOOST_CHECK(!data_desc_id_per_part[0]);
    BOOST_CHECK(!data_desc_id_per_part[1]);
    BOOST_CHECK_EQUAL(map[0], 0);
    BOOST_CHECK_EQUAL(map[1], 1);

    std::vector bands{part_a_bands, part_b_bands};
    HandleData bda_data("test.ms", "DATA", "MODEL_DATA",
                        StorageManagerType::Default, kTemporaryDirectory,
                        ranges_per_part, true, false, kPolsOut, kSelection,
                        bands, 10, true, [](HandleData) {});
    bda_data.metadata_indices_ = MakeMetaFilesMap(bda_data.channels_).first;
    FileWriter reordering_writer(bda_data, polarizations, data_desc_id_per_part,
                                 kStartTime);

    reordering_writer.WriteMetaRow(10.0, 11.0, 12.0, 1.1, kDataDesc2, 3, 4,
                                   kFieldId);
    reordering_writer.WriteMetaRow(13.0, 14.0, 15.0, 1.2, kDataDesc1, 4, 5,
                                   kFieldId);
    reordering_writer.WriteMetaRow(16.0, 17.0, 18.0, 1.3, kDataDesc1, 6, 7,
                                   kFieldId);
    reordering_writer.WriteMetaRow(-19.0, -20.0, -21.0, 1.4, kDataDesc2, 8, 9,
                                   kFieldId);

    reordering_writer.WriteDataRow(kRow1Data, kRow1Data, kRowWeights,
                                   kRowFlagsDDI1, kDataDesc2);
    reordering_writer.WriteDataRow(kRow2Data, kRow2Data, kRowWeights,
                                   kRowFlagsDDI2, kDataDesc1);
    reordering_writer.WriteDataRow(kRow3Data, kRow3Data, kRowWeights,
                                   kRowFlagsDDI2, kDataDesc1);
    reordering_writer.WriteDataRow(kRow4Data, kRow4Data, kRowWeights,
                                   kRowFlagsDDI1, kDataDesc2);

    reordering_writer.UpdateMetaHeaders();
    reordering_writer.UpdatePartHeaders(true);
  }

  auto check_metadata_header = [](std::ifstream& file, double start_time,
                                  size_t row_count) {
    MetaHeader meta_header;
    meta_header.Read(file);
    BOOST_CHECK(file.good());
    BOOST_CHECK_EQUAL(meta_header.start_time, start_time);
    BOOST_CHECK_EQUAL(meta_header.selected_row_count, row_count);
    BOOST_CHECK_EQUAL(meta_header.filename_length, 7);

    std::vector<char> ms_path(meta_header.filename_length);
    file.read(ms_path.data(), meta_header.filename_length);
    BOOST_CHECK(file.good());
    std::string ms_path_str(ms_path.begin(), ms_path.end());
    BOOST_CHECK_EQUAL(ms_path_str, "test.ms");
  };

  auto check_metadata_row = [](std::ifstream& file,
                               const BdaMetaRecord& expected) {
    BdaMetaRecord meta_record;
    meta_record.Read(file);
    BOOST_CHECK_EQUAL(meta_record.u, expected.u);
    BOOST_CHECK_EQUAL(meta_record.v, expected.v);
    BOOST_CHECK_EQUAL(meta_record.w, expected.w);
    BOOST_CHECK_EQUAL(meta_record.time, expected.time);
    BOOST_CHECK_EQUAL(meta_record.antenna1, expected.antenna1);
    BOOST_CHECK_EQUAL(meta_record.antenna2, expected.antenna2);
    BOOST_CHECK_EQUAL(meta_record.field_id, expected.field_id);
    BOOST_CHECK_EQUAL(meta_record.data_desc_id, expected.data_desc_id);
  };

  auto check_data_header = [](std::ifstream& file, size_t max_channel_count) {
    PartHeader part_header;
    part_header.Read(file);
    BOOST_CHECK(file.good());
    BOOST_CHECK_EQUAL(part_header.max_channel_count, max_channel_count);
    BOOST_CHECK(part_header.has_model);
  };

  auto check_data_contents =
      [](std::ifstream& data_file, std::ifstream& model_file,
         std::initializer_list<std::complex<float>> expected_values) {
        const size_t n_channels = expected_values.size();
        std::vector<std::complex<float>> data_buffer(n_channels);
        data_file.read(reinterpret_cast<char*>(data_buffer.data()),
                       n_channels * sizeof(std::complex<float>));
        BOOST_CHECK(data_file.good());
        BOOST_CHECK_EQUAL_COLLECTIONS(data_buffer.begin(), data_buffer.end(),
                                      expected_values.begin(),
                                      expected_values.end());
        std::vector<std::complex<float>> model_buffer(n_channels);

        model_file.read(reinterpret_cast<char*>(model_buffer.data()),
                        n_channels * sizeof(std::complex<float>));
        BOOST_CHECK(model_file.good());
        BOOST_CHECK_EQUAL_COLLECTIONS(model_buffer.begin(), model_buffer.end(),
                                      expected_values.begin(),
                                      expected_values.end());
      };

  auto check_weights_contents =
      [](std::ifstream& file, std::initializer_list<float> expected_values) {
        const size_t n_channels = expected_values.size();
        std::vector<float> weight_buffer(n_channels);
        file.read(reinterpret_cast<char*>(weight_buffer.data()),
                  n_channels * sizeof(float));
        BOOST_CHECK(file.good());

        std::vector<float>::const_iterator weight_iter = weight_buffer.begin();
        for (double expected_value : expected_values) {
          // Weights are pre-multiplied by a factor 4
          BOOST_CHECK_CLOSE_FRACTION(*weight_iter, expected_value * 4.0, 1e-6);
          ++weight_iter;
        }
      };

  for (size_t meta_file_index = 0; meta_file_index != 2; ++meta_file_index) {
    std::ifstream meta_file_0("tmp/test.ms-meta-000" +
                              std::to_string(meta_file_index) + ".tmp");
    check_metadata_header(meta_file_0, kStartTime, 4);
    check_metadata_row(meta_file_0, BdaMetaRecord{10.0, 11.0, 12.0, 1.1, 3, 4,
                                                  kFieldId, kDataDesc2});
    check_metadata_row(meta_file_0, BdaMetaRecord{13.0, 14.0, 15.0, 1.2, 4, 5,
                                                  kFieldId, kDataDesc1});
    check_metadata_row(meta_file_0, BdaMetaRecord{16.0, 17.0, 18.0, 1.3, 6, 7,
                                                  kFieldId, kDataDesc1});
    check_metadata_row(meta_file_0, BdaMetaRecord{-19.0, -20.0, -21.0, 1.4, 8,
                                                  9, kFieldId, kDataDesc2});
  }
  // Test the data contents
  std::ifstream data0_file("tmp/test.ms-part0000-I.tmp");
  std::ifstream model0_file("tmp/test.ms-part0000-I-m.tmp");
  check_data_header(data0_file, 3);
  check_data_contents(data0_file, model0_file, {101.0f, 102.0f});
  check_data_contents(data0_file, model0_file, {201.0f, 202.0f, 203.0f});
  check_data_contents(data0_file, model0_file, {301.0f, 302.0f, 303.0f});
  check_data_contents(data0_file, model0_file, {401.0f, 402.0f});

  std::ifstream weight0_file("tmp/test.ms-part0000-I-w.tmp");
  check_weights_contents(weight0_file, {0.0, 0.2});
  check_weights_contents(weight0_file, {0.1, 0.0, 0.3});
  check_weights_contents(weight0_file, {0.1, 0.0, 0.3});
  check_weights_contents(weight0_file, {0.0, 0.2});

  std::ifstream data1_file("tmp/test.ms-part0001-I.tmp");
  std::ifstream model1_file("tmp/test.ms-part0001-I-m.tmp");
  check_data_header(data1_file, 3);
  check_data_contents(data1_file, model1_file, {103.0f});
  check_data_contents(data1_file, model1_file, {204.0f, 205.0f, 206.0f});
  check_data_contents(data1_file, model1_file, {304.0f, 305.0f, 306.0f});
  check_data_contents(data1_file, model1_file, {403.0f});

  std::ifstream weight1_file("tmp/test.ms-part0001-I-w.tmp");
  check_weights_contents(weight1_file, {0.3});
  check_weights_contents(weight1_file, {0.0, 0.5, 0.6});
  check_weights_contents(weight1_file, {0.0, 0.5, 0.6});
  check_weights_contents(weight1_file, {0.3});
}

BOOST_FIXTURE_TEST_CASE(zeros_model_creation, FixtureDirectory) {
  {
    boost::filesystem::create_directory(kTemporaryDirectory);
    HandleData handle_no_model_column(
        "test.ms", "DATA", "MODEL_DATA", StorageManagerType::Default,
        kTemporaryDirectory, MakeRegularChannelMap(kChannelRanges), false,
        false, kPolsOut, kSelection, std::vector{kBands}, 6, true,
        [](HandleData) {});
    handle_no_model_column.metadata_indices_ =
        MakeMetaFilesMap(handle_no_model_column.channels_).first;
    FileWriter reordering_writer(handle_no_model_column,
                                 kMsPolarizationsPerDataDescId,
                                 kDataDescIdPerPart, kStartTime);

    reordering_writer.WriteMetaRow(0.1, 0.2, 0.3, 1.1, 0, 0, 1, 1);

    reordering_writer.WriteDataRow(kTestData.data(), kTestData.data(),
                                   kTestWeights.data(), kTestFlags, 0);

    reordering_writer.UpdateMetaHeaders();
    reordering_writer.UpdatePartHeaders(true);
    reordering_writer.PopulateModel(true, [](size_t, size_t) {});
  }

  constexpr size_t kMaxChannels = 2;
  std::ifstream model_file("tmp/test.ms-part0000-I-b0-m.tmp");
  std::vector<std::complex<float>> model_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_model_buffer(kMaxChannels,
                                                         {0.0f, 0.0f});

  model_file.read(reinterpret_cast<char*>(model_buffer.data()),
                  2 * sizeof(std::complex<float>));
  BOOST_CHECK_EQUAL_COLLECTIONS(model_buffer.begin(), model_buffer.end(),
                                expected_model_buffer.begin(),
                                expected_model_buffer.end());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace schaapcommon::reordering
