// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_REORDERING_H_
#define SCHAAPCOMMON_REORDERING_H_

#include <string>
#include <map>
#include <vector>
#include <memory>

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/optionalnumber.h>
#include <aocommon/polarization.h>

namespace schaapcommon::reordering {

/**
 * @file
 * Contains several helper functions for WSClean's data reordering. See
 * @ref FileWriter for a description of the reordering format.
 */

// We will create some efficiently packed structs to fetch data with 1 read.
// This will reduce the count of file-reads that are made.
// We do not apply this on the classes/structs themselves, because this may
//  reduce the data access performance.
namespace details {
#pragma pack(push, 1)
struct MetaRecordBuffer {
  double u;
  double v;
  double w;
  double time;
  uint16_t antenna1;
  uint16_t antenna2;
  uint16_t field_id;
};
struct BdaMetaRecordBuffer {
  double u;
  double v;
  double w;
  double time;
  uint16_t antenna1;
  uint16_t antenna2;
  uint16_t field_id;
  uint16_t data_desc_id;
};
struct PartHeaderBuffer {
  uint64_t max_channel_count;
  bool has_model;
};
struct MetaHeaderBuffer {
  double start_time;
  uint64_t selected_row_count;
  uint32_t filename_length;
  uint32_t data_desc_id;
  bool has_single_data_desc_id;
};
#pragma pack(pop)
}  // namespace details

struct MetaHeader {
  double start_time = 0.0;
  uint64_t selected_row_count = 0;
  uint32_t filename_length = 0;
  /// Set only if the metadata file covers only one data desc id.
  aocommon::OptionalNumber<size_t> data_desc_id;

  void Read(std::istream& str) {
    details::MetaHeaderBuffer meta_header_buffer;
    str.read(reinterpret_cast<char*>(&meta_header_buffer),
             sizeof(details::MetaHeaderBuffer));
    start_time = meta_header_buffer.start_time;
    selected_row_count = meta_header_buffer.selected_row_count;
    filename_length = meta_header_buffer.filename_length;
    if (meta_header_buffer.has_single_data_desc_id)
      data_desc_id = meta_header_buffer.data_desc_id;
  }
  void Write(std::ostream& str) const {
    const details::MetaHeaderBuffer meta_header_buffer{
        start_time, selected_row_count, filename_length,
        static_cast<uint32_t>(data_desc_id.ValueOr(0)),
        data_desc_id.HasValue()};
    str.write(reinterpret_cast<const char*>(&meta_header_buffer),
              sizeof(details::MetaHeaderBuffer));
  }
  static constexpr size_t BINARY_SIZE =
      sizeof(start_time) + sizeof(selected_row_count) +
      sizeof(filename_length) + sizeof(uint32_t) + sizeof(bool);
  static_assert(BINARY_SIZE == sizeof(details::MetaHeaderBuffer));
  static_assert(BINARY_SIZE == 25);
};

struct MetaRecord {
  double u = 0.0, v = 0.0, w = 0.0, time = 0.0;
  uint16_t antenna1 = 0, antenna2 = 0, field_id = 0;
  static constexpr size_t BINARY_SIZE =
      sizeof(double) * 4 + sizeof(uint16_t) * 3;
  static_assert(BINARY_SIZE == 38);
  void Read(std::istream& str) {
    details::MetaRecordBuffer meta_record_buffer;
    str.read(reinterpret_cast<char*>(&meta_record_buffer),
             sizeof(details::MetaRecordBuffer));

    u = meta_record_buffer.u;
    v = meta_record_buffer.v;
    w = meta_record_buffer.w;
    time = meta_record_buffer.time;
    antenna1 = meta_record_buffer.antenna1;
    antenna2 = meta_record_buffer.antenna2;
    field_id = meta_record_buffer.field_id;
  }
  void Write(std::ostream& str) const {
    const details::MetaRecordBuffer meta_record_buffer{
        u, v, w, time, antenna1, antenna2, field_id};
    str.write(reinterpret_cast<const char*>(&meta_record_buffer),
              sizeof(details::MetaRecordBuffer));
  }
};

struct BdaMetaRecord {
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;
  double time = 0.0;
  uint16_t antenna1 = 0;
  uint16_t antenna2 = 0;
  uint16_t field_id = 0;
  uint16_t data_desc_id = 0;
  static constexpr size_t BINARY_SIZE =
      sizeof(double) * 4 + sizeof(uint16_t) * 4;
  static_assert(BINARY_SIZE == 40);
  void Read(std::istream& str) {
    details::BdaMetaRecordBuffer meta_record_buffer;
    str.read(reinterpret_cast<char*>(&meta_record_buffer),
             sizeof(details::BdaMetaRecordBuffer));

    u = meta_record_buffer.u;
    v = meta_record_buffer.v;
    w = meta_record_buffer.w;
    time = meta_record_buffer.time;
    antenna1 = meta_record_buffer.antenna1;
    antenna2 = meta_record_buffer.antenna2;
    field_id = meta_record_buffer.field_id;
    data_desc_id = meta_record_buffer.data_desc_id;
  }
  void Write(std::ostream& str) const {
    const details::BdaMetaRecordBuffer meta_record_buffer{
        u, v, w, time, antenna1, antenna2, field_id, data_desc_id};
    str.write(reinterpret_cast<const char*>(&meta_record_buffer),
              sizeof(details::BdaMetaRecordBuffer));
  }
  bool operator==(const BdaMetaRecord&) const = default;
};

struct PartHeader {
  uint64_t max_channel_count = 0;
  bool has_model = false;
  static constexpr size_t BINARY_SIZE =
      sizeof(max_channel_count) + sizeof(has_model);
  static_assert(BINARY_SIZE == 9);
  void Read(std::istream& str) {
    details::PartHeaderBuffer part_header_buffer;
    str.read(reinterpret_cast<char*>(&part_header_buffer),
             sizeof(details::PartHeaderBuffer));
    max_channel_count = part_header_buffer.max_channel_count;
    has_model = part_header_buffer.has_model;
  }
  void Write(std::ostream& str) const {
    const details::PartHeaderBuffer part_header_buffer{max_channel_count,
                                                       has_model};
    str.write(reinterpret_cast<const char*>(&part_header_buffer),
              sizeof(details::PartHeaderBuffer));
  }
};

std::string GetFilenamePrefix(const std::string& ms_path,
                              const std::string& temp_dir);
std::string GetPartPrefix(const std::string& ms_path, size_t part_index,
                          aocommon::PolarizationEnum pol,
                          const std::string& temp_dir);
std::string GetMetaFilename(const std::string& ms_path,
                            const std::string& temp_dir, size_t meta_file_id);

template <typename NumType>
bool IsCFinite(const std::complex<NumType>& c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}

void ExtractData(std::complex<float>* dest, size_t start_channel,
                 size_t end_channel,
                 const std::set<aocommon::PolarizationEnum>& pols_in,
                 const std::complex<float>* data,
                 aocommon::PolarizationEnum pol_out);

template <typename NumType>
void ExtractWeights(NumType* dest, size_t start_channel, size_t end_channel,
                    const std::set<aocommon::PolarizationEnum>& pols_in,
                    const std::complex<float>* data, const float* weights,
                    const bool* flags, aocommon::PolarizationEnum pol_out);

template <bool add>
void StoreData(std::complex<float>* dest, size_t start_channel,
               size_t end_channel,
               const std::set<aocommon::PolarizationEnum>& pols_dest,
               const std::complex<float>* source,
               aocommon::PolarizationEnum pol_source);

void StoreWeights(float* dest, size_t start_channel, size_t end_channel,
                  const std::set<aocommon::PolarizationEnum>& pols_dest,
                  const float* source, aocommon::PolarizationEnum pol_source);

}  // namespace schaapcommon::reordering

#endif
