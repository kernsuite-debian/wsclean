// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "filewriter.h"

#include <fcntl.h>
#include <functional>
#include <stdexcept>

#include <aocommon/logger.h>
#include <aocommon/polarization.h>

using aocommon::Logger;

namespace schaapcommon::reordering {

/**
 * Create an uninitialised file of exactly @p size bytes.
 * Most modern filesystems will treat the space in this file as zero
 * initialised, but use some form of sparse representation rather than incuring
 * time to fill it. On some filesystems the space may actually be completely
 * uninitialised.
 */
void AllocateFile(const std::string& filename, size_t size,
                  size_t permissions) {
  const int file = open(filename.c_str(), O_RDWR | O_CREAT, permissions);
  if (file == -1) {
    throw std::runtime_error("Error creating " + filename);
  }
  if (size > 0) {
    if (lseek(file, std::max<size_t>(1, size) - 1, SEEK_SET) == -1) {
      close(file);
      remove(filename.c_str());
      throw std::runtime_error("Error reserving space for " + filename);
    }
    constexpr char kDummy = 0;
    if (write(file, &kDummy, 1) == -1) {
      close(file);
      remove(filename.c_str());
      throw std::runtime_error("Error writing " + filename);
    }
  }
  close(file);
}

FileWriter::FileWriter(
    const HandleData& data,
    const std::map<size_t, std::set<aocommon::PolarizationEnum>>&
        ms_polarizations_per_data_desc_id,
    const std::vector<aocommon::OptionalNumber<size_t>>& data_desc_id_per_file,
    double start_time)
    : data_(data),
      ms_polarizations_per_data_desc_id_(ms_polarizations_per_data_desc_id),
      start_time_(start_time),
      files_(data_.channels_.size() * data_.polarizations_.size()),
      max_channels_(GetMaxChannels(data_.channels_)) {
  Logger::Debug << "Reordering in " << data_.channels_.size() << " channels:\n";
  for (const aocommon::VectorMap<ChannelRange>& ranges : data_.channels_) {
    Logger::Debug << "  *";
    for (const ChannelRange& range : ranges) {
      if (!range.Empty()) {
        Logger::Debug << ' ' << range.data_desc_id << ':' << range.start << '-'
                      << range.end;
      }
    }
    Logger::Debug << '\n';
  }

  size_t file_index = 0;
  for (size_t part = 0; part != data_.channels_.size(); ++part) {
    for (aocommon::PolarizationEnum p : data_.polarizations_) {
      ReorderedDataFiles& file = files_[file_index];
      const std::string part_prefix =
          GetPartPrefix(data_.ms_path_, part, p, data_.temporary_directory_);
      file.data = std::make_unique<std::ofstream>(part_prefix + ".tmp");
      file.weight = std::make_unique<std::ofstream>(part_prefix + "-w.tmp");
      if (data_.initial_model_required_) {
        file.model = std::make_unique<std::ofstream>(part_prefix + "-m.tmp");
      }
      file.data->seekp(PartHeader::BINARY_SIZE, std::ios::beg);

      ++file_index;
    }
  }

  meta_files_.resize(data_desc_id_per_file.size());
  data_desc_id_to_meta_file_indices_ =
      MakeDataDescIdToMetaFileMap(data_.metadata_indices_);

  for (const std::pair<const size_t, size_t>& part_and_meta_file_index :
       data_.metadata_indices_) {
    MetaFileData& meta_file_data = meta_files_[part_and_meta_file_index.second];
    meta_file_data.filename_ =
        GetMetaFilename(data_.ms_path_, data_.temporary_directory_,
                        part_and_meta_file_index.second);
    meta_file_data.data_desc_id =
        data_desc_id_per_file[part_and_meta_file_index.second];
    meta_file_data.file = std::fstream(
        meta_file_data.filename_,
        std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!meta_file_data.file.good()) {
      throw std::runtime_error("Error opening temporary meta file " +
                               meta_file_data.filename_);
    }

    // Skip over the header ; will be written later.
    meta_file_data.file.seekp(reordering::MetaHeader::BINARY_SIZE +
                              data_.ms_path_.size());
    if (meta_file_data.file.fail()) {
      throw std::runtime_error("Error seeking through temporary meta file " +
                               meta_file_data.filename_);
    }
  }
}

aocommon::VectorMap<std::vector<FileWriter::MetaFileData*>>
FileWriter::MakeDataDescIdToMetaFileMap(
    const std::map<size_t, size_t>& meta_file_indices) {
  std::map<size_t, std::set<size_t>> data_desc_id_to_meta_indices;
  for (const std::pair<const size_t, size_t>& item : meta_file_indices) {
    const size_t part = item.first;
    const size_t meta_file_index = item.second;
    const aocommon::VectorMap<schaapcommon::reordering::ChannelRange>& ranges =
        data_.channels_[part];
    for (const schaapcommon::reordering::ChannelRange& range : ranges) {
      if (!range.Empty()) {
        data_desc_id_to_meta_indices[range.data_desc_id].insert(
            meta_file_index);
      }
    }
  }
  // Convert from temporary map<.., set> structure
  aocommon::VectorMap<std::vector<MetaFileData*>> result;
  for (const std::pair<const size_t, std::set<size_t>>& data_desc_item :
       data_desc_id_to_meta_indices) {
    const size_t data_desc_id = data_desc_item.first;
    const std::set<size_t>& meta_indices = data_desc_item.second;
    std::vector<MetaFileData*>& new_list = result.AlwaysEmplace(data_desc_id);
    new_list.reserve(meta_indices.size());
    for (size_t meta_index : meta_indices) {
      new_list.emplace_back(&meta_files_[meta_index]);
    }
  }
  return result;
}

void FileWriter::WriteMetaRow(double u, double v, double w, double time,
                              uint32_t data_desc_id, uint32_t antenna1,
                              uint32_t antenna2, uint32_t field_id) {
  const std::vector<MetaFileData*>& meta_file_list =
      data_desc_id_to_meta_file_indices_[data_desc_id];
  for (MetaFileData* file_data : meta_file_list) {
    if (file_data->data_desc_id.HasValue()) {
      reordering::MetaRecord meta;
      meta.u = u;
      meta.v = v;
      meta.w = w;
      meta.time = time;
      meta.antenna1 = antenna1;
      meta.antenna2 = antenna2;
      meta.field_id = field_id;
      ++selected_rows_total_;
      meta.Write(file_data->file);
    } else {
      reordering::BdaMetaRecord meta;
      meta.u = u;
      meta.v = v;
      meta.w = w;
      meta.time = time;
      meta.antenna1 = antenna1;
      meta.antenna2 = antenna2;
      meta.field_id = field_id;
      meta.data_desc_id = data_desc_id;
      ++selected_rows_total_;
      meta.Write(file_data->file);
    }
    ++file_data->selected_row_count;
    if (!file_data->file.good()) {
      throw std::runtime_error("Error writing to temporary file");
    }
  }
}

void FileWriter::WriteDataRow(const std::complex<float>* data_array,
                              const std::complex<float>* model_array,
                              const float* weight_spectrum_array,
                              const bool* flag_array, size_t data_desc_id) {
  const size_t polarizations_per_file =
      aocommon::Polarization::GetVisibilityCount(*data_.polarizations_.begin());
  std::vector<std::complex<float>> data_buffer(polarizations_per_file *
                                               max_channels_);
  std::vector<float> weight_buffer(polarizations_per_file * max_channels_);

  size_t file_index = 0;
  for (const aocommon::VectorMap<ChannelRange>& ranges : data_.channels_) {
    if (ContainsDataDescId(ranges, data_desc_id)) {
      const size_t part_start_ch = ranges[data_desc_id].start;
      const size_t part_end_ch = ranges[data_desc_id].end;
      const std::set<aocommon::PolarizationEnum>& ms_polarizations =
          ms_polarizations_per_data_desc_id_.find(data_desc_id)->second;

      for (aocommon::PolarizationEnum p : data_.polarizations_) {
        ReorderedDataFiles& f = files_[file_index];
        reordering::ExtractData(data_buffer.data(), part_start_ch, part_end_ch,
                                ms_polarizations, data_array, p);
        f.data->write(reinterpret_cast<char*>(data_buffer.data()),
                      (part_end_ch - part_start_ch) *
                          sizeof(std::complex<float>) * polarizations_per_file);
        if (!f.data->good()) {
          throw std::runtime_error("Error writing to temporary data file");
        }

        if (data_.initial_model_required_) {
          reordering::ExtractData(data_buffer.data(), part_start_ch,
                                  part_end_ch, ms_polarizations, model_array,
                                  p);
          f.model->write(reinterpret_cast<char*>(data_buffer.data()),
                         (part_end_ch - part_start_ch) *
                             sizeof(std::complex<float>) *
                             polarizations_per_file);
          if (!f.model->good()) {
            throw std::runtime_error(
                "Error writing to temporary model data file");
          }
        }

        reordering::ExtractWeights(weight_buffer.data(), part_start_ch,
                                   part_end_ch, ms_polarizations, data_array,
                                   weight_spectrum_array, flag_array, p);
        f.weight->write(reinterpret_cast<char*>(weight_buffer.data()),
                        (part_end_ch - part_start_ch) * sizeof(float) *
                            polarizations_per_file);
        if (!f.weight->good()) {
          throw std::runtime_error("Error writing to temporary weights file");
        }
        ++file_index;
      }
    } else {
      file_index += data_.polarizations_.size();
    }
  }
}

void FileWriter::UpdateMetaHeaders() {
  // Rewrite meta header to include selected row count
  for (MetaFileData& file_data : meta_files_) {
    // Data narrowed for reordered file
    reordering::MetaHeader meta_header;
    meta_header.start_time = start_time_;
    meta_header.selected_row_count = file_data.selected_row_count;
    meta_header.filename_length = data_.ms_path_.size();
    meta_header.data_desc_id = file_data.data_desc_id;
    file_data.file.seekp(0, std::ios::beg);
    if (!file_data.file.good()) {
      throw std::runtime_error("Error seeking through temporary meta file " +
                               file_data.filename_);
    }
    meta_header.Write(file_data.file);
    file_data.file.write(data_.ms_path_.c_str(), data_.ms_path_.size());
    if (!file_data.file.good()) {
      throw std::runtime_error("Error writing header to temporary meta file" +
                               file_data.filename_);
    }
  }
}

void FileWriter::UpdatePartHeaders(bool include_model) {
  size_t file_index = 0;
  for (const aocommon::VectorMap<ChannelRange>& ranges : data_.channels_) {
    reordering::PartHeader header;
    header.max_channel_count = GetMaxChannels(ranges);
    header.has_model = include_model;

    for ([[maybe_unused]] const aocommon::PolarizationEnum& pol :
         data_.polarizations_) {
      ReorderedDataFiles& file = files_[file_index];
      file.data->seekp(0, std::ios::beg);
      header.Write(*file.data);
      if (!file.data->good()) {
        throw std::runtime_error("Error writing to temporary data file");
      }
      ++file_index;
    }
  }
}

void FileWriter::PopulateModel(
    bool zero_initialize,
    std::function<void(size_t progress, size_t total)> update_progress) {
  const size_t polarizations_per_row =
      aocommon::Polarization::GetVisibilityCount(*data_.polarizations_.begin());
  std::vector<std::complex<float>> data_buffer(
      polarizations_per_row * max_channels_, {0.0, 0.0});

  for (size_t part = 0; part != data_.channels_.size(); ++part) {
    const size_t meta_file_index = data_.metadata_indices_.find(part)->second;
    MetaFileData& meta_file_data = meta_files_[meta_file_index];
    const size_t selected_row_count = meta_file_data.selected_row_count;
    const aocommon::MultiBandData& bands = data_.bands_per_part_[part];
    size_t n_visibilities = 0;
    if (meta_file_data.data_desc_id.HasValue()) {
      n_visibilities = selected_row_count * max_channels_;
    } else {
      MetaHeader header;
      meta_file_data.file.seekg(0, std::ios::beg);
      header.Read(meta_file_data.file);
      meta_file_data.file.seekg(header.filename_length, std::ios::cur);
      for (size_t row = 0; row != selected_row_count; ++row) {
        BdaMetaRecord record;
        record.Read(meta_file_data.file);
        n_visibilities += bands[record.data_desc_id].ChannelCount();
      }
    }

    for (const aocommon::PolarizationEnum& pol : data_.polarizations_) {
      const std::string part_prefix = reordering::GetPartPrefix(
          data_.ms_path_, part, pol, data_.temporary_directory_);
      const std::string model_filename = part_prefix + "-m.tmp";

      AllocateFile(model_filename, n_visibilities * polarizations_per_row *
                                       sizeof(std::complex<float>));
    }
    update_progress(part, data_.channels_.size() + 1);
  }
}

}  // namespace schaapcommon::reordering
