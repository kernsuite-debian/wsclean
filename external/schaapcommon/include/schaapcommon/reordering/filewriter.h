// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_REORDERED_FILE_WRITER_
#define SCHAAPCOMMON_REORDERED_FILE_WRITER_

#include "handledata.h"
#include "reordering.h"

#include <fstream>
#include <cstddef>
#include <vector>
#include <set>
#include <string>
#include <span>
#include <memory>

#include <aocommon/polarization.h>

namespace schaapcommon::reordering {

/**
 * Create an uninitialised file of exactly @p size bytes.
 * Most modern filesystems will treat the space in this file as zero
 * initialised, but use some form of sparse representation rather than incuring
 * time to fill it. On some filesystems the space may actually be completely
 * uninitialised.
 */
void AllocateFile(const std::string& filename, size_t size,
                  size_t permissions = 0664);

/**
 * Class that writes the temporary files for WSClean, containing the "reordered"
 * visibilities. Each output channel and possibly polarization is stored in
 * separate data files. The model and weight values are also stored in separate
 * files. A collection of files for the same output channel is called a 'part'
 * in the code. Each part has an associated meta data file, which may be shared
 * over multiple parts when possible. A part is that what's necessary to grid
 * one polarization of one output channel.
 *
 * Each data file starts with a header (see @ref PartHeader). It is followed by
 * the visibilities, which are stored as one flattened array, where the
 * visibilities of a row are immediately followed by the visibilities of the
 * next row. The model and weight files do not start with a header.
 *
 * Different rows may have different number of visibilities, particularly when
 * baseline-dependent averaging is used. The number of visibilities in a row can
 * be determined from the data_desc_id value in the meta data.
 *
 * In the case that so-called "instrumental" polarizations are used (e.g.
 * LOFAR's XX, XY, YX, YY correlations), for which diagonal or full-Jones
 * corrections need to be applied to correct for beam or gain solutions, a
 * single data, model data or weight file will contain multiple polarizations.
 * In that case, polarization is the quickest varying dimension, followed by
 * frequency (which together form a 'row), then followed by baseline, time, etc.
 *
 * The meta file also starts with a header (see @ref MetaHeader), followed by a
 * variable-sized string containing the measurement set path. After that,
 * the meta file contains a record for every row in the data file. For regular
 * data (not BDA), the @ref MetaRecord is used, otherwise @ref BdaMetaRecord is
 * used (which is the same, except for adding a data_desc_id field).
 */
class FileWriter {
 private:
  struct ReorderedDataFiles {
    std::unique_ptr<std::ofstream> data;
    std::unique_ptr<std::ofstream> weight;
    std::unique_ptr<std::ofstream> model;
  };

 public:
  /**
   * @param data_desc_id_per_metafile has an optional value per metadata file,
   * that is set to the data_desc_id if the metadata file has only one
   * data_desc_id. If unset, the metadata file covers multiple data_desc_ids
   * (e.g. BDA data).
   */
  FileWriter(const HandleData& data,
             const std::map<size_t, std::set<aocommon::PolarizationEnum>>&
                 ms_polarizations_per_data_desc_id,
             const std::vector<aocommon::OptionalNumber<size_t>>&
                 data_desc_id_per_metafile,
             double start_time);

  void WriteMetaRow(double u, double v, double w, double time,
                    uint32_t data_desc_id, uint32_t antenna1, uint32_t antenna2,
                    uint32_t field_id);

  void WriteDataRow(const std::complex<float>* data_array,
                    const std::complex<float>* model_array,
                    const float* weight_spectrum_array, const bool* flag_array,
                    size_t data_desc_id);
  void UpdateMetaHeaders();
  void UpdatePartHeaders(bool include_model);
  /**
   * Reserve space for model file.
   * When @p zero_initialize is true, fill the model with zeros.
   * This is needed for add-assigning to the model.
   * When @p zero_initialize is false, reserve the space but leave
   * uninitialised. This is faster if add-assigning is not necessarry.
   */
  void PopulateModel(
      bool zero_initialize,
      std::function<void(size_t progress, size_t total)> update_progress);

  ~FileWriter() = default;

  size_t MetaFileCount() const { return meta_files_.size(); }

 private:
  struct MetaFileData {
    std::string filename_;
    aocommon::OptionalNumber<size_t> data_desc_id;
    size_t selected_row_count = 0;
    std::fstream file;
  };

  aocommon::VectorMap<std::vector<MetaFileData*>> MakeDataDescIdToMetaFileMap(
      const std::map<size_t, size_t>& meta_file_indices);

  const HandleData& data_;
  std::map<size_t, std::set<aocommon::PolarizationEnum>>
      ms_polarizations_per_data_desc_id_;
  double start_time_;

  // Ordered as files[pol x channelpart]
  std::vector<ReorderedDataFiles> files_;

  std::vector<MetaFileData> meta_files_;
  // Maps data desc id to the list of meta data files that need to be written.
  aocommon::VectorMap<std::vector<MetaFileData*>>
      data_desc_id_to_meta_file_indices_;
  size_t max_channels_;
  size_t selected_rows_total_ = 0;
};

}  // namespace schaapcommon::reordering

#endif
