#ifndef AOCOMMON_MULTIBANDDATA_H_
#define AOCOMMON_MULTIBANDDATA_H_

#include <aocommon/banddata.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <algorithm>
#include <stdexcept>

namespace aocommon {
/**
 * Contains information about a set of bands. This follows the CASA Measurement
 * Set model; one MultiBandData instance can contain the band information
 * contained in the CASA Measurement Set.
 */
class MultiBandData {
 public:
  using iterator = std::vector<BandData>::iterator;
  using const_iterator = std::vector<BandData>::const_iterator;

  /**
   * Construct an empty MultiBandData.
   */
  MultiBandData() = default;

  /**
   * Construct a MultiBandData from a Measurement Set.
   * @param ms A measurement set. MultiBandData reads the spectral window table
   * and the data description table of this measurement set.
   */
  explicit MultiBandData(const casacore::MeasurementSet& ms)
      : MultiBandData(ms.spectralWindow(), ms.dataDescription()) {}

  /**
   * Construct a MultiBandData from the Measurement Set tables.
   * @param spw_table The spectral window table of a measurement set.
   * @param data_desc_table The data description table of a measurement set.
   */
  MultiBandData(const casacore::MSSpectralWindow& spw_table,
                const casacore::MSDataDescription& data_desc_table)
      : data_desc_to_band_(data_desc_table.nrow()),
        band_data_(spw_table.nrow()) {
    for (size_t spw = 0; spw != band_data_.size(); ++spw) {
      band_data_[spw] = BandData(spw_table, spw);
    }

    casacore::ScalarColumn<int> spw_column(
        data_desc_table,
        casacore::MSDataDescription::columnName(
            casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));
    for (size_t id = 0; id != data_desc_to_band_.size(); ++id)
      data_desc_to_band_[id] = spw_column(id);
  }

  /**
   * Construct a MultiBandData from another instance but only select a part of
   * each band data. This function also works when not all bands have the
   * same number of channels. If end_channel is larger than the number of
   * channels for one of the bands, the band is selected up to its last channel.
   * @param source Other instance that will be partially copied.
   * @param start_channel Start of channel range to initialize this instance
   * with.
   * @param end_channel End of channel range (exclusive) to initialize this
   * instance with.
   */
  MultiBandData(const MultiBandData& source, size_t start_channel,
                size_t end_channel)
      : data_desc_to_band_(source.data_desc_to_band_),
        band_data_(source.BandCount()) {
    for (size_t spw = 0; spw != source.BandCount(); ++spw) {
      // In case end_channel is beyond the nr of channels in this band,
      // set end_channel to the last channel of this band.
      const size_t band_end_channel =
          std::min(source.band_data_[spw].ChannelCount(), end_channel);
      if (start_channel > band_end_channel)
        throw std::runtime_error(
            "Invalid band selection: MultiBandData constructed with "
            "start_channel=" +
            std::to_string(start_channel) + ", nr of channels is " +
            std::to_string(band_end_channel) + ", source bandwidth = " +
            std::to_string(source.LowestFrequency() / 1e6) + " - " +
            std::to_string(source.HighestFrequency() / 1e6) + " MHz.");
      band_data_[spw] =
          BandData(source.band_data_[spw], start_channel, band_end_channel);
    }
  }

  /**
   * Index operator to retrieve a band data given a data_desc_id.
   * @param data_desc_id A valid data description ID for which the band is
   * returned.
   * @returns The BandData for the requested band.
   */
  const BandData& operator[](size_t data_desc_id) const {
    return band_data_[data_desc_to_band_[data_desc_id]];
  }

  /**
   * Get number of bands stored.
   * @returns Number of bands.
   */
  size_t BandCount() const { return band_data_.size(); }

  /**
   * Returns the unique number of data description IDs.
   * @returns Unique number of data desc IDs.
   */
  size_t DataDescCount() const { return data_desc_to_band_.size(); }

  /**
   * Get lowest frequency.
   * @returns The channel frequency of the channel with lowest frequency.
   */
  double LowestFrequency() const {
    if (band_data_.empty()) return 0.0;
    double freq = band_data_[0].LowestFrequency();
    for (size_t i = 0; i != band_data_.size(); ++i)
      freq = std::min(freq, band_data_[i].LowestFrequency());
    return freq;
  }

  /**
   * Get centre frequency.
   * @returns (BandStart() + BandEnd()) * 0.5.
   */
  double CentreFrequency() const { return (BandStart() + BandEnd()) * 0.5; }

  /**
   * Get highest frequency.
   * @returns The channel frequency of the channel with highest frequency.
   */
  double HighestFrequency() const {
    if (band_data_.empty()) return 0.0;
    double freq = band_data_[0].HighestFrequency();
    for (size_t i = 0; i != band_data_.size(); ++i)
      freq = std::max(freq, band_data_[i].HighestFrequency());
    return freq;
  }

  /**
   * Get total bandwidth covered.
   * @returns BandEnd() - BandStart().
   */
  double Bandwidth() const { return BandEnd() - BandStart(); }

  /**
   * Get the start frequency of the lowest frequency channel.
   * @return Start of covered bandwidth.
   */
  double BandStart() const {
    if (band_data_.empty()) return 0.0;
    double freq = std::min(band_data_[0].BandStart(), band_data_[0].BandEnd());
    for (size_t i = 0; i != band_data_.size(); ++i)
      freq = std::min(
          freq, std::min(band_data_[i].BandStart(), band_data_[i].BandEnd()));
    return freq;
  }

  /**
   * Get the end frequency of the highest frequency channel.
   * @return End of covered bandwidth.
   */
  double BandEnd() const {
    if (band_data_.empty()) return 0.0;
    double freq = std::max(band_data_[0].BandStart(), band_data_[0].BandEnd());
    for (size_t i = 0; i != band_data_.size(); ++i)
      freq = std::max(
          freq, std::max(band_data_[i].BandStart(), band_data_[i].BandEnd()));
    return freq;
  }

  /**
   * Map a data_desc_id to the corresponding band index.
   * @param data_desc_id A data_desc_id as e.g. used in a main table.
   * @returns The band index, which is equal to the row index in the spw
   * table that describes the band in a measurement set.
   */
  size_t GetBandIndex(size_t data_desc_id) const {
    return data_desc_to_band_[data_desc_id];
  }

  /**
   * Compose a list of dataDescIds that are used in the measurement set.
   * "Used" here means it is references from the main table.
   * @param main_table the measurement set.
   * @returns Set of used dataDescIds.
   */
  std::set<size_t> GetUsedDataDescIds(
      casacore::MeasurementSet& main_table) const {
    // If there is only one band, we assume it is used so as to avoid
    // scanning through the measurement set
    std::set<size_t> used_data_desc_ids;
    if (band_data_.size() == 1)
      used_data_desc_ids.insert(0);
    else {
      casacore::ScalarColumn<int> dataDescIdCol(
          main_table, casacore::MeasurementSet::columnName(
                          casacore::MSMainEnums::DATA_DESC_ID));
      for (size_t row = 0; row != main_table.nrow(); ++row) {
        size_t data_desc_id = dataDescIdCol(row);
        if (used_data_desc_ids.find(data_desc_id) == used_data_desc_ids.end())
          used_data_desc_ids.insert(data_desc_id);
      }
    }
    return used_data_desc_ids;
  }

  /**
   * Adds a new band at the end of the list of bands.
   * The band will be linked to the first available data_desc_id, which
   * is the number returned by @ref DataDescCount().
   * @returns the data_desc_id of this band.
   */
  size_t AddBand(const BandData& data) {
    const size_t data_desc_id = data_desc_to_band_.size();
    const size_t band_id = band_data_.size();
    data_desc_to_band_.emplace_back(band_id);
    band_data_.emplace_back(data);
    return data_desc_id;
  }

  iterator begin() { return band_data_.begin(); }
  const_iterator begin() const { return band_data_.begin(); }

  iterator end() { return band_data_.end(); }
  const_iterator end() const { return band_data_.end(); }

 private:
  std::vector<size_t> data_desc_to_band_;
  std::vector<BandData> band_data_;
};

}  // namespace aocommon
#endif  // AOCOMMON_MULTIBANDDATA_H_
