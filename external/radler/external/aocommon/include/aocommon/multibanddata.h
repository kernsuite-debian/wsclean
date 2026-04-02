#ifndef AOCOMMON_MULTIBANDDATA_H_
#define AOCOMMON_MULTIBANDDATA_H_

#include <aocommon/banddata.h>

// Because current Casacore versions aren't suporting C++20 yet in the Ubuntu
// versions, the casacore overloads are removed using this macro. This prevents
// having to compile Casacore in CI only for a few overloads that are not tested
// anyway. This check can be removed once Casacore support is improved.
#ifndef DISABLE_CASACORE_IN_BANDDATA
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#endif

#include <algorithm>
#include <cassert>
#include <set>
#include <stdexcept>

#include "optionalnumber.h"
#include "vectormap.h"

#include "io/serialistream.h"
#include "io/serialostream.h"

namespace aocommon {

/**
 * Iterator over a multibanddata that provides const BandData references.
 */
class MultiBandDataConstIterator {
 public:
  using value_type = const BandData;
  using pointer = const BandData*;
  using reference = const BandData&;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;

  MultiBandDataConstIterator& operator++() {
    assert(iterator_ != end_);
    ++iterator_;
    while (iterator_ != end_ && !iterator_->first) ++iterator_;
    return *this;
  }

  MultiBandDataConstIterator operator++(int) {
    assert(iterator_ != end_);
    MultiBandDataConstIterator before(*this);
    ++iterator_;
    while (iterator_ != end_ && !iterator_->first) ++iterator_;
    return before;
  }

  const BandData& operator*() const { return iterator_->second; }
  const BandData* operator->() const { return &iterator_->second; }

  bool operator==(const MultiBandDataConstIterator& rhs) const {
    return iterator_ == rhs.iterator_;
  }

 protected:
  using parent_const_iterator = VectorMap<
      std::pair<aocommon::OptionalNumber<size_t>, BandData>>::const_iterator;
  MultiBandDataConstIterator(parent_const_iterator iter,
                             parent_const_iterator end_iter)
      : iterator_(iter), end_(end_iter) {
    while (iterator_ != end_ && !iterator_->first) ++iterator_;
  }

  parent_const_iterator iterator_;
  parent_const_iterator end_;
  friend class MultiBandData;
};

/**
 * Iterator over a multibanddata that provides (writable) BandData references.
 */
class MultiBandDataIterator final : public MultiBandDataConstIterator {
 public:
  using value_type = BandData;
  using pointer = BandData*;
  using reference = BandData&;

  BandData& operator*() const {
    // Because MultiBandDataConstIterator is used as base, iterator_ is a
    // const iterator. It's however certain that it points to a non-const
    // object, because the constructor of MultiBandDataIterator takes
    // a non-const iterator.
    return const_cast<BandData&>(iterator_->second);
  }
  BandData* operator->() const {
    return const_cast<BandData*>(&iterator_->second);
  }

 private:
  using parent_iterator = VectorMap<
      std::pair<aocommon::OptionalNumber<size_t>, BandData>>::iterator;
  MultiBandDataIterator(parent_iterator iter, parent_iterator end_iter)
      : MultiBandDataConstIterator(iter, end_iter) {}
  friend class MultiBandData;
};

/**
 * Iterator over a multibanddata that provides (only) the data desc ids.
 */
class DataDescIdIterator {
 public:
  using value_type = size_t;
  using pointer = size_t*;
  // 'reference' should be the type returned by operator* (not necessarily a
  // reference for an input iterator).
  using reference = size_t;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;

  DataDescIdIterator& operator++() {
    assert(iterator_ != data_->end());
    ++iterator_;
    while (iterator_ != data_->end() && !iterator_->first) ++iterator_;
    return *this;
  }

  DataDescIdIterator operator++(int) {
    assert(iterator_ != data_->end());
    DataDescIdIterator before(*this);
    ++iterator_;
    while (iterator_ != data_->end() && !iterator_->first) ++iterator_;
    return before;
  }

  size_t operator*() const {
    assert(iterator_ != data_->end());
    assert(iterator_->first);
    return iterator_ - data_->begin();
  }

  bool operator==(const DataDescIdIterator& rhs) const {
    return iterator_ == rhs.iterator_;
  }

 private:
  friend class DataDescIdRange;
  using element_type = std::pair<aocommon::OptionalNumber<size_t>, BandData>;
  using parent_const_iterator = VectorMap<element_type>::const_iterator;
  DataDescIdIterator(parent_const_iterator iter,
                     const VectorMap<element_type>& data)
      : iterator_(iter), data_(&data) {
    while (iterator_ != data_->end() && !iterator_->first) ++iterator_;
  }

  parent_const_iterator iterator_;
  const VectorMap<element_type>* data_;
};

/**
 * Class that can provide the begin and end iterators for iterating over the
 * data desc ids of a multiband. This is used in @ref
 * MultiBandData::DataDescIds().
 */
class DataDescIdRange {
 public:
  using iterator = DataDescIdIterator;
  using const_iterator = DataDescIdIterator;

  DataDescIdIterator begin() const {
    return DataDescIdIterator(data_->begin(), *data_);
  }

  DataDescIdIterator end() const {
    return DataDescIdIterator(data_->end(), *data_);
  }

 private:
  friend class MultiBandData;
  using element_type = DataDescIdIterator::element_type;
  DataDescIdRange(const VectorMap<element_type>& data) : data_(&data) {}

  const VectorMap<element_type>* data_;
};

/**
 * Contains information about a set of bands. This follows the CASA Measurement
 * Set model; one MultiBandData instance can contain the band information
 * contained in the CASA Measurement Set.
 *
 * The interface allows "missing data descriptions", e.g. the class supports
 * that only data desc id 0 and 3 are defined. Before accessing a band by its
 * data desc id, the caller should check using @ref HasDataDescId() if the data
 * desc id exists.
 *
 * The casacore measurement set does not allow such 'missing' data desc ids,
 * because data descriptions are indexed by its row numbers (and an 'empty' row
 * is not allowed). The X-radio format may however specify non-consecutive
 * band/data description ids. Missing data desc ids could also occur because
 * some data desc ids are not selected in the wsclean call.
 */
class MultiBandData {
 public:
  using iterator = MultiBandDataIterator;
  using const_iterator = MultiBandDataConstIterator;

  /**
   * Construct an empty MultiBandData.
   */
  MultiBandData() = default;

#ifndef DISABLE_CASACORE_IN_BANDDATA
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
      : band_data_(data_desc_table.nrow()) {
    casacore::ScalarColumn<int> spw_column(
        data_desc_table,
        casacore::MSDataDescription::columnName(
            casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));
    for (size_t id = 0; id != band_data_.Size(); ++id) {
      const size_t spw = spw_column(id);
      band_data_[id] = std::pair(spw, BandData(spw_table, spw));
    }
  }
#endif

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
      : band_data_(source.band_data_.Size()) {
    for (size_t data_desc_id : source.DataDescIds()) {
      const element_type& source_band = source.band_data_[data_desc_id];
      // In case end_channel is beyond the nr of channels in this band,
      // set end_channel to the last channel of this band.
      const size_t band_end_channel =
          std::min(source_band.second.ChannelCount(), end_channel);
      if (start_channel > band_end_channel)
        throw std::runtime_error(
            "Invalid band selection: MultiBandData constructed with "
            "start_channel=" +
            std::to_string(start_channel) + ", nr of channels is " +
            std::to_string(band_end_channel) + ", source bandwidth = " +
            std::to_string(source.LowestFrequency() / 1e6) + " - " +
            std::to_string(source.HighestFrequency() / 1e6) + " MHz.");
      band_data_[data_desc_id] = std::pair(
          source_band.first,
          BandData(source_band.second, start_channel, band_end_channel));
    }
  }

  /**
   * Index operator to retrieve a band data given a data_desc_id.
   * @param data_desc_id A valid data description ID for which the band is
   * returned.
   * @returns The BandData for the requested band.
   */
  const BandData& operator[](size_t data_desc_id) const {
    assert(HasDataDescId(data_desc_id));
    return band_data_[data_desc_id].second;
  }

  /**
   * Get number of bands stored.
   * @returns Number of bands.
   */
  size_t BandCount() const {
    return std::count_if(band_data_.begin(), band_data_.end(),
                         [](const element_type& element) -> bool {
                           return element.first.HasValue();
                         });
  }

  /**
   * Returns the largest data description id, or zero if empty.
   * @returns Unique number of data desc IDs.
   */
  size_t HighestDataDescId() const {
    return std::max<size_t>(1, band_data_.Size()) - 1;
  }

  /**
   * Returns the largest data description id, or zero if empty.
   * @returns Unique number of data desc IDs.
   */
  size_t HighestBandId() const {
    constexpr auto compare = [](const element_type& lhs,
                                const element_type& rhs) -> bool {
      return lhs.first < rhs.first;
    };
    return band_data_.Empty() ? 0
                              : *std::max_element(band_data_.begin(),
                                                  band_data_.end(), compare)
                                     ->first;
  }

  /**
   * Get lowest frequency. Returns zero when empty.
   * @returns The channel frequency of the channel with lowest frequency.
   */
  double LowestFrequency() const {
    double freq = std::numeric_limits<double>::max();
    for (size_t i = 0; i != band_data_.Size(); ++i) {
      if (band_data_[i].first)
        freq = std::min(freq, band_data_[i].second.LowestFrequency());
    }
    return freq == std::numeric_limits<double>::max() ? 0.0 : freq;
  }

  /**
   * Get smallest wavelength. Returns zero when empty.
   */
  double ShortestWavelength() const {
    double wavelength = std::numeric_limits<double>::max();
    for (size_t i = 0; i != band_data_.Size(); ++i) {
      if (band_data_[i].first)
        wavelength =
            std::min(wavelength, band_data_[i].second.SmallestWavelength());
    }
    return wavelength == std::numeric_limits<double>::max() ? 0.0 : wavelength;
  }

  /**
   * Get centre frequency.
   * @returns (BandStart() + BandEnd()) * 0.5.
   */
  double CentreFrequency() const { return (BandStart() + BandEnd()) * 0.5; }

  /**
   * Get highest frequency. Returns zero when empty.
   * @returns The channel frequency of the channel with highest frequency.
   */
  double HighestFrequency() const {
    double freq = 0.0;
    for (const element_type& band : band_data_) {
      freq = std::max(freq, band.second.HighestFrequency());
    }
    return freq;
  }

  /**
   * Get longest wavelength. Returns zero when empty.
   */
  double LongestWavelength() const {
    double wavelength = 0.0;
    for (const element_type& band : band_data_) {
      wavelength = std::max(wavelength, band.second.LongestWavelength());
    }
    return wavelength;
  }

  /**
   * Get total bandwidth covered.
   * @returns BandEnd() - BandStart().
   */
  double Bandwidth() const { return BandEnd() - BandStart(); }

  /**
   * Get the start frequency of the lowest frequency channel, or zero when
   * empty.
   * @return Start of covered bandwidth.
   */
  double BandStart() const {
    double freq = std::numeric_limits<double>::max();
    for (size_t i = 0; i != band_data_.Size(); ++i) {
      if (band_data_[i].first)
        freq = std::min(freq, std::min(band_data_[i].second.BandStart(),
                                       band_data_[i].second.BandEnd()));
    }
    return freq == std::numeric_limits<double>::max() ? 0.0 : freq;
  }

  /**
   * Get the end frequency of the highest frequency channel, or zero when empty.
   * @return End of covered bandwidth.
   */
  double BandEnd() const {
    double freq = 0.0;
    for (size_t i = 0; i != band_data_.Size(); ++i)
      freq = std::max(freq, std::max(band_data_[i].second.BandStart(),
                                     band_data_[i].second.BandEnd()));
    return freq;
  }

  /**
   * Map a data_desc_id to the corresponding band index.
   * @param data_desc_id A data_desc_id as e.g. used in a main table.
   * @returns The band index, which is equal to the row index in the spw
   * table that describes the band in a measurement set.
   */
  size_t GetBandIndex(size_t data_desc_id) const {
    assert(HasDataDescId(data_desc_id));
    return *band_data_[data_desc_id].first;
  }

  /**
   * Returns true if this multibanddata has a band associated with the
   * specified data desc id.
   */
  bool HasDataDescId(size_t data_desc_id) const {
    return data_desc_id < band_data_.Size() &&
           band_data_[data_desc_id].first.HasValue();
  }

  size_t MaxBandChannels() const {
    size_t max_channels = 0;
    for (const element_type& band : band_data_)
      max_channels = std::max(max_channels, band.second.ChannelCount());
    return max_channels;
  }

  /**
   * Returns the data desc id of the band with the largest number of channels.
   * If there are multiple, it will return the first one.
   */
  size_t DataDescIdWithMaxChannels() const {
    size_t max_channel_data_desc_id = 0;
    size_t max_channels = 0;
    for (size_t data_desc_id : DataDescIds()) {
      if (band_data_[data_desc_id].second.ChannelCount() > max_channels) {
        max_channel_data_desc_id = data_desc_id;
        max_channels = band_data_[data_desc_id].second.ChannelCount();
      }
    }
    return max_channel_data_desc_id;
  }

#ifndef DISABLE_CASACORE_IN_BANDDATA
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
    if (band_data_.Size() == 1)
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
#endif

  /**
   * Adds a new band at the end of the list of bands.
   * The band will be linked to the first available data_desc_id.
   * @returns the data_desc_id of this band.
   */
  size_t AddBand(const BandData& data) {
    const size_t data_desc_id = band_data_.Size();
    const size_t band_id = band_data_.Empty() ? 0 : HighestBandId() + 1;
    band_data_.EmplaceBack(band_id, data);
    return data_desc_id;
  }

  /**
   * Add or replace a band and associate it with a specified data_desc_id.
   */
  void SetBand(size_t data_desc_id, const BandData& data) {
    const size_t band_id = band_data_.Empty() ? 0 : HighestBandId() + 1;
    band_data_.AlwaysEmplace(data_desc_id, element_type(band_id, data));
  }

  iterator begin() { return iterator(band_data_.begin(), band_data_.end()); }
  const_iterator begin() const {
    return const_iterator(band_data_.begin(), band_data_.end());
  }

  iterator end() { return iterator(band_data_.end(), band_data_.end()); }
  const_iterator end() const {
    return const_iterator(band_data_.end(), band_data_.end());
  }

  /**
   * Returns a class that can be used to iterate over the DataDescIds in the
   * MultiBandData. The typical usage is inside a ranged for construct:
   *
   *     const MultiBandData bands = ...
   *     for(size_t data_desc_id : bands.DataDescIds()) {
   *       // process data_desc_id
   *     }
   */
  DataDescIdRange DataDescIds() const { return DataDescIdRange(band_data_); }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt64(band_data_.Size());
    for (const element_type& item : band_data_) {
      stream.Bool(item.first.HasValue());
      if (item.first) {
        stream.UInt64(*item.first);
        item.second.Serialize(stream);
      }
    }
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    band_data_.Clear();
    const size_t n = stream.UInt64();
    for (size_t i = 0; i != n; ++i) {
      element_type& item = band_data_.EmplaceBack();
      if (stream.Bool()) {
        item.first = stream.UInt64();
        item.second.Unserialize(stream);
      }
    }
  }

 private:
  /**
   * This map is indexed by data_desc_id. Bands that have not
   * been set are left default constructed. The first value maps
   * the data_desc_id to the band index.
   */
  using element_type = std::pair<aocommon::OptionalNumber<size_t>, BandData>;
  VectorMap<element_type> band_data_;
};

}  // namespace aocommon
#endif  // AOCOMMON_MULTIBANDDATA_H_
