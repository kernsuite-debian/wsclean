#ifndef AOCOMMON_BANDDATA_H_
#define AOCOMMON_BANDDATA_H_

#include <stdexcept>
#include <vector>

// Because current Casacore versions aren't suporting C++20 yet in the Ubuntu
// versions, the casacore overloads are removed using this macro. This prevents
// having to compile Casacore in CI only for a few overloads that are not tested
// anyway. This check can be removed once Casacore support is improved.
#ifndef DISABLE_CASACORE_IN_BANDDATA
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#endif

#include "constants.h"
#include "io/serialistream.h"
#include "io/serialostream.h"

namespace aocommon {

/** Holds the meta data of a channel. */
class ChannelInfo {
 public:
  /** Construct a channel.
   * @param frequency Channel frequency in Hz.
   * @param width Channel width in Hz.
   */
  constexpr ChannelInfo(double frequency, double width)
      : _frequency(frequency), _width(width) {}

  /** Whether the frequency of the lhs is less than that of the rhs.
   * @param rhs ChannelInfo to compare with.
   * @returns lhs.Frequency() < rhs.Frequency()
   */
  constexpr bool operator<(const ChannelInfo& rhs) const {
    return _frequency < rhs._frequency;
  }

  /** Whether the frequency of the lhs is greater than that of the rhs.
   * @param rhs ChannelInfo to compare with.
   * @returns lhs.Frequency() > rhs.Frequency()
   */
  constexpr bool operator>(const ChannelInfo& rhs) const {
    return _frequency > rhs._frequency;
  }

  /** Whether the frequencies of lhs and rhs are the same. The channel width is
   * ignored.
   * @param rhs ChannelInfo to compare with
   * @returns lhs.Frequency() == rhs.Frequency()
   */
  constexpr bool operator==(const ChannelInfo& rhs) const {
    return _frequency == rhs._frequency;
  }

  /** Frequency of channel in Hz. */
  constexpr double Frequency() const { return _frequency; }
  /** Width of channel in Hz. */
  constexpr double Width() const { return _width; }

 private:
  double _frequency, _width;
};

/**
 * Contains information about a single band ("spectral window").
 * A band consists of a sequence of contiguous channels.
 */
class BandData {
 public:
  using iterator = std::vector<double>::iterator;
  using const_iterator = std::vector<double>::const_iterator;
  using reverse_iterator = std::vector<double>::reverse_iterator;
  using const_reverse_iterator = std::vector<double>::const_reverse_iterator;

  /**
   * Construct an empty instance.
   */
  BandData()
      : _channelFrequencies(), _frequencyStep(0.0), _referenceFrequency(0.0) {}

#ifndef DISABLE_CASACORE_IN_BANDDATA
  /**
   * Construct an instance from a spectral window table. The spectral window
   * table can only have a single entry, otherwise an exception is thrown.
   * @param spwTable The CASA Measurement Set spectral window table.
   */
  explicit BandData(const casacore::MSSpectralWindow& spwTable) {
    if (spwTable.nrow() != 1)
      throw std::runtime_error("Set should have exactly one spectral window");

    initFromTable(spwTable, 0);
  }

  /**
   * Construct an instance from a specified entry of a spectral window table.
   * @param spwTable The CASA Measurement Set spectral window table.
   * @param bandIndex The entry index of the spectral window table.
   */
  BandData(const casacore::MSSpectralWindow& spwTable, size_t bandIndex) {
    initFromTable(spwTable, bandIndex);
  }
#endif

  /**
   * Construct a new instance from a part of another band.
   * @param source Instance that is partially copied.
   * @param startChannel Start of range of channels that are copied.
   * @param endChannel End of range, exclusive.
   */
  BandData(const BandData& source, size_t startChannel, size_t endChannel)
      : _frequencyStep(source._frequencyStep),
        _referenceFrequency(source._referenceFrequency) {
    if (endChannel < startChannel)
      throw std::runtime_error("Invalid band specification");
    _channelFrequencies =
        std::vector<double>(source._channelFrequencies.begin() + startChannel,
                            source._channelFrequencies.begin() + endChannel);
    if (_channelFrequencies.empty())
      throw std::runtime_error("No channels in set");
  }

  /**
   * Construct a banddata from an array with channel infos.
   */
  BandData(const std::vector<ChannelInfo>& channels, double referenceFrequency)
      : _referenceFrequency(referenceFrequency) {
    initFromArray(channels);
  }

  /** Iterator over frequencies, pointing to first channel */
  iterator begin() { return _channelFrequencies.begin(); }
  /** Iterator over frequencies, pointing past last channel */
  iterator end() { return _channelFrequencies.end(); }
  /** Constant iterator over frequencies, pointing to first channel */
  const_iterator begin() const { return _channelFrequencies.begin(); }
  /** Constant iterator over frequencies, pointing to last channel */
  const_iterator end() const { return _channelFrequencies.end(); }

  /** Reverse iterator over frequencies, pointing to last channel */
  reverse_iterator rbegin() { return _channelFrequencies.rbegin(); }

  /** Reverse iterator over frequencies, pointing past first channel */
  reverse_iterator rend() { return _channelFrequencies.rend(); }

  /** Constant reverse iterator over frequencies, pointing to last channel */
  const_reverse_iterator rbegin() const { return _channelFrequencies.rbegin(); }

  /** Constant reverse iterator over frequencies, pointing past first channel */
  const_reverse_iterator rend() const { return _channelFrequencies.rend(); }

  /**
   * Assign new frequencies to this instance. The reference frequency
   * remains unmodified.
   * @param channelCount Number of channels.
   * @param frequencies Array of @p channelCount doubles containing the channel
   * frequencies.
   */
  void Set(size_t channelCount, const double* frequencies) {
    _channelFrequencies.assign(frequencies, frequencies + channelCount);
  }

  /** Retrieve number of channels in this band.
   * @returns Number of channels.
   */
  size_t ChannelCount() const { return _channelFrequencies.size(); }

  /** Get the frequency in Hz of a specified channel.
   * @param channelIndex Zero-indexed channel index.
   */
  double ChannelFrequency(size_t channelIndex) const {
    return _channelFrequencies[channelIndex];
  }

  /** Get the channelwidth in Hz of a specified channel.
   * @param channelIndex Zero-indexed channel index.
   */
  double ChannelWidth(size_t /*channelIndex*/) const { return _frequencyStep; }

  /** Get information of a specified channel.
   * @param channelIndex Zero-indexed channel index.
   */
  ChannelInfo Channel(size_t channelIndex) const {
    return ChannelInfo(_channelFrequencies[channelIndex], _frequencyStep);
  }

  /** Get the wavelength in m of a specified channel.
   * @param channelIndex Zero-indexed channel index.
   */
  double ChannelWavelength(size_t channelIndex) const {
    return kSpeedOfLight / _channelFrequencies[channelIndex];
  }

  /**
   * Get the frequency of the last channel.
   * In case the frequencies are stored in reverse channel order, the frequency
   * of the first channel is returned.
   * @returns Highest frequency.
   */
  double HighestFrequency() const {
    return _channelFrequencies.empty()      ? 0.0
           : lastChannel() > firstChannel() ? lastChannel()
                                            : firstChannel();
  }

  /**
   * Get the frequency of the first channel.
   * In case the frequencies are stored in reverse channel order, the frequency
   * of the last channel is returned.
   * @returns Lowest frequency.
   */
  double LowestFrequency() const {
    return _channelFrequencies.empty()
               ? 0
               : (firstChannel() < lastChannel() ? firstChannel()
                                                 : lastChannel());
  }

  /** Get the centre frequency.
   * @returns 0.5 * (HighestFrequency + LowestFrequency)
   */
  double CentreFrequency() const {
    return (HighestFrequency() + LowestFrequency()) * 0.5;
  }

  /**
   * @brief Get the reference frequency in Hz as stored in the spectral window
   * table. Can be slightly different from centre frequency.
   */
  double ReferenceFrequency() const { return _referenceFrequency; }

  /** Convert a frequency to a wavelength.
   * @param frequencyHz Frequency in Hz.
   * @returns Wavelength in m.
   */
  static double FrequencyToLambda(double frequencyHz) {
    return kSpeedOfLight / frequencyHz;
  }

  /** Get the wavelength of the central channel.
   * @returns Central channel wavelength.
   */
  double CentreWavelength() const {
    return kSpeedOfLight / ((HighestFrequency() + LowestFrequency()) * 0.5);
  }

  /** Get the distance between channels in Hz.
   * @returns Distance between channels.
   */
  double FrequencyStep() const { return _frequencyStep; }

  /** Get the wavelength of the first channel.
   * @returns longest wavelength. */
  double LongestWavelength() const {
    return _channelFrequencies.empty() ? 0 : kSpeedOfLight / LowestFrequency();
  }

  /**
   * Get the wavelength of the last channel.
   * @returns smallest wavelength.
   */
  double SmallestWavelength() const {
    return _channelFrequencies.empty() ? 0 : kSpeedOfLight / HighestFrequency();
  }

  /** Get the start of the frequency range covered by this band.
   * @returns Start of the band in Hz.
   */
  double BandStart() const { return LowestFrequency() - FrequencyStep() * 0.5; }
  /** Get the end of the frequency range covered by this band.
   * @returns End of the band in Hz. */
  double BandEnd() const { return HighestFrequency() + FrequencyStep() * 0.5; }

  /** Get the total bandwidth covered by this band.
   * @returns Bandwidth in Hz. */
  double Bandwidth() const {
    return HighestFrequency() - LowestFrequency() + FrequencyStep();
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.Vector(_channelFrequencies)
        .Double(_frequencyStep)
        .Double(_referenceFrequency);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.Vector(_channelFrequencies)
        .Double(_frequencyStep)
        .Double(_referenceFrequency);
  }

 private:
#ifndef DISABLE_CASACORE_IN_BANDDATA
  void initFromTable(const casacore::MSSpectralWindow& spwTable,
                     size_t bandIndex) {
    casacore::ScalarColumn<int> numChanCol(
        spwTable, casacore::MSSpectralWindow::columnName(
                      casacore::MSSpectralWindowEnums::NUM_CHAN));
    int n_channels;
    numChanCol.get(bandIndex, n_channels);
    if (n_channels == 0) throw std::runtime_error("No channels in set");

    casacore::ArrayColumn<double> chanFreqCol(
        spwTable, casacore::MSSpectralWindow::columnName(
                      casacore::MSSpectralWindowEnums::CHAN_FREQ));
    casacore::ArrayColumn<double> chanWidthCol(
        spwTable, casacore::MSSpectralWindow::columnName(
                      casacore::MSSpectralWindowEnums::CHAN_WIDTH));
    casacore::Array<double> channelFrequencies, channelWidths;
    chanFreqCol.get(bandIndex, channelFrequencies, true);
    chanWidthCol.get(bandIndex, channelWidths, true);

    _channelFrequencies.resize(n_channels);
    size_t index = 0;
    for (casacore::Array<double>::const_iterator i = channelFrequencies.begin();
         i != channelFrequencies.end(); ++i) {
      _channelFrequencies[index] = *i;
      ++index;
    }
    _frequencyStep = 0.0;
    index = 0;
    for (casacore::Array<double>::const_iterator i = channelWidths.begin();
         i != channelWidths.end(); ++i) {
      _frequencyStep += *i;
      ++index;
    }
    _frequencyStep /= double(index);

    casacore::ScalarColumn<double> referenceFrequencyColumn(
        spwTable, casacore::MSSpectralWindow::columnName(
                      casacore::MSSpectralWindowEnums::REF_FREQUENCY));
    _referenceFrequency = referenceFrequencyColumn(bandIndex);
  }
#endif

  void initFromArray(const std::vector<ChannelInfo>& channels) {
    _channelFrequencies.resize(channels.size());
    size_t index = 0;
    _frequencyStep = 0.0;
    for (const ChannelInfo& channel : channels) {
      _channelFrequencies[index] = channel.Frequency();
      _frequencyStep += channel.Width();
      ++index;
    }
    _frequencyStep /= double(index);
  }

  double firstChannel() const { return _channelFrequencies.front(); }
  double lastChannel() const { return _channelFrequencies.back(); }

  std::vector<double> _channelFrequencies;
  double _frequencyStep;
  double _referenceFrequency;
};

}  // namespace aocommon

#endif
