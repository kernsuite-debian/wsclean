#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include "msprovider.h"

#include "../structures/msselection.h"

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <fstream>
#include <string>
#include <map>

class PartitionedMSReader;

class PartitionedMS final : public MSProvider {
  friend class PartitionedMSReader;

 public:
  class Handle;

  struct ChannelRange {
    int dataDescId;
    size_t start, end;
    bool operator<(const ChannelRange& rhs) const {
      if (dataDescId < rhs.dataDescId) return true;
      if (dataDescId > rhs.dataDescId) return false;
      if (start < rhs.start) return true;
      if (start > rhs.start) return false;
      return end < rhs.end;
    }
  };

  PartitionedMS(const Handle& handle, size_t partIndex,
                aocommon::PolarizationEnum polarization, size_t bandIndex);

  virtual ~PartitionedMS();

  PartitionedMS(const PartitionedMS&) = delete;
  PartitionedMS& operator=(const PartitionedMS&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  SynchronizedMS MS() override {
    return SynchronizedMS(_handle._data->_msPath.data());
  }

  const std::string& DataColumnName() override {
    return _handle._data->_dataColumnName;
  }

  void NextOutputRow() override;

  void ResetWritePosition() override { _currentOutputRow = 0; };

  void WriteModel(const std::complex<float>* buffer, bool addToMS) override;

  void ReopenRW() override {}

  double StartTime() override { return _metaHeader.startTime; }

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override;

  aocommon::PolarizationEnum Polarization() override { return _polarization; }

  size_t NChannels() override { return _partHeader.channelCount; }
  size_t NPolarizations() override { return _polarizationCountInFile; }
  size_t NAntennas() override { return _handle._data->_nAntennas; }

  size_t DataDescId() override { return _partHeader.dataDescId; }

  static Handle Partition(const string& msPath,
                          const std::vector<ChannelRange>& channels,
                          class MSSelection& selection,
                          const string& dataColumnName, bool includeModel,
                          bool initialModelRequired,
                          const class Settings& settings);

  class Handle {
    // PartitionedMSReader is a friend of Handle
    // in order to access the _data member.
    friend class PartitionedMSReader;

   public:
    Handle() = default;

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    friend class PartitionedMS;

   private:
    struct HandleData {
      HandleData() : _isCopy(false) {}

      HandleData(const std::string& msPath, const string& dataColumnName,
                 const std::string& temporaryDirectory,
                 const std::vector<ChannelRange>& channels,
                 bool initialModelRequired, bool modelUpdateRequired,
                 const std::set<aocommon::PolarizationEnum>& polarizations,
                 const MSSelection& selection, size_t nAntennas)
          : _msPath(msPath),
            _dataColumnName(dataColumnName),
            _temporaryDirectory(temporaryDirectory),
            _channels(channels),
            _initialModelRequired(initialModelRequired),
            _modelUpdateRequired(modelUpdateRequired),
            _polarizations(polarizations),
            _selection(selection),
            _nAntennas(nAntennas),
            _isCopy(false) {}

      ~HandleData();

      std::string _msPath, _dataColumnName, _temporaryDirectory;
      std::vector<ChannelRange> _channels;
      bool _initialModelRequired, _modelUpdateRequired;
      std::set<aocommon::PolarizationEnum> _polarizations;
      MSSelection _selection;
      size_t _nAntennas;
      bool _isCopy;

      void Serialize(aocommon::SerialOStream& stream) const;
      void Unserialize(aocommon::SerialIStream& stream);
    };
    std::shared_ptr<HandleData> _data;

    Handle(const std::string& msPath, const string& dataColumnName,
           const std::string& temporaryDirectory,
           const std::vector<ChannelRange>& channels, bool initialModelRequired,
           bool modelUpdateRequired,
           const std::set<aocommon::PolarizationEnum>& polarizations,
           const MSSelection& selection, size_t nAntennas)
        : _data(new HandleData(msPath, dataColumnName, temporaryDirectory,
                               channels, initialModelRequired,
                               modelUpdateRequired, polarizations, selection,
                               nAntennas)) {}
  };

 private:
  static void unpartition(const Handle::HandleData& handle);

  static void getDataDescIdMap(
      std::map<size_t, size_t>& dataDescIds,
      const std::vector<PartitionedMS::ChannelRange>& channels);

  const Handle _handle;
  const size_t _partIndex;
  char* _modelFileMap;
  size_t _currentOutputRow;
  std::unique_ptr<std::ofstream> _modelDataFile;
  int _fd;
  const aocommon::PolarizationEnum _polarization;
  const size_t _polarizationCountInFile;

  struct MetaHeader {
    uint64_t selectedRowCount;
    uint32_t filenameLength;
    double startTime;
  } _metaHeader;
  struct MetaRecord {
    double u, v, w, time;
    uint16_t antenna1, antenna2, fieldId;
    static constexpr size_t BINARY_SIZE = 8 * 4 + 2 * 3;
    void read(std::istream& str) {
      str.read(reinterpret_cast<char*>(&u), sizeof(double));
      str.read(reinterpret_cast<char*>(&v), sizeof(double));
      str.read(reinterpret_cast<char*>(&w), sizeof(double));
      str.read(reinterpret_cast<char*>(&time), sizeof(double));
      str.read(reinterpret_cast<char*>(&antenna1), sizeof(uint16_t));
      str.read(reinterpret_cast<char*>(&antenna2), sizeof(uint16_t));
      str.read(reinterpret_cast<char*>(&fieldId), sizeof(uint16_t));
    }
    void write(std::ostream& str) const {
      str.write(reinterpret_cast<const char*>(&u), sizeof(double));
      str.write(reinterpret_cast<const char*>(&v), sizeof(double));
      str.write(reinterpret_cast<const char*>(&w), sizeof(double));
      str.write(reinterpret_cast<const char*>(&time), sizeof(double));
      str.write(reinterpret_cast<const char*>(&antenna1), sizeof(uint16_t));
      str.write(reinterpret_cast<const char*>(&antenna2), sizeof(uint16_t));
      str.write(reinterpret_cast<const char*>(&fieldId), sizeof(uint16_t));
    }
  };
  struct PartHeader {
    uint64_t channelCount;
    uint64_t channelStart;
    uint32_t dataDescId;
    bool hasModel;
  } _partHeader;

  static std::string getFilenamePrefix(const std::string& msPath,
                                       const std::string& tempDir);
  static std::string getPartPrefix(const std::string& msPath, size_t partIndex,
                                   aocommon::PolarizationEnum pol,
                                   size_t dataDescId,
                                   const std::string& tempDir);
  static std::string getMetaFilename(const std::string& msPath,
                                     const std::string& tempDir,
                                     size_t dataDescId);
};

#endif
