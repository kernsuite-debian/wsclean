#include "partitionedms.h"
#include "msreaders/partitionedmsreader.h"

#include "averagingmsrowprovider.h"
#include "directmsrowprovider.h"
#include "msrowprovider.h"
#include "noisemsrowprovider.h"

#include "../main/progressbar.h"
#include "../main/settings.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>

#include <aocommon/logger.h>

#include <boost/filesystem/path.hpp>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

using aocommon::Logger;

/**
 * MAP_NORESERVE is unsuported AND not defined on hurd-i386, so
 * assign it to zero in this case.
 */
#ifndef MAP_NORESERVE
#define MAP_NORESERVE 0
#endif

namespace {
struct PartitionFiles {
  std::unique_ptr<std::ofstream> data;
  std::unique_ptr<std::ofstream> weight;
  std::unique_ptr<std::ofstream> model;
};

std::map<size_t, std::vector<aocommon::PolarizationEnum>>
GetMSPolarizationsPerDataDescId(
    const std::vector<PartitionedMS::ChannelRange>& ranges,
    casacore::MeasurementSet& ms) {
  std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId;
  for (const PartitionedMS::ChannelRange& range : ranges) {
    msPolarizationsPerDataDescId.emplace(
        range.dataDescId,
        PartitionedMS::GetMSPolarizations(range.dataDescId, ms));
  }
  return msPolarizationsPerDataDescId;
}

size_t GetMaxChannels(
    const std::vector<PartitionedMS::ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const PartitionedMS::ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

}  // namespace

PartitionedMS::PartitionedMS(const Handle& handle, size_t partIndex,
                             aocommon::PolarizationEnum polarization,
                             size_t dataDescId)
    : _handle(handle),
      _partIndex(partIndex),
      _dataDescId(dataDescId),
      _currentOutputRow(0),
      _polarization(polarization),
      _polarizationCountInFile(
          aocommon::Polarization::GetVisibilityCount(_polarization)) {
  std::ifstream metaFile(getMetaFilename(
      handle._data->_msPath, handle._data->_temporaryDirectory, dataDescId));
  if (!metaFile) {
    throw std::runtime_error("Error opening meta file for ms " +
                             handle._data->_msPath + ", dataDescId " +
                             std::to_string(dataDescId));
  }

  _metaHeader.Read(metaFile);
  std::vector<char> msPath(_metaHeader.filenameLength + 1, char(0));
  metaFile.read(msPath.data(), _metaHeader.filenameLength);
  Logger::Info << "Opening reordered part " << partIndex << " spw "
               << dataDescId << " for " << msPath.data() << '\n';
  std::string partPrefix =
      getPartPrefix(msPath.data(), partIndex, polarization, dataDescId,
                    handle._data->_temporaryDirectory);

  std::ifstream dataFile(partPrefix + ".tmp", std::ios::in);
  if (!dataFile.good())
    throw std::runtime_error("Error opening temporary data file '" +
                             partPrefix + ".tmp'");
  _partHeader.Read(dataFile);
  if (!dataFile.good())
    throw std::runtime_error("Error reading header from file '" + partPrefix +
                             ".tmp'");

  if (_partHeader.hasModel) {
    const size_t length =
        _partHeader.channelCount * _metaHeader.selectedRowCount *
        _polarizationCountInFile * sizeof(std::complex<float>);
    _modelFile = MappedFile(partPrefix + "-m.tmp", length);
  }
  metaFile.close();
  dataFile.close();
}

PartitionedMS::~PartitionedMS() {}

std::unique_ptr<MSReader> PartitionedMS::MakeReader() {
  std::unique_ptr<MSReader> reader(new PartitionedMSReader(this));
  return reader;
}

void PartitionedMS::NextOutputRow() { ++_currentOutputRow; }

void PartitionedMS::WriteModel(const std::complex<float>* buffer,
                               bool addToMS) {
#ifndef NDEBUG
  if (!_partHeader.hasModel)
    throw std::runtime_error("Partitioned MS initialized without model");
#endif
  size_t rowLength = _partHeader.channelCount * _polarizationCountInFile *
                     sizeof(std::complex<float>);
  std::complex<float>* modelWritePtr = reinterpret_cast<std::complex<float>*>(
      _modelFile.Data() + rowLength * _currentOutputRow);

  // In case the value was not sampled in this pass, it has been set to infinite
  // and should not overwrite the current value in the set.
  if (addToMS) {
    for (size_t i = 0; i != _partHeader.channelCount * _polarizationCountInFile;
         ++i) {
      if (std::isfinite(buffer[i].real())) modelWritePtr[i] += buffer[i];
    }
  } else {
    for (size_t i = 0; i != _partHeader.channelCount * _polarizationCountInFile;
         ++i) {
      if (std::isfinite(buffer[i].real())) modelWritePtr[i] = buffer[i];
    }
  }
}

std::string PartitionedMS::getFilenamePrefix(const std::string& msPathStr,
                                             const std::string& tempDir) {
  boost::filesystem::path prefixPath;
  if (tempDir.empty())
    prefixPath = msPathStr;
  else {
    std::string msPathCopy(msPathStr);
    while (!msPathCopy.empty() && *msPathCopy.rbegin() == '/')
      msPathCopy.resize(msPathCopy.size() - 1);
    boost::filesystem::path msPath(msPathCopy);
    prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
  }
  std::string prefix(prefixPath.string());
  while (!prefix.empty() && *prefix.rbegin() == '/')
    prefix.resize(prefix.size() - 1);
  return prefix;
}

std::string PartitionedMS::getPartPrefix(const std::string& msPathStr,
                                         size_t partIndex,
                                         aocommon::PolarizationEnum pol,
                                         size_t dataDescId,
                                         const std::string& tempDir) {
  std::string prefix = getFilenamePrefix(msPathStr, tempDir);

  std::ostringstream partPrefix;
  partPrefix << prefix << "-part";
  if (partIndex < 1000) partPrefix << '0';
  if (partIndex < 100) partPrefix << '0';
  if (partIndex < 10) partPrefix << '0';
  partPrefix << partIndex;
  partPrefix << "-";
  partPrefix << aocommon::Polarization::TypeToShortString(pol);
  partPrefix << "-b" << dataDescId;
  return partPrefix.str();
}

string PartitionedMS::getMetaFilename(const string& msPathStr,
                                      const std::string& tempDir,
                                      size_t dataDescId) {
  std::string prefix = getFilenamePrefix(msPathStr, tempDir);

  std::ostringstream s;
  s << prefix << "-spw" << dataDescId << "-parted-meta.tmp";
  return s.str();
}

/*
 * When partitioned:
 * One global file stores:
 * - Metadata:
 *   * Number of selected rows
 *   * Filename length + string
 *   * [ UVW, dataDescId ]
 * The binary parts store the following information:
 * - Number of channels
 * - Start channel in MS
 * - Total weight in part
 * - Data    (single polarization, as requested)
 * - Weights (single)
 * - Model, optionally
 */
PartitionedMS::Handle PartitionedMS::Partition(
    const string& msPath, const std::vector<ChannelRange>& channels,
    MSSelection& selection, const string& dataColumnName, bool includeModel,
    bool initialModelRequired, const Settings& settings) {
  const bool modelUpdateRequired = settings.modelUpdateRequired;
  std::set<aocommon::PolarizationEnum> polsOut;
  if (settings.gridderType == GridderType::IDG) {
    if (settings.polarizations.size() == 1) {
      if ((settings.ddPsfGridWidth > 1 || settings.ddPsfGridHeight > 1) &&
          settings.gridWithBeam) {
        polsOut.insert(aocommon::Polarization::StokesI);
      } else {
        polsOut.insert(aocommon::Polarization::DiagonalInstrumental);
      }
    } else {
      polsOut.insert(aocommon::Polarization::Instrumental);
    }
  } else if (settings.diagonalSolutions) {
    polsOut.insert(aocommon::Polarization::DiagonalInstrumental);
  } else {
    polsOut = settings.polarizations;
  }
  const size_t polarizationsPerFile =
      aocommon::Polarization::GetVisibilityCount(*polsOut.begin());
  const std::string& temporaryDirectory = settings.temporaryDirectory;

  const size_t channelParts = channels.size();

  if (channelParts != 1) {
    Logger::Debug << "Partitioning in " << channels.size() << " channels:";
    for (size_t i = 0; i != channels.size(); ++i)
      Logger::Debug << ' ' << channels[i].dataDescId << ':' << channels[i].start
                    << '-' << channels[i].end;
  }
  Logger::Debug << '\n';

  // Ordered as files[pol x channelpart]
  std::vector<PartitionFiles> files(channelParts * polsOut.size());

  const size_t max_channels = GetMaxChannels(channels);

  // Each data desc id needs a separate meta file because they can have
  // different uvws and other info.
  size_t fileIndex = 0;
  for (size_t part = 0; part != channelParts; ++part) {
    for (aocommon::PolarizationEnum p : polsOut) {
      PartitionFiles& f = files[fileIndex];
      std::string partPrefix = getPartPrefix(
          msPath, part, p, channels[part].dataDescId, temporaryDirectory);
      f.data = std::make_unique<std::ofstream>(partPrefix + ".tmp");
      f.weight = std::make_unique<std::ofstream>(partPrefix + "-w.tmp");
      if (initialModelRequired)
        f.model = std::make_unique<std::ofstream>(partPrefix + "-m.tmp");
      f.data->seekp(PartHeader::BINARY_SIZE, std::ios::beg);

      ++fileIndex;
    }
  }

  // This maps dataDescId to spw index.
  const std::map<size_t, size_t> selectedDataDescIds =
      getDataDescIdMap(channels);

  std::unique_ptr<MsRowProviderBase> rowProvider;
  if (settings.baselineDependentAveragingInWavelengths == 0.0) {
    if (settings.simulateNoise) {
      std::unique_ptr<NoiseMSRowProvider> noiseRowProvider(
          new NoiseMSRowProvider(msPath, selection, selectedDataDescIds,
                                 dataColumnName, initialModelRequired));
      if (settings.simulatedBaselineNoiseFilename.empty())
        noiseRowProvider->SetNoiseLevel(settings.simulatedNoiseStdDev);
      else
        noiseRowProvider->SetNoiseBaselineFile(
            settings.simulatedBaselineNoiseFilename);
      rowProvider = std::move(noiseRowProvider);
    } else
      rowProvider = MakeMsRowProvider(msPath, selection, selectedDataDescIds,
                                      dataColumnName, initialModelRequired);
  } else {
    if (initialModelRequired)
      throw std::runtime_error(
          "Baseline-dependent averaging is enabled together with a mode that "
          "requires the model data (e.g. -continue or -subtract-model). This "
          "is not possible.");
    rowProvider = std::make_unique<AveragingMSRowProvider>(
        settings.baselineDependentAveragingInWavelengths, msPath, selection,
        selectedDataDescIds, settings.fieldIds[0], dataColumnName,
        initialModelRequired);
  }

  const std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId =
          GetMSPolarizationsPerDataDescId(channels, rowProvider->Ms());
  const size_t nAntennas = rowProvider->Ms().antenna().nrow();
  const aocommon::MultiBandData bands(rowProvider->Ms());

  if (settings.parallelReordering == 1)
    Logger::Info << "Reordering " << msPath << " into " << channelParts << " x "
                 << polsOut.size() << " parts.\n";

  // Write header of meta file, one meta file for each data desc id
  // TODO rather than writing we can just skip and write later
  std::vector<std::unique_ptr<std::ofstream>> metaFiles(
      selectedDataDescIds.size());
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t dataDescId = p.first;
    const size_t spwIndex = p.second;
    std::string metaFilename =
        getMetaFilename(msPath, temporaryDirectory, dataDescId);
    metaFiles[spwIndex] = std::make_unique<std::ofstream>(metaFilename);
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = 0;  // not yet known
    metaHeader.filenameLength = msPath.size();
    metaHeader.startTime = rowProvider->StartTime();
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
    if (!metaFiles[spwIndex]->good())
      throw std::runtime_error("Error writing to temporary file " +
                               metaFilename);
  }

  // Write actual data
  std::vector<std::complex<float>> dataBuffer(polarizationsPerFile *
                                              max_channels);
  std::vector<float> weightBuffer(polarizationsPerFile * max_channels);

  casacore::Array<std::complex<float>> dataArray;
  casacore::Array<std::complex<float>> modelArray;
  casacore::Array<float> weightSpectrumArray;
  casacore::Array<bool> flagArray;

  std::unique_ptr<ProgressBar> progress1;
  if (settings.parallelReordering == 1)
    progress1.reset(new ProgressBar("Reordering"));

  size_t selectedRowsTotal = 0;
  aocommon::UVector<size_t> selectedRowCountPerSpwIndex(
      selectedDataDescIds.size(), 0);
  while (!rowProvider->AtEnd()) {
    if (progress1)
      progress1->SetProgress(rowProvider->CurrentProgress(),
                             rowProvider->TotalProgress());

    MetaRecord meta;

    double time;
    uint32_t dataDescId, antenna1, antenna2, fieldId;
    rowProvider->ReadData(dataArray, flagArray, weightSpectrumArray, meta.u,
                          meta.v, meta.w, dataDescId, antenna1, antenna2,
                          fieldId, time);
    meta.antenna1 = antenna1;
    meta.antenna2 = antenna2;
    meta.fieldId = fieldId;
    meta.time = time;
    const size_t spwIndex = selectedDataDescIds.find(dataDescId)->second;
    ++selectedRowCountPerSpwIndex[spwIndex];
    ++selectedRowsTotal;
    std::ofstream& metaFile = *metaFiles[spwIndex];
    meta.Write(metaFile);
    if (!metaFile.good())
      throw std::runtime_error("Error writing to temporary file");

    if (initialModelRequired) rowProvider->ReadModel(modelArray);

    fileIndex = 0;
    for (size_t part = 0; part != channelParts; ++part) {
      if (channels[part].dataDescId == int(dataDescId)) {
        const size_t partStartCh = channels[part].start;
        const size_t partEndCh = channels[part].end;
        const std::vector<aocommon::PolarizationEnum>& msPolarizations =
            msPolarizationsPerDataDescId.find(dataDescId)->second;

        for (aocommon::PolarizationEnum p : polsOut) {
          PartitionFiles& f = files[fileIndex];
          CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
                   dataArray, p);
          f.data->write(reinterpret_cast<char*>(dataBuffer.data()),
                        (partEndCh - partStartCh) *
                            sizeof(std::complex<float>) * polarizationsPerFile);
          if (!f.data->good())
            throw std::runtime_error("Error writing to temporary data file");

          if (initialModelRequired) {
            CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
                     modelArray, p);
            f.model->write(reinterpret_cast<char*>(dataBuffer.data()),
                           (partEndCh - partStartCh) *
                               sizeof(std::complex<float>) *
                               polarizationsPerFile);
            if (!f.model->good())
              throw std::runtime_error(
                  "Error writing to temporary model data file");
          }

          CopyWeights(weightBuffer.data(), partStartCh, partEndCh,
                      msPolarizations, dataArray, weightSpectrumArray,
                      flagArray, p);
          f.weight->write(
              reinterpret_cast<char*>(weightBuffer.data()),
              (partEndCh - partStartCh) * sizeof(float) * polarizationsPerFile);
          if (!f.weight->good())
            throw std::runtime_error("Error writing to temporary weights file");
          ++fileIndex;
        }
      } else {
        fileIndex += polsOut.size();
      }
    }

    rowProvider->NextRow();
  }
  progress1.reset();
  Logger::Debug << "Total selected rows: " << selectedRowsTotal << '\n';
  rowProvider->OutputStatistics();

  // Rewrite meta headers to include selected row count
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t spwIndex = p.second;
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = selectedRowCountPerSpwIndex[spwIndex];
    metaHeader.filenameLength = msPath.size();
    metaHeader.startTime = rowProvider->StartTime();
    metaFiles[spwIndex]->seekp(0);
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
  }

  // Write header to parts and write empty model files (if requested)
  PartHeader header;
  header.hasModel = includeModel;
  fileIndex = 0;
  dataBuffer.assign(max_channels * polarizationsPerFile, 0.0);
  std::unique_ptr<ProgressBar> progress2;
  if (includeModel && !initialModelRequired && settings.parallelReordering == 1)
    progress2 =
        std::make_unique<ProgressBar>("Initializing model visibilities");
  for (size_t part = 0; part != channelParts; ++part) {
    header.channelStart = channels[part].start,
    header.channelCount = channels[part].end - header.channelStart;
    header.dataDescId = channels[part].dataDescId;
    for (std::set<aocommon::PolarizationEnum>::const_iterator p =
             polsOut.begin();
         p != polsOut.end(); ++p) {
      PartitionFiles& f = files[fileIndex];
      f.data->seekp(0, std::ios::beg);
      header.Write(*f.data);
      if (!f.data->good())
        throw std::runtime_error("Error writing to temporary data file");

      f.data.reset();
      f.weight.reset();
      f.model.reset();
      ++fileIndex;

      // If model is requested, fill model file with zeros
      if (includeModel && !initialModelRequired) {
        std::string partPrefix = getPartPrefix(
            msPath, part, *p, header.dataDescId, temporaryDirectory);
        std::ofstream modelFile(partPrefix + "-m.tmp");
        const size_t selectedRowCount = selectedRowCountPerSpwIndex
            [selectedDataDescIds.find(channels[part].dataDescId)->second];
        for (size_t i = 0; i != selectedRowCount; ++i) {
          modelFile.write(reinterpret_cast<char*>(dataBuffer.data()),
                          header.channelCount * sizeof(std::complex<float>) *
                              polarizationsPerFile);
          if (progress2)
            progress2->SetProgress(part * selectedRowCount + i,
                                   channelParts * selectedRowCount);
        }
      }
    }
  }
  progress2.reset();

  return Handle(msPath, dataColumnName, temporaryDirectory, channels,
                initialModelRequired, modelUpdateRequired, polsOut, selection,
                bands, nAntennas);
}

void PartitionedMS::unpartition(
    const PartitionedMS::Handle::HandleData& handle) {
  const std::set<aocommon::PolarizationEnum> pols = handle._polarizations;

  const std::map<size_t, size_t> dataDescIds =
      getDataDescIdMap(handle._channels);

  std::vector<MetaHeader> metaHeaders(dataDescIds.size());
  for (const std::pair<const size_t, size_t>& dataDescId : dataDescIds) {
    std::ifstream metaFile(getMetaFilename(
        handle._msPath, handle._temporaryDirectory, dataDescId.first));
    MetaHeader& metaHeader = metaHeaders[dataDescId.second];
    metaHeader.Read(metaFile);
    std::vector<char> msPath(metaHeader.filenameLength + 1, char(0));
    metaFile.read(msPath.data(), metaHeader.filenameLength);
  }

  ChannelRange firstRange = handle._channels[0];
  std::ifstream firstDataFile(
      getPartPrefix(handle._msPath, 0, *pols.begin(), firstRange.dataDescId,
                    handle._temporaryDirectory) +
          ".tmp",
      std::ios::in);
  if (!firstDataFile.good())
    throw std::runtime_error("Error opening temporary data file");
  PartHeader firstPartHeader;
  firstPartHeader.Read(firstDataFile);
  if (!firstDataFile.good())
    throw std::runtime_error("Error reading from temporary data file");

  if (firstPartHeader.hasModel) {
    const size_t channelParts = handle._channels.size();

    // Open the temporary files
    std::vector<std::unique_ptr<std::ifstream>> modelFiles(channelParts *
                                                           pols.size());
    size_t fileIndex = 0;
    for (size_t part = 0; part != channelParts; ++part) {
      size_t dataDescId = handle._channels[part].dataDescId;
      for (aocommon::PolarizationEnum p : pols) {
        std::string partPrefix = getPartPrefix(
            handle._msPath, part, p, dataDescId, handle._temporaryDirectory);
        modelFiles[fileIndex] =
            std::make_unique<std::ifstream>(partPrefix + "-m.tmp");
        if (!*modelFiles[fileIndex])
          throw std::runtime_error("Error opening temporary model data file '" +
                                   partPrefix + "-m.tmp' for reading");
        ++fileIndex;
      }
    }

    casacore::MeasurementSet ms(handle._msPath, casacore::Table::Update);
    const std::map<size_t, std::vector<aocommon::PolarizationEnum>>
        msPolarizationsPerDataDescId =
            GetMSPolarizationsPerDataDescId(handle._channels, ms);
    InitializeModelColumn(ms);
    casacore::ScalarColumn<int> antenna1Column(
        ms, ms.columnName(casacore::MSMainEnums::ANTENNA1));
    casacore::ScalarColumn<int> antenna2Column(
        ms, ms.columnName(casacore::MSMainEnums::ANTENNA2));
    casacore::ScalarColumn<int> fieldIdColumn(
        ms, ms.columnName(casacore::MSMainEnums::FIELD_ID));
    casacore::ScalarColumn<double> timeColumn(
        ms, ms.columnName(casacore::MSMainEnums::TIME));
    casacore::ScalarColumn<int> dataDescIdColumn(
        ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
    casacore::ArrayColumn<casacore::Complex> dataColumn(ms,
                                                        handle._dataColumnName);
    casacore::ArrayColumn<casacore::Complex> modelColumn(
        ms, ms.columnName(casacore::MSMainEnums::MODEL_DATA));
    casacore::ArrayColumn<double> uvwColumn(
        ms, ms.columnName(casacore::MSMainEnums::UVW));

    const casacore::IPosition shape(dataColumn.shape(0));
    const size_t max_channels = GetMaxChannels(handle._channels);

    const size_t polarizationsPerFile =
        aocommon::Polarization::GetVisibilityCount(*pols.begin());
    std::vector<std::complex<float>> modelDataBuffer(max_channels *
                                                     polarizationsPerFile);
    casacore::Array<std::complex<float>> modelDataArray(shape);

    ProgressBar progress(std::string("Writing changed model back to ") +
                         handle._msPath);
    size_t startRow, endRow;
    GetRowRange(ms, handle._selection, startRow, endRow);
    size_t timestep =
        handle._selection.HasInterval() ? handle._selection.IntervalStart() : 0;
    double time = timeColumn(startRow);
    size_t selectedRowCountForDebug = 0;
    for (size_t row = startRow; row != endRow; ++row) {
      progress.SetProgress(row - startRow, endRow - startRow);
      const int a1 = antenna1Column(row);
      const int a2 = antenna2Column(row);
      const int fieldId = fieldIdColumn(row);
      const size_t dataDescId = dataDescIdColumn(row);
      casacore::Vector<double> uvw = uvwColumn(row);

      if (time != timeColumn(row)) {
        ++timestep;
        time = timeColumn(row);
      }
      if (handle._selection.IsSelected(fieldId, timestep, a1, a2, uvw)) {
        std::map<size_t, size_t>::const_iterator dataDescIdIter =
            dataDescIds.find(dataDescId);
        if (dataDescIdIter != dataDescIds.end()) {
          modelColumn.get(row, modelDataArray, true);
          size_t fileIndex = 0;
          for (size_t part = 0; part != channelParts; ++part) {
            const size_t partDataDescId = handle._channels[part].dataDescId;
            if (dataDescId == partDataDescId) {
              const size_t partStartCh = handle._channels[part].start;
              const size_t partEndCh = handle._channels[part].end;
              const std::vector<aocommon::PolarizationEnum>& msPolarizations =
                  msPolarizationsPerDataDescId.find(dataDescId)->second;
              for (aocommon::PolarizationEnum p : pols) {
                modelFiles[fileIndex]->read(
                    reinterpret_cast<char*>(modelDataBuffer.data()),
                    (partEndCh - partStartCh) * polarizationsPerFile *
                        sizeof(std::complex<float>));
                if (!modelFiles[fileIndex]->good())
                  throw std::runtime_error(
                      "Error reading from temporary model data file");
                ReverseCopyData<false>(modelDataArray, partStartCh, partEndCh,
                                       msPolarizations, modelDataBuffer.data(),
                                       p);

                ++fileIndex;
              }
            } else {
              fileIndex += pols.size();
            }
          }
          modelColumn.put(row, modelDataArray);
          selectedRowCountForDebug++;
        }
      }
    }
    progress.SetProgress(endRow - startRow, endRow - startRow);

    Logger::Debug << "Row count during unpartitioning: "
                  << selectedRowCountForDebug << '\n';
  }
}

PartitionedMS::Handle::HandleData::~HandleData() {
  if (!_isCopy) {
    // We can't throw inside destructor, so catch potential exceptions that
    // occur during writing the measurement sets.
    try {
      // Skip writing back data if we are in the middle of handling an exception
      // (stack unwinding)
      if (std::uncaught_exceptions()) {
        Logger::Info << "An exception occurred, writing back will be skipped. "
                        "Cleaning up...\n";
      } else {
        if (_modelUpdateRequired) PartitionedMS::unpartition(*this);
        Logger::Info << "Cleaning up temporary files...\n";
      }

      std::set<size_t> removedMetaFiles;
      for (size_t part = 0; part != _channels.size(); ++part) {
        for (aocommon::PolarizationEnum p : _polarizations) {
          std::string prefix =
              getPartPrefix(_msPath, part, p, _channels[part].dataDescId,
                            _temporaryDirectory);
          std::remove((prefix + ".tmp").c_str());
          std::remove((prefix + "-w.tmp").c_str());
          std::remove((prefix + "-m.tmp").c_str());
        }
        const size_t dataDescId = _channels[part].dataDescId;
        if (removedMetaFiles.count(dataDescId) == 0) {
          removedMetaFiles.insert(dataDescId);
          std::string metaFile =
              getMetaFilename(_msPath, _temporaryDirectory, dataDescId);
          std::remove(metaFile.c_str());
        }
      }
    } catch (std::exception& exception) {
      Logger::Error << "Error occurred while finishing IO task: "
                    << exception.what()
                    << "\nMeasurement set might not have been updated.\n";
    }
  }
}

void PartitionedMS::MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) {
  const MSSelection& selection = _handle._data->_selection;
  const std::map<size_t, size_t> dataDescIds =
      getDataDescIdMap(_handle._data->_channels);
  std::set<size_t> dataDescIdSet;
  for (std::map<size_t, size_t>::const_iterator i = dataDescIds.begin();
       i != dataDescIds.end(); ++i)
    dataDescIdSet.insert(i->first);
  size_t startRow, endRow;
  SynchronizedMS ms = MS();
  GetRowRangeAndIDMap(*ms, selection, startRow, endRow, dataDescIdSet,
                      idToMSRow);
}

std::map<size_t, size_t> PartitionedMS::getDataDescIdMap(
    const std::vector<PartitionedMS::ChannelRange>& channels) {
  std::map<size_t, size_t> dataDescIds;
  size_t spwIndex = 0;
  for (const PartitionedMS::ChannelRange& range : channels) {
    if (dataDescIds.count(range.dataDescId) == 0) {
      dataDescIds.emplace(range.dataDescId, spwIndex);
      ++spwIndex;
    }
  }
  return dataDescIds;
}

void PartitionedMS::Handle::Serialize(aocommon::SerialOStream& stream) const {
  stream.Ptr(_data);
}

void PartitionedMS::Handle::Unserialize(aocommon::SerialIStream& stream) {
  stream.Ptr(_data);
}

void PartitionedMS::Handle::HandleData::Serialize(
    aocommon::SerialOStream& stream) const {
  stream.String(_msPath)
      .String(_dataColumnName)
      .String(_temporaryDirectory)
      .UInt64(_channels.size());
  for (const ChannelRange& range : _channels) {
    stream.UInt64(range.dataDescId).UInt64(range.start).UInt64(range.end);
  }
  stream.Bool(_initialModelRequired)
      .Bool(_modelUpdateRequired)
      .UInt64(_polarizations.size());
  for (aocommon::PolarizationEnum p : _polarizations) stream.UInt32(p);
  _selection.Serialize(stream);
  stream.UInt64(_nAntennas);
}

void PartitionedMS::Handle::HandleData::Unserialize(
    aocommon::SerialIStream& stream) {
  _isCopy = true;
  stream.String(_msPath).String(_dataColumnName).String(_temporaryDirectory);
  _channels.resize(stream.UInt64());
  for (ChannelRange& range : _channels) {
    stream.UInt64(range.dataDescId).UInt64(range.start).UInt64(range.end);
  }
  stream.Bool(_initialModelRequired).Bool(_modelUpdateRequired);
  size_t nPol = stream.UInt64();
  _polarizations.clear();
  for (size_t i = 0; i != nPol; ++i)
    _polarizations.emplace((aocommon::PolarizationEnum)stream.UInt32());
  _selection.Unserialize(stream);
  stream.UInt64(_nAntennas);
}
