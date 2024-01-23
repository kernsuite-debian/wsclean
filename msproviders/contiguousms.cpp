#include "contiguousms.h"
#include "msreaders/contiguousmsreader.h"

#include <aocommon/logger.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/TableLocker.h>

#include <memory>

ContiguousMS::ContiguousMS(const string& msPath,
                           const std::string& dataColumnName,
                           const MSSelection& selection,
                           aocommon::PolarizationEnum outputPolarization,
                           size_t dataDescId, bool useMPI)
    : _currentOutputRow(0),
      _currentOutputTimestep(0),
      _currentOutputTime(0.0),
      _dataDescId(dataDescId),
      _useMPI(useMPI),
      _nAntenna(0),
      _isModelColumnPrepared(false),
      _selection(selection),
      _outputPolarization(outputPolarization),
      _msPath(msPath),
      _dataColumnName(dataColumnName) {
  open();
}

void ContiguousMS::open() {
  aocommon::Logger::Info << "Opening " << _msPath << ", spw " << _dataDescId
                         << " with contiguous MS reader.\n";

  _ms = SynchronizedMS(_msPath, _useMPI ? casacore::TableLock::UserNoReadLocking
                                        : casacore::TableLock::DefaultLocking);

  _antenna1Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  _antenna2Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  _fieldIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  _dataDescIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID));
  _timeColumn = casacore::ScalarColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  _uvwColumn = casacore::ArrayColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  _dataColumn = casacore::ArrayColumn<casacore::Complex>(*_ms, _dataColumnName);
  _flagColumn = casacore::ArrayColumn<bool>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG));

  _inputPolarizations = GetMSPolarizations(_dataDescId, *_ms);

  const casacore::IPosition shape(_dataColumn.shape(0));
  _dataArray = casacore::Array<std::complex<float>>(shape);
  _weightSpectrumArray = casacore::Array<float>(shape);
  _imagingWeightSpectrumArray = casacore::Array<float>(shape);
  _flagArray = casacore::Array<bool>(shape);
  _bandData = aocommon::MultiBandData(*_ms);

  if (_bandData.BandCount() > 1) {
    throw std::runtime_error(
        "This set contains multiple spws, and can therefore not be opened "
        "directly due to possible synchronization issues between spws. You can "
        "force reordering of the measurement by adding -reorder to the command "
        "line.");
  }
  _nAntenna = _ms->antenna().nrow();

  _msHasWeightSpectrum = OpenWeightSpectrumColumn(*_ms, _weightSpectrumColumn);
  if (!_msHasWeightSpectrum) {
    casacore::IPosition scalarShape(1, shape[0]);
    _weightScalarArray = casacore::Array<float>(scalarShape);
    _weightScalarColumn.reset(new casacore::ArrayColumn<float>(
        *_ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
  }

  GetRowRangeAndIDMap(*_ms, _selection, _startRow, _endRow,
                      std::set<size_t>{size_t(_dataDescId)}, _idToMSRow);
  ResetWritePosition();
}

std::unique_ptr<MSReader> ContiguousMS::MakeReader() {
  std::unique_ptr<MSReader> reader(new ContiguousMSReader(this));
  return reader;
}

void ContiguousMS::ResetWritePosition() {
  _currentOutputRow = _startRow - 1;
  _currentOutputTime = 0.0;
  // TODO: something similar needed in the ContiguousMSReader class?
  if (_selection.HasInterval())
    _currentOutputTimestep = _selection.IntervalStart() - 1;
  else
    _currentOutputTimestep = -1;
  NextOutputRow();
}

void ContiguousMS::NextOutputRow() {
  int fieldId, a1, a2, dataDescId;
  casacore::Vector<double> uvw;
  do {
    ++_currentOutputRow;
    if (_currentOutputRow >= _endRow) return;

    fieldId = _fieldIdColumn(_currentOutputRow);
    a1 = _antenna1Column(_currentOutputRow);
    a2 = _antenna2Column(_currentOutputRow);
    uvw = _uvwColumn(_currentOutputRow);
    dataDescId = _dataDescIdColumn(_currentOutputRow);
    if (_currentOutputTime != _timeColumn(_currentOutputRow)) {
      ++_currentOutputTimestep;
      _currentOutputTime = _timeColumn(_currentOutputRow);
    }
  } while (
      !_selection.IsSelected(fieldId, _currentOutputTimestep, a1, a2, uvw) ||
      (dataDescId != _dataDescId));
}

double ContiguousMS::StartTime() {
  return casacore::MEpoch::ScalarColumn(
             *_ms, casacore::MS::columnName(casacore::MS::TIME))(_startRow)
      .getValue()
      .get();
}

size_t ContiguousMS::NChannels() {
  if (_selection.HasChannelRange())
    return _selection.ChannelRangeEnd() - _selection.ChannelRangeStart();
  else
    return _bandData[_dataDescId].ChannelCount();
}

size_t ContiguousMS::NPolarizations() {
  switch (_outputPolarization) {
    case aocommon::Polarization::Instrumental:
      return 4;
    case aocommon::Polarization::DiagonalInstrumental:
      return 2;
    default:
      return 1;
  }
}

void ContiguousMS::prepareModelColumn() {
  InitializeModelColumn(*_ms);
  _modelColumn = casacore::ArrayColumn<casacore::Complex>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
  const casacore::IPosition shape(_modelColumn.shape(0));
  _modelArray = casacore::Array<std::complex<float>>(shape);
  _isModelColumnPrepared = true;
}

void ContiguousMS::WriteModel(const std::complex<float>* buffer, bool addToMS) {
  std::unique_ptr<casacore::TableLocker> lock;
  if (_useMPI) {
    // When using different MPI processes, automatic casacore locks do not work.
    // -> Use UserNoReadLocking with explicit write locks.
    lock = std::make_unique<casacore::TableLocker>(*_ms,
                                                   casacore::FileLocker::Write);

    // This resync() is required for synchronizing different MPI processes.
    _ms->resync();
  }

  if (!_isModelColumnPrepared) prepareModelColumn();

  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }

  _modelColumn.get(_currentOutputRow, _modelArray);
  if (addToMS) {
    ReverseCopyData<true>(_modelArray, startChannel, endChannel,
                          _inputPolarizations, buffer, _outputPolarization);
  } else {
    ReverseCopyData<false>(_modelArray, startChannel, endChannel,
                           _inputPolarizations, buffer, _outputPolarization);
  }
  _modelColumn.put(_currentOutputRow, _modelArray);
}

void ContiguousMS::MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) {
  idToMSRow = _idToMSRow;
}
