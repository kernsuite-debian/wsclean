#ifndef MSPROVIDERS_MSREADERS_PARTITIONEDMSREADER_H
#define MSPROVIDERS_MSREADERS_PARTITIONEDMSREADER_H

#include "msreader.h"

#include <aocommon/uvector.h>

#include <fstream>

class PartitionedMS;

class PartitionedMSReader final : public MSReader {
 public:
  PartitionedMSReader(PartitionedMS* partitionedMS);
  virtual ~PartitionedMSReader(){};

  size_t RowId() const override { return _currentInputRow; }

  bool CurrentRowAvailable() override;

  void NextInputRow() override;

  void ReadMeta(double& u, double& v, double& w) override;

  void ReadMeta(MSProvider::MetaData& metaData) override;

  void ReadData(std::complex<float>* buffer) override;

  void ReadModel(std::complex<float>* buffer) override;

  void ReadWeights(float* buffer) override;

  void WriteImagingWeights(const float* buffer) override;

 private:
  size_t _currentInputRow;

  // Chunkoffset counts the amount of data rows we are ahead or behind.
  // Positive values mean we are lagging, whereas negative values mean we are
  //  ahead of the current time step.
  long _readPtrRowOffset, _metaPtrRowOffset, _weightPtrRowOffset;

  std::ifstream _metaFile, _weightFile, _dataFile;

  aocommon::UVector<float> _imagingWeightBuffer;
  std::unique_ptr<std::fstream> _imagingWeightsFile;
};

#endif
