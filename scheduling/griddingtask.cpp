#include "griddingtask.h"

#include "../idg/averagebeam.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>
#include <schaapcommon/facets/facet.h>

GriddingTask::GriddingTask() = default;
GriddingTask::GriddingTask(GriddingTask&& source) noexcept = default;
GriddingTask::~GriddingTask() noexcept = default;
GriddingTask& GriddingTask::operator=(GriddingTask&& source) noexcept = default;

void GriddingTask::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(operation)
      .Bool(imagePSF)
      .Bool(subtractModel)
      .UInt32(polarization)
      .Bool(verbose)
      .Ptr(cache)
      .Ptr(averageBeam)
      .Bool(storeImagingWeights)
      .Ptr(imageWeights)
      .Ptr(facet)
      .UInt64(facetIndex)
      .UInt64(facetGroupIndex);

  // msList
  stream.UInt64(msList.size());
  for (const std::unique_ptr<MSDataDescription>& dataDesc : msList)
    dataDesc->Serialize(stream);

  stream.ObjectVector(modelImages)
      .Object(observationInfo)
      .Double(l_shift)
      .Double(m_shift);
}

void GriddingTask::Unserialize(aocommon::SerialIStream& stream) {
  operation = static_cast<Operation>(stream.UInt32());
  stream.Bool(imagePSF).Bool(subtractModel);
  polarization = static_cast<aocommon::PolarizationEnum>(stream.UInt32());
  stream.Bool(verbose)
      .Ptr(cache)
      .Ptr(averageBeam)
      .Bool(storeImagingWeights)
      .Ptr(imageWeights)
      .Ptr(facet)
      .UInt64(facetIndex)
      .UInt64(facetGroupIndex);

  // msList
  msList.resize(stream.UInt64());
  for (std::unique_ptr<MSDataDescription>& dataDesc : msList)
    dataDesc = MSDataDescription::Unserialize(stream);

  stream.ObjectVector(modelImages)
      .Object(observationInfo)
      .Double(l_shift)
      .Double(m_shift);
}
