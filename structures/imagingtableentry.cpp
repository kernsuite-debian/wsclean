#include "imagingtableentry.h"

ImagingTableEntry::ImagingTableEntry()
    : index(0),
      lowestFrequency(0.0),
      highestFrequency(0.0),
      bandStartFrequency(0.0),
      bandEndFrequency(0.0),
      inputChannelCount(0),
      polarization(aocommon::Polarization::StokesI),
      facetGroupIndex(0),
      facetIndex(0),
      facet(nullptr),
      centreShiftX(0),
      centreShiftY(0),
      outputChannelIndex(0),
      outputIntervalIndex(0),
      msData(),
      squaredDeconvolutionIndex(0),
      joinedGroupIndex(0),
      imageCount(0),
      tmpFilePrefix(),
      imageWeight(0.0) {}
