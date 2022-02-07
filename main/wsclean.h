#ifndef WSCLEAN_H
#define WSCLEAN_H

#include <aocommon/polarization.h>
#include <schaapcommon/facets/facet.h>

#include "../deconvolution/deconvolution.h"

#include "../scheduling/griddingresult.h"

#include "../io/cachedimageset.h"
#include "../io/wscfitswriter.h"

#include "../structures/imagingtable.h"
#include "../structures/msselection.h"
#include "../structures/multibanddata.h"
#include "../structures/observationinfo.h"
#include "../structures/outputchannelinfo.h"
#include "../structures/weightmode.h"

#include "../gridding/msgridderbase.h"

#include "../msproviders/partitionedms.h"

#include "stopwatch.h"
#include "settings.h"

#include <set>

namespace schaapcommon {
namespace facets {
class FacetImage;
}
}  // namespace schaapcommon

class PrimaryBeam;
class WSClean {
 public:
  WSClean();
  ~WSClean();

  Settings& GetSettings() { return _settings; }
  const Settings& GetSettings() const { return _settings; }
  void ResetSettings() { _settings = Settings(); }
  void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }

  void RunClean();

  /**
   * Entry point for performing a single prediction for an existing model image.
   *
   * In case of a facet-based prediction, the provided model images are assumed
   * to have the same size, so that the image size of the full image can be
   * inferred from the first entry in the _imagingTable in an early stage.
   */
  void RunPredict();

 private:
  void runIndependentGroup(ImagingTable& groupTable,
                           std::unique_ptr<PrimaryBeam>& primaryBeam);
  void saveRestoredImagesForGroup(
      const ImagingTable& table,
      std::unique_ptr<PrimaryBeam>& primaryBeam) const;
  void predictGroup(const ImagingTable& groupTable);

  void runFirstInversions(ImagingTable& groupTable,
                          std::unique_ptr<PrimaryBeam>& primaryBeam,
                          bool requestPolarizationsAtOnce,
                          bool parallelizePolarizations);
  void runSingleFirstInversion(ImagingTableEntry& entry,
                               std::unique_ptr<PrimaryBeam>& primaryBeam);
  void runMajorIterations(ImagingTable& groupTable,
                          std::unique_ptr<PrimaryBeam>& primaryBeam,
                          bool requestPolarizationsAtOnce,
                          bool parallelizePolarizations);

  void performReordering(bool isPredictMode);

  ObservationInfo getObservationInfo() const;
  /**
   * Add the phase shift of a facet to an ObservationInfo object.
   * @param entry entry. If its facet is null, nothing happens.
   * @param observationInfo shiftL and shiftM of this object are updated.
   */
  void applyFacetPhaseShift(const ImagingTableEntry& entry,
                            ObservationInfo& observationInfo) const;
  std::shared_ptr<ImageWeights> initializeImageWeights(
      const ImagingTableEntry& entry,
      std::vector<std::unique_ptr<class MSDataDescription>>& msList);
  void initializeMFSImageWeights();
  void initializeMSList(
      const ImagingTableEntry& entry,
      std::vector<std::unique_ptr<MSDataDescription>>& msList);
  void resetModelColumns(const ImagingTable& groupTable);
  void resetModelColumns(const ImagingTableEntry& entry);
  void storeAndCombineXYandYX(CachedImageSet& dest, size_t joinedChannelIndex,
                              const ImagingTableEntry& entry,
                              aocommon::PolarizationEnum polarization,
                              bool isImaginary, const float* image);
  bool selectChannels(MSSelection& selection, size_t msIndex, size_t dataDescId,
                      const ImagingTableEntry& entry);
  MSSelection selectInterval(MSSelection& fullSelection, size_t intervalIndex);

  void makeImagingTable(size_t outputIntervalIndex);
  void makeImagingTableEntry(const std::vector<aocommon::ChannelInfo>& channels,
                             size_t outIntervalIndex, size_t outChannelIndex,
                             ImagingTableEntry& entry);
  void makeImagingTableEntryChannelSettings(
      const std::vector<aocommon::ChannelInfo>& channels,
      size_t outIntervalIndex, size_t outChannelIndex, size_t nOutChannels,
      ImagingTableEntry& entry);
  void addFacetsToImagingTable(ImagingTableEntry& templateEntry);
  void addPolarizationsToImagingTable(ImagingTableEntry& templateEntry);
  std::unique_ptr<class ImageWeightCache> createWeightCache();

  void multiplyImage(double factor, double* image) const;

  /**
   * Initializes full-size model images for the given entry. Depending on the
   * settings, this might load existing images from disk or initialize
   * them to zero.
   */
  void initializeModelImages(const ImagingTableEntry& entry,
                             aocommon::PolarizationEnum polarization,
                             size_t nFacetGroups);
  void readExistingModelImages(const ImagingTableEntry& entry,
                               aocommon::PolarizationEnum polarization,
                               size_t nFacetGroups);
  /**
   * Override the image settings given a FitsReader object.
   * The boolean return value indicates whether the gridder needs
   * to be reset.
   */
  bool overrideImageSettings(const FitsReader& reader);
  GriddingResult loadExistingImage(ImagingTableEntry& entry, bool isPSF);
  void loadExistingPSF(ImagingTableEntry& entry);
  void loadExistingDirty(ImagingTableEntry& entry, bool updateBeamInfo);

  void imagePSF(ImagingTableEntry& entry);
  void imagePSFCallback(ImagingTableEntry& entry, struct GriddingResult& result,
                        bool writeBeamImage);

  void imageMain(ImagingTableEntry& entry, bool isFirstInversion,
                 bool updateBeamInfo);
  void imageMainCallback(ImagingTableEntry& entry,
                         struct GriddingResult& result, bool updateBeamInfo,
                         bool isInitialInversion);

  void predict(const ImagingTableEntry& entry);

  void saveUVImage(const Image& image, const ImagingTableEntry& entry,
                   bool isImaginary, const std::string& prefix) const;

  void processFullPSF(Image& image, const ImagingTableEntry& entry);

  /**
   * @brief Stitch facets for all FacetGroups
   *
   * @param table Imaging table
   * @param cachedImage CachedImages
   * @param writeDirty Write dirty image?
   * @param writePSF Write PSF image?
   */
  void stitchFacets(const ImagingTable& table, CachedImageSet& imageCache,
                    bool writeDirty, bool isPSF);

  /**
   * @brief Stitch facet for a single (Facet)Group
   */
  void stitchSingleGroup(const ImagingTable& facetGroup, size_t imageIndex,
                         CachedImageSet& imageCache, bool writeDirty,
                         bool isPSF, Image& fullImage,
                         schaapcommon::facets::FacetImage& facetImage,
                         size_t nFacetGroups);
  /**
   * Partition model image into facets and save them into fits files
   */
  void partitionModelIntoFacets(const ImagingTable& table, bool isPredictOnly);

  /**
   * Partition image into facets for a single (Facet)Group
   */
  void partitionSingleGroup(const ImagingTable& facetGroup, size_t imageIndex,
                            CachedImageSet& imageCache, const Image& fullImage,
                            schaapcommon::facets::FacetImage& facetImage,
                            bool isPredictOnly);

  void writeFirstResidualImages(const ImagingTable& groupTable) const;
  void writeModelImages(const ImagingTable& groupTable) const;

  double minTheoreticalBeamSize(const ImagingTable& table) const;

  void makeBeam();

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    bool isImaginary, bool isModel,
                                    bool isFullImage) const;

  WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry,
                                    aocommon::PolarizationEnum polarization,
                                    bool isImaginary, bool isModel,
                                    bool isFullImage) const;
  /**
   * @brief Apply the H5 solution to the (restored) image and save as -pb.fits
   * file. Method is only invoked in case no beam corrections are applied.
   */
  void correctImagesH5(aocommon::FitsWriter& writer, const ImagingTable& table,
                       const ImageFilename& imageName,
                       const std::string& filenameKind) const;

  /**
   * @brief Compute the total amount of MSProviders that will be generated.
   * This number is needed to initialize the writer locks in the prediction
   * tasks, which are set via a call to _griddingTaskManager->Start(). The
   * number of @p MSProviders is the acummulated number of bands per MS.
   *
   * @return size_t Number of MSProviders
   */
  size_t getMaxNrMSProviders() const {
    size_t msCount = 0;
    for (const auto& msBand : _msBands) msCount += msBand.DataDescCount();
    return msCount;
  }

  MSSelection _globalSelection;
  std::string _commandLine;

  Settings _settings;

  std::vector<OutputChannelInfo> _infoPerChannel;
  OutputChannelInfo _infoForMFS;
  std::map<size_t, std::unique_ptr<MetaDataCache>> _msGridderMetaCache;

  std::unique_ptr<class GriddingTaskManager> _griddingTaskManager;
  std::unique_ptr<class ImageWeightCache> _imageWeightCache;
  Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
  bool _isFirstInversion;
  size_t _majorIterationNr;
  CachedImageSet _psfImages, _modelImages, _residualImages;
  std::vector<PartitionedMS::Handle> _partitionedMSHandles;
  std::vector<MultiBandData> _msBands;
  Deconvolution _deconvolution;
  ImagingTable _imagingTable;
  ObservationInfo _observationInfo;
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> _facets;
  double _lastStartTime;
};

#endif
