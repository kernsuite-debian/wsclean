#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "gridmode.h"

#include <aocommon/banddata.h>
#include <aocommon/polarization.h>
#include <aocommon/imagecoordinates.h>

#include "../structures/observationinfo.h"
#include "../structures/msselection.h"
#include "../structures/image.h"
#include "../structures/weightmode.h"

#include "visibilityweightingmode.h"

#include "../main/settings.h"

#include "../scheduling/metadatacache.h"
#include "../scheduling/griddingtaskmanager.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <schaapcommon/h5parm/jonesparameters.h>

#include <aocommon/uvector.h>

#include <mutex>
#include <memory>

class MSReader;
namespace schaapcommon {
namespace h5parm {
class H5Parm;
class SolTab;
}  // namespace h5parm
}  // namespace schaapcommon

/**
 * Enum for selecting the entry or entries from the direction dependent gain
 * matrix that are to be used to correct the visibilities during the reading
 * and/or writing operations.
 */
enum class DDGainMatrix { kXX, kYY, kTrace, kFull };

class MSGridderBase {
 public:
  MSGridderBase(const Settings& settings);
  virtual ~MSGridderBase();

  size_t ImageWidth() const { return _imageWidth; }
  size_t ImageHeight() const { return _imageHeight; }
  double PixelSizeX() const { return _pixelSizeX; }
  double PixelSizeY() const { return _pixelSizeY; }
  size_t ActualWGridSize() const { return _actualWGridSize; }

  void ClearMeasurementSetList() {
    _measurementSets.clear();
    _selections.clear();
  }
  class MSProvider& MeasurementSet(size_t index) const {
    return *_measurementSets[index];
  }
  const MSSelection& Selection(size_t index) const {
    return _selections[index];
  }
  size_t MeasurementSetCount() const { return _measurementSets.size(); }
  void AddMeasurementSet(class MSProvider* msProvider,
                         const MSSelection& selection) {
    _measurementSets.push_back(msProvider);
    _selections.push_back(selection);
  }

  const std::string& DataColumnName() const { return _dataColumnName; }
  bool DoImagePSF() const { return _doImagePSF; }
  bool DoSubtractModel() const { return _doSubtractModel; }
  bool SmallInversion() const { return _smallInversion; }
  aocommon::PolarizationEnum Polarization() const { return _polarization; }
  WeightMode Weighting() const { return _weighting; }
  const class ImageWeights* GetImageWeights() const {
    return _precalculatedWeightInfo;
  }
  bool IsComplex() const { return _isComplex; }

  VisibilityWeightingMode GetVisibilityWeightingMode() const {
    return _visibilityWeightingMode;
  }
  bool StoreImagingWeights() const { return _storeImagingWeights; }

  void SetFacetIndex(size_t facetIndex) { _facetIndex = facetIndex; }
  void SetFacetGroupIndex(size_t index) { _facetGroupIndex = index; }
  /**
   * @brief In case of facet-based imaging, the model data in the @param
   * MSProvider is reset to zeros in every major cycle, and predicted data
   * should be add-assigned to the model data (_additivePredict = true) rather
   * than overwriting it. For "standard" imaging, the model data should
   * be overwritten (_additivePredict = false).
   */
  void SetAdditivePredict(bool hasFacets) { _additivePredict = hasFacets; }
  void SetImageWidth(size_t imageWidth) { _imageWidth = imageWidth; }
  void SetImageHeight(size_t imageHeight) { _imageHeight = imageHeight; }
  void SetActualWGridSize(size_t actualWGridSize) {
    _actualWGridSize = actualWGridSize;
  }
  void SetDoImagePSF(bool doImagePSF) { _doImagePSF = doImagePSF; }
  void SetPolarization(aocommon::PolarizationEnum polarization) {
    _polarization = polarization;
  }
  void SetIsComplex(bool isComplex) { _isComplex = isComplex; }
  void SetDoSubtractModel(bool doSubtractModel) {
    _doSubtractModel = doSubtractModel;
  }

  void SetWriterLockManager(WriterLockManager* writerLockManager) {
    _writerLockManager = writerLockManager;
  }

  void SetImageWeights(const class ImageWeights* weights) {
    _precalculatedWeightInfo = weights;
  }

  /**
   * If this is the first gridder iteration, the gridder may output more
   * information.
   */
  bool IsFirstIteration() const { return _isFirstIteration; }
  void SetIsFirstIteration(bool isFirstIteration) {
    _isFirstIteration = isFirstIteration;
  }

  void SetStoreImagingWeights(bool storeImagingWeights) {
    _storeImagingWeights = storeImagingWeights;
  }

  virtual void Invert() = 0;

  virtual void Predict(std::vector<Image>&& images) = 0;

  virtual std::vector<Image> ResultImages() = 0;

  void SetPhaseCentreRA(const double phaseCentreRA) {
    _phaseCentreRA = phaseCentreRA;
  }
  void SetPhaseCentreDec(const double phaseCentreDec) {
    _phaseCentreDec = phaseCentreDec;
  }
  double PhaseCentreRA() const { return _phaseCentreRA; }
  double PhaseCentreDec() const { return _phaseCentreDec; }
  void SetPhaseCentreDL(const double phaseCentreDL) {
    _phaseCentreDL = phaseCentreDL;
  }
  void SetPhaseCentreDM(const double phaseCentreDM) {
    _phaseCentreDM = phaseCentreDM;
  }

  void SetFacetDirectionRA(double facetDirectionRA) {
    _facetDirectionRA = facetDirectionRA;
  }
  void SetFacetDirectionDec(double facetDirectionDec) {
    _facetDirectionDec = facetDirectionDec;
  }

  double FacetDirectionRA() const { return _facetDirectionRA; }
  double FacetDirectionDec() const { return _facetDirectionDec; }
  double PhaseCentreDL() const { return _phaseCentreDL; }
  double PhaseCentreDM() const { return _phaseCentreDM; }

  /**
   * Deallocate any data that is no longer necessary, but all methods
   * will still return results from the imaging, with the exception of
   * ImageReal/ImageResult().
   */
  virtual void FreeImagingData() {}

  virtual size_t ActualInversionWidth() const { return _imageWidth; }
  virtual size_t ActualInversionHeight() const { return _imageHeight; }

  GridMode GetGridMode() const { return _gridMode; }
  void SetGridMode(GridMode gridMode) { _gridMode = gridMode; }

  size_t TrimWidth() const { return _trimWidth; }
  size_t TrimHeight() const { return _trimHeight; }
  bool HasTrimSize() const { return _trimWidth != 0 || _trimHeight != 0; }
  void SetTrimSize(size_t trimWidth, size_t trimHeight) {
    _trimWidth = trimWidth;
    _trimHeight = trimHeight;
  }

  double StartTime() const { return _startTime; }
  bool HasDenormalPhaseCentre() const {
    return _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
  }
  double ImageWeight() const { return _totalWeight; }
  double NormalizationFactor() const { return _totalWeight; }
  double BeamSize() const { return _theoreticalBeamSize; }

  /**
   * This is the sum of the weights as given by the measurement set, before the
   * image weighting is applied.
   */
  double VisibilityWeightSum() const { return _visibilityWeightSum; }
  /**
   * The number of visibilities that were gridded.
   */
  size_t GriddedVisibilityCount() const { return _griddedVisibilityCount; }
  /**
   * The maximum weight, after having applied the imaging weighting.
   */
  double MaxGriddedWeight() const { return _maxGriddedWeight; }
  /**
   * The effective number of visibilities, taking into account imaging weighting
   * and visibility weighting. This number is relative to the "best" visibility:
   * if one visibility with a weight of 10 and 5 visibilities with
   * a weight of 4 were gridded, the effective number of visibilities is
   * (10 + 5 x 4) / 10 = 3
   */
  double EffectiveGriddedVisibilityCount() const {
    return totalWeight() / MaxGriddedWeight();
  }

  void SetMetaDataCache(std::unique_ptr<MetaDataCache> cache) {
    _metaDataCache = std::move(cache);
  }
  std::unique_ptr<MetaDataCache> AcquireMetaDataCache() {
    return std::move(_metaDataCache);
  }

 protected:
  int64_t getAvailableMemory(double memFraction, double absMemLimit);

  struct MSData {
   public:
    MSData();
    MSData(const MSData& source) = delete;
    ~MSData() = default;
    MSData& operator=(const MSData& source) = delete;

    class MSProvider* msProvider;
    size_t msIndex;
    size_t dataDescId;
    aocommon::BandData bandData;
    size_t startChannel, endChannel;
    size_t matchingRows, totalRowsProcessed;
    double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM;
    size_t rowStart, rowEnd;
    double integrationTime;
    std::vector<std::string> antennaNames;

    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(bandData, startChannel, endChannel);
    }
  };

  struct InversionRow {
    double uvw[3];
    size_t rowId;
    std::complex<float>* data;
  };

  /**
   * Initializes MS related data members, i.e. the @param _telescope and the
   * @param _pointResponse data in case a beam is applied on the facets and
   * EveryBeam is available and the @param _predictReader data member in case
   * @param isPredict is true.
   */
  void StartMeasurementSet(const MSGridderBase::MSData& msData,
                           bool isPredict) {
    initializePointResponse(msData);
    _msIndex = msData.msIndex;
    if (isPredict) initializePredictReader(*msData.msProvider);
  }

  /**
   * Read the visibilities from the msprovider, and apply weights and flags.
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc). To
   * read the data, this function requires scratch weight and model buffers to
   * store intermediate values in. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   * @tparam PolarizationCount Normally set to one when imaging a single
   * polarization, but set to 4 for IDG as it images all polarizations at once.
   * @tparam DDGainMatrix Selects which entry or entries in the gain matrix
   * (provided by EveryBeam and/or an h5 solution) file to use for correcting
   * the visibilities. Can be kXX for the XX-entry, kYY for the YY-entry, kTrace
   * for the trace of the gain matrix and kFull to take all entries into
   * account.
   * @param msProvider The measurement set provider
   * @param rowData The resulting weighted data
   * @param curBand The spectral band currently being imaged
   * @param weightBuffer An allocated buffer to store intermediate weights in.
   * After returning from the call, these values will hold the full applied
   * weight (i.e. visibility weight * imaging weight).
   * @param modelBuffer An allocated buffer to store intermediate model data in.
   * @param isSelected Per visibility whether that visibility will be gridded in
   * this pass. When the visibility is not gridded, its weight will not be added
   * to the relevant sums (visibility count, weight sum, etc.).
   */
  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void readAndWeightVisibilities(MSReader& msReader,
                                 const std::vector<std::string>& antennaNames,
                                 InversionRow& rowData,
                                 const aocommon::BandData& curBand,
                                 float* weightBuffer,
                                 std::complex<float>* modelBuffer,
                                 const bool* isSelected);

  /**
   * @brief Write (modelled) visibilities to MS, provides an interface to
   * MSProvider::WriteModel(). Method can be templated on the number of
   * polarizations (1 or 4), and the DDGainMatrix which can be used to
   * select an entry or entries from the gain matrix that should be used for the
   * correction (XX-pol: kXX, YY-pol: kYY, Trace: kTrace, Full Jones: kFull)
   */
  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void writeVisibilities(MSProvider& msProvider,
                         const std::vector<std::string>& antennaNames,
                         const aocommon::BandData& curBand,
                         std::complex<float>* buffer);

  double _maxW, _minW;
  size_t _actualInversionWidth, _actualInversionHeight;
  double _actualPixelSizeX, _actualPixelSizeY;

  virtual size_t getSuggestedWGridSize() const = 0;

  void resetVisibilityCounters() {
    _griddedVisibilityCount = 0;
    _totalWeight = 0.0;
    _maxGriddedWeight = 0.0;
    _visibilityWeightSum = 0.0;
  }

  double totalWeight() const { return _totalWeight; }

  void initializeMSDataVector(std::vector<MSData>& msDataVector);

  std::unique_ptr<struct MetaDataCache> _metaDataCache;

  template <size_t PolarizationCount>
  static void rotateVisibilities(const aocommon::BandData& bandData,
                                 double shiftFactor,
                                 std::complex<float>* dataIter);

  const Settings& _settings;

 private:
  static std::vector<std::string> getAntennaNames(
      const casacore::MSAntenna& msAntenna);

  void resetMetaData() { _hasFrequencies = false; }

  void calculateMSLimits(const aocommon::BandData& selectedBand,
                         double startTime) {
    if (_hasFrequencies) {
      _freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
      _freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
      _bandStart = std::min(_bandStart, selectedBand.BandStart());
      _bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
      _startTime = std::min(_startTime, startTime);
    } else {
      _freqLow = selectedBand.LowestFrequency();
      _freqHigh = selectedBand.HighestFrequency();
      _bandStart = selectedBand.BandStart();
      _bandEnd = selectedBand.BandEnd();
      _startTime = startTime;
      _hasFrequencies = true;
    }
  }

  template <size_t NPolInMSProvider>
  void calculateWLimits(MSGridderBase::MSData& msData);

  void initializeMeasurementSet(MSGridderBase::MSData& msData,
                                MetaDataCache::Entry& cacheEntry,
                                bool isCacheInitialized);

  void calculateOverallMetaData(const MSData* msDataVector);
  bool hasWGridSize() const { return _wGridSize != 0; }
  void initializeBandData(const casacore::MeasurementSet& ms,
                          MSGridderBase::MSData& msData);
  void initializePointResponse(const MSGridderBase::MSData& msData);
  void initializePredictReader(MSProvider& msProvider);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void ApplyConjugatedFacetBeam(MSReader& msReader, InversionRow& rowData,
                                const aocommon::BandData& curBand,
                                float* weightBuffer);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void ApplyConjugatedH5Parm(MSReader& msReader,
                             const std::vector<std::string>& antennaNames,
                             InversionRow& rowData,
                             const aocommon::BandData& curBand,
                             float* weightBuffer);

  double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
  double _facetDirectionRA, _facetDirectionDec;
  size_t _facetIndex;
  /// @p _facetGroupIndex and @p _msIndex in conjunction with the @p
  /// MeasurementSetCount() determine the index in the _writerGroupLocks vector,
  /// having size FacetGroupCount() * MeasurementSetCount(). These variable are
  /// only relevant for prediction.
  size_t _facetGroupIndex;
  size_t _msIndex;
  /// @see SetAdditivePredict()
  bool _additivePredict;
  size_t _imageWidth, _imageHeight;
  size_t _trimWidth, _trimHeight;
  double _pixelSizeX, _pixelSizeY;
  size_t _wGridSize, _actualWGridSize;
  std::vector<MSProvider*> _measurementSets;
  std::string _dataColumnName;
  bool _doImagePSF, _doSubtractModel, _smallInversion;
  double _wLimit;
  const class ImageWeights* _precalculatedWeightInfo;
  aocommon::PolarizationEnum _polarization;
  bool _isComplex;
  WeightMode _weighting;
  bool _isFirstIteration;
  std::vector<MSSelection> _selections;
  VisibilityWeightingMode _visibilityWeightingMode;
  GridMode _gridMode;
  bool _storeImagingWeights;
  double _theoreticalBeamSize;

  bool _hasFrequencies;
  double _freqHigh, _freqLow;
  double _bandStart, _bandEnd;
  double _startTime;

  size_t _griddedVisibilityCount;
  double _totalWeight;
  double _maxGriddedWeight;
  double _visibilityWeightSum;

  aocommon::UVector<float> _scratchWeights;

  std::unique_ptr<MSReader> _predictReader;
  WriterLockManager* _writerLockManager;

#ifdef HAVE_EVERYBEAM
  // _telescope attribute needed to keep the telecope in _pointResponse alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  aocommon::UVector<std::complex<float>> _cachedBeamResponse;
  const everybeam::BeamMode _beamMode;
  const everybeam::BeamNormalisationMode _beamNormalisationMode;
#endif
  std::vector<std::vector<std::complex<float>>> _cachedParmResponse;
  std::vector<std::unique_ptr<schaapcommon::h5parm::H5Parm>> _h5parms;
  std::vector<
      std::pair<schaapcommon::h5parm::SolTab*, schaapcommon::h5parm::SolTab*>>
      _h5SolTabs;
  std::vector<schaapcommon::h5parm::JonesParameters::CorrectType> _correctType;
  std::vector<std::vector<double>> _cachedMSTimes;
  std::vector<size_t> _timeOffset;
};

#endif
