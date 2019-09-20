#ifndef PARALLEL_DECONVOLUTION_H
#define PARALLEL_DECONVOLUTION_H

#include "../fftwmanager.h"
#include "../image.h"
#include "../uvector.h"

#include <memory>
#include <mutex>
#include <vector>

class ParallelDeconvolution
{
public:
	ParallelDeconvolution(const class WSCleanSettings& settings);
	
	~ParallelDeconvolution();
	
	class DeconvolutionAlgorithm& FirstAlgorithm()
	{
		return *_algorithms.front();
	}
	const class DeconvolutionAlgorithm& FirstAlgorithm() const
	{
		return *_algorithms.front();
	}
	
	void SetAllocator(class ImageBufferAllocator* allocator)
	{
		_allocator = allocator;
	}
	
	void SetAlgorithm(std::unique_ptr<class DeconvolutionAlgorithm> algorithm);
	
	void SetRMSFactorImage(class Image&& image);
	
	void SetThreshold(double threshold);
	
	bool IsInitialized() const
	{
		return !_algorithms.empty();
	}
	
	void SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks);
	
	void SetCleanMask(const bool* mask);
	
	void ExecuteMajorIteration(class ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, bool& reachedMajorThreshold);
	
	void FreeDeconvolutionAlgorithms()
	{
		_algorithms.clear(); 
		_mask = nullptr;
	}
	
	void SaveSourceList(class CachedImageSet& modelImages, const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec);
	
	void SavePBSourceList(class CachedImageSet& modelImages, const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const;
	
	class FFTWManager& GetFFTWManager() { return _fftwManager; }
	
private:
	void executeParallelRun(class ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, bool& reachedMajorThreshold);
	
	struct SubImage {
		size_t index, x, y, width, height;
		ao::uvector<bool> mask;
		double peak;
		bool reachedMajorThreshold;
	};
		
	void runSubImage(SubImage& subImg, ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, double majorIterThreshold, bool findPeakOnly, std::mutex* mutex);
	
	void correctChannelForPB(class ComponentList& list, const class ImagingTableEntry& entry) const;
	
	void loadAveragePrimaryBeam(class PrimaryBeamImageSet& beamImages, size_t imageIndex, const class ImagingTable& table) const;
	
	void writeSourceList(ComponentList& componentList, const std::string& filename, long double phaseCentreRA, long double phaseCentreDec) const;
	
	FFTWManager _fftwManager;
	std::vector<std::unique_ptr<class DeconvolutionAlgorithm>> _algorithms;
	size_t _horImages, _verImages;
	const WSCleanSettings& _settings;
	ImageBufferAllocator* _allocator;
	const bool* _mask;
	bool _trackPerScaleMasks, _usePerScaleMasks;
	std::vector<ao::uvector<bool>> _scaleMasks;
	std::unique_ptr<class ComponentList> _componentList;
	Image _rmsImage;
};

#endif

