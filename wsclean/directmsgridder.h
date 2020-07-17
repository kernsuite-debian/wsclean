#ifndef DIRECT_MS_GRIDDER_H
#define DIRECT_MS_GRIDDER_H

#include <aocommon/lane.h>

#include "msgridderbase.h"

template<typename num_t>
class DirectMSGridder : public MSGridderBase
{
public:
	const static size_t num_t_factor = (sizeof(num_t) + sizeof(double) - 1) / sizeof(double);
	
	DirectMSGridder(size_t nThreads);

	virtual void Invert() final override;
	
	virtual void Predict(Image image) final override;
	virtual void Predict(Image /*real*/, Image /*imaginary*/) final override
	{
		throw std::runtime_error("Direct FT imager can not predict complex images");
	}
	
	virtual Image ImageRealResult() final override { return std::move(_image); }
	virtual Image ImageImaginaryResult() final override {
		throw std::runtime_error("Direct FT imager can not make complex images");
	}
	virtual size_t getSuggestedWGridSize() const override final { return 1; }
	
private:
	struct InversionSample {
		num_t uInLambda, vInLambda, wInLambda;
		std::complex<float> sample;
	};
	size_t _nThreads;
	Image _image;
	num_t* _sqrtLMTable;
	std::vector<num_t*> _layers;
	aocommon::Lane<InversionSample> _inversionLane;
	
	void invertMeasurementSet(const MSData& msData, class ProgressBar& progress, size_t msIndex);
	void inversionWorker(size_t layer);
	void gridSample(const InversionSample& sample, size_t layer);
	void initializeSqrtLMLookupTable();
	
	num_t* allocate()
	{
		return new num_t[ImageWidth() * ImageHeight()];
	}
	void freeImg(num_t* ptr)
	{
		delete[] ptr;
	}
};

#endif
