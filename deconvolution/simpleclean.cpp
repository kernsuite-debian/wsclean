#include "simpleclean.h"

#include "../imagecoordinates.h"

#include "../lane.h"
#include "../areaset.h"

#include <boost/thread/thread.hpp>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include <iostream>
#include <limits>

double SimpleClean::FindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = width * height;
	
	const size_t horBorderSize = round(width*borderRatio), verBorderSize = round(height*borderRatio);
	size_t xiStart = horBorderSize, xiEnd = width - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(endY, height - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;

	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index = yi*width + xiStart;
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = image[index];
			if(std::isfinite(value))
			{
				if(allowNegativeComponents) value = std::fabs(value);
				if(value > peakMax)
				{
					peakIndex = index;
					peakMax = std::fabs(value);
				}
			}
			++value;
			++index;
		}
	}
	if(peakIndex == width * height)
	{
		x = width; y = height;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else {
		x = peakIndex % width;
		y = peakIndex / width;
		return image[x + y*width];
	}
}

double SimpleClean::FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const bool* cleanMask, double borderRatio)
{
	double peakMax = std::numeric_limits<double>::min();
	x = width; y = height;
	
	const size_t horBorderSize = round(width*borderRatio), verBorderSize = round(height*borderRatio);
	size_t xiStart = horBorderSize, xiEnd = width - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(endY, height - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	
	for(size_t yi=startY; yi!=endY; ++yi)
	{
		const double *imgIter = &image[yi*width+xiStart];
		const bool* cleanMaskPtr = &cleanMask[yi*width+xiStart];
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = *imgIter;
			if(std::isfinite(value))
			{
				if(allowNegativeComponents) value = std::fabs(value);
				if(value > peakMax && *cleanMaskPtr)
				{
					x = xi;
					y = yi;
					peakMax = std::fabs(value);
				}
			}
			++imgIter;
			++cleanMaskPtr;
		}
	}
	if(y == height)
		return std::numeric_limits<double>::quiet_NaN();
	else
		return image[x + y*width];
}

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
template<bool AllowNegativeComponent>
double SimpleClean::FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, double borderRatio)
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = 0;
	
	__m256d mPeakMax = _mm256_set1_pd(peakMax);
	
	const size_t horBorderSize = floor(width*borderRatio), verBorderSize = floor(height*borderRatio);
	size_t xiStart = horBorderSize, xiEnd = width - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(endY, height - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index = yi*width + xiStart;
		const double* const endPtr = image + yi*width + xiEnd - 4;
		const double *i=image + index;
		for(; i<endPtr; i+=4)
		{
			__m256d val = _mm256_loadu_pd(i);
			if(AllowNegativeComponent) {
				__m256d negVal = _mm256_sub_pd(_mm256_set1_pd(0.0), val);
				val = _mm256_max_pd(val, negVal);
			}
			int mask = _mm256_movemask_pd(_mm256_cmp_pd(val, mPeakMax, _CMP_GT_OQ));
			if(mask != 0)
			{
				for(size_t di=0; di!=4; ++di)
				{
					double value = i[di];
					if(AllowNegativeComponent) value = std::fabs(value);
					if(value > peakMax)
					{
						peakIndex = index+di;
						peakMax = std::fabs(i[di]);
						mPeakMax = _mm256_set1_pd(peakMax);
					}
				}
			}
			index+=4;
		}
		for(; i!=endPtr+4; ++i)
		{
			double value = *i;
			if(AllowNegativeComponent) value = std::fabs(value);
			if(value > peakMax)
			{
				peakIndex = index;
				peakMax = std::fabs(*i);
			}
			++index;
		}
	}
	x = peakIndex % width;
	y = peakIndex / width;
	return image[x + y*width];
}

template
double SimpleClean::FindPeakAVX<false>(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, double borderRatio);
template
double SimpleClean::FindPeakAVX<true>(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, double borderRatio);
#else
#warning "Not using AVX optimized version of FindPeak()!"
#endif // __AVX__

void SimpleClean::SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor)
{
	size_t startX, startY, endX, endY;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > 0)
		startY = offsetY;
	else
		startY = 0;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = y + height/2;
	if(endY > height) endY = height;
	
	for(size_t ypos = startY; ypos != endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos++)
		{
			// I've SSE-ified this, but it didn't improve speed at all :-/
			// (Compiler probably already did it)
			*imageIter -= (*psfIter * factor);
			//*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			++imageIter;
			++psfIter;
		}
	}
}

void SimpleClean::PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + height/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			*imageIter = *imageIter - (*psfIter * factor);
			*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}

void SimpleClean::PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - psfWidth/2, offsetY = (int) y - psfHeight/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = std::min(x + psfWidth/2, imgWidth);
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + psfHeight/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * imgWidth + startX;
		const double *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			*imageIter = *imageIter - (*psfIter * factor);
			*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}

#if defined __AVX__ && defined USE_INTRINSICS
void SimpleClean::PartialSubtractImageAVX(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - psfWidth/2, offsetY = (int) y - psfHeight/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = std::min(x + psfWidth/2, imgWidth);
	
	size_t unAlignedCount = (endX - startX) % 4;
	endX -= unAlignedCount;
	
	endY = std::min(y + psfHeight/2, endY);
	
	const __m256d mFactor = _mm256_set1_pd(-factor);
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * imgWidth + startX;
		const double *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=4)
		{
			__m256d
				imgVal = _mm256_loadu_pd(imageIter),
				psfVal = _mm256_loadu_pd(psfIter);
#ifdef __AVX2__
			_mm256_storeu_pd(imageIter, _mm256_fmadd_pd(psfVal, mFactor, imgVal));
#else
			_mm256_storeu_pd(imageIter, _mm256_add_pd(imgVal, _mm256_mul_pd(psfVal, mFactor)));
#endif
			imageIter+=4;
			psfIter+=4;
		}
		for(size_t xpos = endX; xpos!=endX + unAlignedCount; ++xpos)
		{
			*imageIter -= *psfIter * factor;
			++imageIter;
			++psfIter;
		}
	}
}
#endif

void SimpleClean::ExecuteMajorIterationST(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height)
{
	size_t componentX, componentY;
	double peak = FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents, 0, height, CleanBorderRatio());
	Logger::Info << "Initial peak: " << peak << '\n';
	while(fabs(peak) > _threshold && _iterationNumber < _maxIter)
	{
		if(_iterationNumber % 10 == 0)
			Logger::Info << "Iteration " << _iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		SubtractImage(dataImage, psfImage, width, height, componentX, componentY, _gain * peak);
		modelImage[componentX + componentY*width] += _gain * peak;
		
		peak = FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents, 0, height, CleanBorderRatio());
		++_iterationNumber;
	}
	Logger::Info << "Stopped on peak " << peak << '\n';
}

void SimpleClean::ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain)
{
	if(_stopOnNegativeComponent)
		_allowNegativeComponents = true;
	
	size_t componentX=0, componentY=0;
	double peak = _cleanMask==0 ?
		FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents, 0, height, CleanBorderRatio()) :
		FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents, 0, height, _cleanMask, CleanBorderRatio());
	Logger::Info << "Initial peak: " << peak << '\n';
	double firstThreshold = _threshold, stopGainThreshold = fabs(peak*(1.0-_mGain));
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		Logger::Info << "Next major iteration at: " << stopGainThreshold << '\n';
	}
	else if(_mGain != 1.0) {
		Logger::Info << "Major iteration threshold reached global threshold of " << _threshold << ": final major iteration.\n";
	}

	std::vector<ao::lane<CleanTask>*> taskLanes(_threadCount);
	std::vector<ao::lane<CleanResult>*> resultLanes(_threadCount);
	boost::thread_group threadGroup;
	for(size_t i=0; i!=_threadCount; ++i)
	{
		taskLanes[i] = new ao::lane<CleanTask>(1);
		resultLanes[i] = new ao::lane<CleanResult>(1);
		CleanThreadData cleanThreadData;
		cleanThreadData.imgWidth = width;
		cleanThreadData.imgHeight = height;
		cleanThreadData.dataImage = dataImage;
		cleanThreadData.psfWidth = width;
		cleanThreadData.psfHeight = height;
		cleanThreadData.psfImage = psfImage;
		cleanThreadData.startY = (height*i)/_threadCount;
		cleanThreadData.endY = height*(i+1)/_threadCount;
		threadGroup.add_thread(new boost::thread(&SimpleClean::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	while(fabs(peak) > firstThreshold && _iterationNumber < _maxIter && (peak >= 0.0 || !_stopOnNegativeComponent))
	{
		if(
			(_iterationNumber <= 100 && _iterationNumber % 10 == 0) ||
			(_iterationNumber <= 1000 && _iterationNumber % 100 == 0) ||
			_iterationNumber % 1000 == 0)
			Logger::Info << "Iteration " << _iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		
		CleanTask task;
		task.cleanCompX = componentX;
		task.cleanCompY = componentY;
		task.peakLevel = peak;
		for(size_t i=0; i!=_threadCount; ++i)
			taskLanes[i]->write(task);
		
		modelImage[componentX + componentY*width] += _gain * peak;
		
		peak = 0.0;
		for(size_t i=0; i!=_threadCount; ++i)
		{
			CleanResult result;
			resultLanes[i]->read(result);
			if(std::isfinite(result.peakLevel) && fabs(result.peakLevel) >= fabs(peak))
			{
				peak = result.peakLevel;
				componentX = result.nextPeakX;
				componentY = result.nextPeakY;
			}
		}
		
		++_iterationNumber;
	}
	for(size_t i=0; i!=_threadCount; ++i)
		taskLanes[i]->write_end();
	threadGroup.join_all();
	for(size_t i=0; i!=_threadCount; ++i)
	{
		delete taskLanes[i];
		delete resultLanes[i];
	}
	Logger::Info << "Stopped on peak " << peak << '\n';
	reachedStopGain = (fabs(peak) <= stopGainThreshold) && (peak != 0.0);
}

void SimpleClean::cleanThreadFunc(ao::lane<CleanTask> *taskLane, ao::lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	while(taskLane->read(task))
	{
		PartialSubtractImage(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, cleanData.psfImage, cleanData.psfWidth, cleanData.psfHeight, task.cleanCompX, task.cleanCompY, _gain * task.peakLevel, cleanData.startY, cleanData.endY);
		
		CleanResult result;
		if(_cleanMask == 0)
			result.peakLevel = FindPeak(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY, CleanBorderRatio());
		else
			result.peakLevel = FindPeak(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY, _cleanMask, CleanBorderRatio());
		
		resultLane->write(result);
	}
}

