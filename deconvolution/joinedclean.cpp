#include "joinedclean.h"

#include "../lane.h"

#include <boost/thread/thread.hpp>

template<typename ImageSetType>
void JoinedClean<ImageSetType>::ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedStopGain)
{
	if(this->_stopOnNegativeComponent)
		this->_allowNegativeComponents = true;
	_width = width;
	_height = height;
	_curPeakValues.resize(dataImage.ImageCount());
	
	size_t componentX=0, componentY=0;
	findPeak(dataImage, componentX, componentY);
	Logger::Info << "Initial peak: " << peakDescription(dataImage, componentX, componentY) << '\n';
	
	size_t peakIndex = componentX + componentY*_width;
	double peakNormalized = dataImage.JoinedValueNormalized(peakIndex);
	double firstThreshold = this->_threshold;
	double stopGainThreshold = std::fabs(peakNormalized)*(1.0-this->_mGain);
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		Logger::Info << "Next major iteration at: " << stopGainThreshold << '\n';
	}
	else if(this->_mGain != 1.0) {
		Logger::Info << "Major iteration threshold reached global threshold of " << this->_threshold << ": final major iteration.\n";
	}

	size_t cpuCount = this->_threadCount;
	std::vector<ao::lane<CleanTask>*> taskLanes(cpuCount);
	std::vector<ao::lane<CleanResult>*> resultLanes(cpuCount);
	boost::thread_group threadGroup;
	for(size_t i=0; i!=cpuCount; ++i)
	{
		taskLanes[i] = new ao::lane<CleanTask>(1);
		resultLanes[i] = new ao::lane<CleanResult>(1);
		CleanThreadData cleanThreadData;
		cleanThreadData.dataImage = &dataImage;
		cleanThreadData.psfImages = psfImages;
		cleanThreadData.startY = (height*i)/cpuCount;
		cleanThreadData.endY = height*(i+1)/cpuCount;
		threadGroup.add_thread(new boost::thread(&JoinedClean::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	
	while(fabs(peakNormalized) > firstThreshold && this->_iterationNumber < this->_maxIter && !(dataImage.IsComponentNegative(peakIndex) && this->_stopOnNegativeComponent))
	{
		if(this->_iterationNumber <= 10 ||
			(this->_iterationNumber <= 100 && this->_iterationNumber % 10 == 0) ||
			(this->_iterationNumber <= 1000 && this->_iterationNumber % 100 == 0) ||
			this->_iterationNumber % 1000 == 0)
			Logger::Info << "Iteration " << this->_iterationNumber << ": " << peakDescription(dataImage, componentX, componentY) << '\n';
		
		CleanTask task;
		task.cleanCompX = componentX;
		task.cleanCompY = componentY;
		typename ImageSetType::Value peakValues = dataImage.Get(peakIndex);
		
		for(size_t i=0; i!=dataImage.ImageCount(); ++i)
			_curPeakValues[i] = peakValues.GetValue(i);
		
		this->PerformSpectralFit(_curPeakValues.data());
		
		for(size_t i=0; i!=dataImage.ImageCount(); ++i)
			_curPeakValues[i] *= this->_gain;
			
		for(size_t i=0; i!=cpuCount; ++i)
			taskLanes[i]->write(task);
		
		modelImage.AddComponent(peakIndex, _curPeakValues.data());
		
		double peakUnnormalized = 0.0;
		for(size_t i=0; i!=cpuCount; ++i)
		{
			CleanResult result;
			resultLanes[i]->read(result);
			if(std::isfinite(result.peakLevelUnnormalized) && result.peakLevelUnnormalized >= peakUnnormalized)
			{
				peakUnnormalized = result.peakLevelUnnormalized;
				componentX = result.nextPeakX;
				componentY = result.nextPeakY;
			}
		}
		peakIndex = componentX + componentY*_width;
		peakNormalized = dataImage.JoinedValueNormalized(peakIndex);
		
		++this->_iterationNumber;
	}
	for(size_t i=0; i!=cpuCount; ++i)
		taskLanes[i]->write_end();
	threadGroup.join_all();
	for(size_t i=0; i!=cpuCount; ++i)
	{
		delete taskLanes[i];
		delete resultLanes[i];
	}
	Logger::Info << "Stopped on peak " << peakNormalized << '\n';
	reachedStopGain = std::fabs(peakNormalized) <= stopGainThreshold && (peakNormalized != 0.0);
}

template<typename ImageSetType>
void JoinedClean<ImageSetType>::findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY) const
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = _width * _height;
	
	const size_t
		horBorderSize = floor(_width*this->CleanBorderRatio()),
		verBorderSize = floor(_height*this->CleanBorderRatio());
	size_t xiStart = horBorderSize, xiEnd = _width - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(stopY, _height - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index=yi*_width + xiStart;
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = image.AbsJoinedValue(index);
			if(std::isfinite(value))
			{
				if(value > peakMax)
				{
					peakIndex = index;
					peakMax = value;
				}
			}
			++index;
		}
	}
	if(peakIndex == _width * _height)
	{
		x = _width;
		y = _height;
	}
	else {
		x = peakIndex % _width;
		y = peakIndex / _width;
	}
}

template<typename ImageSetType>
void JoinedClean<ImageSetType>::findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY, const bool* mask) const
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = _width * _height;
	
	const size_t
		horBorderSize = floor(_width*this->CleanBorderRatio()),
		verBorderSize = floor(_height*this->CleanBorderRatio());
	size_t xiStart = horBorderSize, xiEnd = _width - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(stopY, _height - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index=yi*_width + xiStart;
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			if(mask[index])
			{
				double value = image.AbsJoinedValue(index);
				if(std::isfinite(value))
				{
					if(value > peakMax)
					{
						peakIndex = index;
						peakMax = value;
					}
				}
			}
			++index;
		}
	}
	if(peakIndex == _width * _height)
	{
		x = _width;
		y = _height;
	}
	else {
		x = peakIndex % _width;
		y = peakIndex / _width;
	}
}

template<typename ImageSetType>
void JoinedClean<ImageSetType>::cleanThreadFunc(ao::lane<CleanTask> *taskLane, ao::lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	while(taskLane->read(task))
	{
		for(size_t i=0; i!=cleanData.dataImage->ImageCount(); ++i)
		{
			subtractImage(cleanData.dataImage->GetImage(i), cleanData.psfImages[ImageSetType::PSFIndex(i)], task.cleanCompX, task.cleanCompY, _curPeakValues[i], cleanData.startY, cleanData.endY);
		}
		
		CleanResult result;
		if(this->_cleanMask == 0)
			findPeak(*cleanData.dataImage, result.nextPeakX, result.nextPeakY, cleanData.startY, cleanData.endY);
		else
			findPeak(*cleanData.dataImage, result.nextPeakX, result.nextPeakY, cleanData.startY, cleanData.endY, this->_cleanMask);
		if(result.nextPeakX < _width)
			result.peakLevelUnnormalized = cleanData.dataImage->AbsJoinedValue(result.nextPeakX + result.nextPeakY*_width);
		else
			result.peakLevelUnnormalized = std::numeric_limits<double>::quiet_NaN();
		
		resultLane->write(result);
	}
}

template<typename ImageSetType>
std::string JoinedClean<ImageSetType>::peakDescription(const ImageSetType& image, size_t& x, size_t& y)
{
	std::ostringstream str;
	size_t index = x + y*_width;
	double peak = image.JoinedValueNormalized(index);
	str << peak << " Jy at " << x << "," << y;
	return str.str();
}

template class JoinedClean<deconvolution::SingleImageSet>;
template class JoinedClean<deconvolution::PolarizedImageSet<2>>;
template class JoinedClean<deconvolution::PolarizedImageSet<4>>;

template class JoinedClean<deconvolution::MultiImageSet<deconvolution::SingleImageSet>>;
template class JoinedClean<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<2>>>;
template class JoinedClean<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<4>>>;
