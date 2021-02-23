#include "wsmsgridder.h"

#include "logger.h"

#include "../imageweights.h"
#include "../buffered_lane.h"
#include "../fftresampler.h"
#include "../image.h"

#include "../msproviders/msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <fftw3.h>

#include <stdexcept>

WSMSGridder::WSMSGridder(size_t threadCount, double memFraction, double absMemLimit) :
	MSGridderBase(),
	_cpuCount(threadCount),
	_laneBufferSize(std::max<size_t>(_cpuCount*2,1024))
{
	_memSize = getAvailableMemory(memFraction, absMemLimit);
	
	// We do this once here. WStackingGridder does this too, but by default only for the float
	// variant of fftw. FFTResampler does double fft's multithreaded, hence this needs to be done
	// here too.
	fftw_make_planner_thread_safe();
}

void WSMSGridder::countSamplesPerLayer(MSData& msData)
{
	aocommon::UVector<size_t> sampleCount(ActualWGridSize(), 0);
	size_t total = 0;
	msData.matchingRows = 0;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		double uInM, vInM, wInM;
		size_t dataDescId;
		msData.msProvider->ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& bandData(msData.bandData[dataDescId]);
		for(size_t ch=msData.startChannel; ch!=msData.endChannel; ++ch)
		{
			double w = wInM / bandData.ChannelWavelength(ch);
			size_t wLayerIndex = _gridder->WToLayer(w);
			if(wLayerIndex < ActualWGridSize())
			{
				++sampleCount[wLayerIndex];
				++total;
			}
		}
		++msData.matchingRows;
		msData.msProvider->NextRow();
	}
	Logger::Debug << "Visibility count per layer: ";
	for(size_t& count : sampleCount)
	{
		Logger::Debug << count << ' ';
	}
	Logger::Debug << "\nTotal nr. of visibilities to be gridded: " << total << '\n';
}

size_t WSMSGridder::getSuggestedWGridSize() const
{
	size_t wWidth, wHeight;
	if(HasNWSize()) {
		wWidth = NWWidth(); wHeight = NWHeight();
	}
	else {
		wWidth = TrimWidth(); wHeight = TrimHeight();
	}
	double
		maxL = wWidth * PixelSizeX() * 0.5 + fabs(PhaseCentreDL()),
		maxM = wHeight * PixelSizeY() * 0.5 + fabs(PhaseCentreDM()),
		lmSq = maxL * maxL + maxM * maxM;
	double cMinW = IsComplex() ? -_maxW : _minW;
	double radiansForAllLayers;
	if(lmSq < 1.0)
		radiansForAllLayers = 2 * M_PI * (_maxW - cMinW) * (1.0 - sqrt(1.0 - lmSq));
	else
		radiansForAllLayers = 2 * M_PI * (_maxW - cMinW);
	size_t suggestedGridSize = size_t(ceil(radiansForAllLayers * NWFactor()));
	if(suggestedGridSize == 0) suggestedGridSize = 1;
	if(suggestedGridSize < _cpuCount)
	{
		// When nwlayers is lower than the nr of cores, we cannot parallellize well. 
		// However, we don't want extra w-layers if we are low on mem, as that might slow down the process
		double memoryRequired = double(_cpuCount) * double(sizeof(GridderType::num_t))*double(_actualInversionWidth*_actualInversionHeight);
		if(4.0 * memoryRequired < double(_memSize))
		{
			Logger::Info <<
				"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
				"cores (" << _cpuCount << "). Changing suggested number of w-layers to " << _cpuCount << ".\n";
			suggestedGridSize = _cpuCount;
		}
		else {
			Logger::Info <<
				"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
				"cores (" << _cpuCount << "), but there is not enough memory available to increase the number of w-layers.\n"
				"Not all cores can be used efficiently.\n";
		}
	}
	if(Verbose())
		Logger::Info << "Suggested number of w-layers: " << ceil(suggestedGridSize) << '\n';
	return suggestedGridSize;
}

void WSMSGridder::gridMeasurementSet(MSData &msData)
{
	const MultiBandData selectedBand(msData.SelectedBand());
	_gridder->PrepareBand(selectedBand);
	aocommon::UVector<std::complex<float>> modelBuffer(selectedBand.MaxChannels());
	aocommon::UVector<float> weightBuffer(selectedBand.MaxChannels());
	aocommon::UVector<bool> isSelected(selectedBand.MaxChannels());
	
	// Samples of the same w-layer are collected in a buffer
	// before they are written into the lane. This is done because writing
	// to a lane is reasonably slow; it requires holding a mutex. Without
	// these buffers, writing the lane was a bottleneck and multithreading
	// did not help. I think.
	std::vector<lane_write_buffer<InversionWorkSample>> bufferedLanes(_cpuCount);
	size_t bufferSize = std::max<size_t>(8u, _inversionCPULanes[0].capacity()/8);
	bufferSize = std::min<size_t>(128, std::min(bufferSize, _inversionCPULanes[0].capacity()));
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		bufferedLanes[i].reset(&_inversionCPULanes[i], bufferSize);
	}
	
	InversionRow newItem;
	aocommon::UVector<std::complex<float>> newItemData(selectedBand.MaxChannels());
	newItem.data = newItemData.data();
			
	size_t rowsRead = 0;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBand[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_gridder->IsInLayerRange(w1, w2))
		{
			newItem.uvw[0] = uInMeters;
			newItem.uvw[1] = vInMeters;
			newItem.uvw[2] = wInMeters;
			newItem.dataDescId = dataDescId;
			
			// Any visibilities that are not gridded in this pass
			// should not contribute to the weight sum, so set these
			// to have zero weight.
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double w = newItem.uvw[2] / curBand.ChannelWavelength(ch);
				isSelected[ch] = _gridder->IsInLayerRange(w);
			}
	
			readAndWeightVisibilities<1>(*msData.msProvider, newItem, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());
			
			InversionWorkSample sampleData;
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double wavelength = curBand.ChannelWavelength(ch);
				sampleData.sample = newItem.data[ch];
				sampleData.uInLambda = newItem.uvw[0] / wavelength;
				sampleData.vInLambda = newItem.uvw[1] / wavelength;
				sampleData.wInLambda = newItem.uvw[2] / wavelength;
				size_t cpu = _gridder->WToLayer(sampleData.wInLambda) % _cpuCount;
				bufferedLanes[cpu].write(sampleData);
			}
			
			++rowsRead;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(lane_write_buffer<InversionWorkSample>& buflane : bufferedLanes)
		buflane.write_end();
	
	if(Verbose())
		Logger::Info << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsRead;
}

void WSMSGridder::startInversionWorkThreads(size_t maxChannelCount)
{
	_inversionCPULanes.resize(_cpuCount);
	_threadGroup.clear();
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		_inversionCPULanes[i].resize(maxChannelCount * _laneBufferSize);
		set_lane_debug_name(_inversionCPULanes[i], "Work lane (buffered) containing individual visibility samples");
		_threadGroup.emplace_back(&WSMSGridder::workThreadPerSample, this, &_inversionCPULanes[i]);
	}
}

void WSMSGridder::finishInversionWorkThreads()
{
	for(std::thread& thrd : _threadGroup)
		thrd.join();
	_threadGroup.clear();
	_inversionCPULanes.clear();
}

void WSMSGridder::workThreadPerSample(aocommon::Lane<InversionWorkSample>* workLane)
{
	size_t bufferSize = std::max<size_t>(8u, workLane->capacity()/8);
	bufferSize = std::min<size_t>(128,std::min(bufferSize, workLane->capacity()));
	lane_read_buffer<InversionWorkSample> buffer(workLane, bufferSize);
	InversionWorkSample sampleData;
	while(buffer.read(sampleData))
	{
		_gridder->AddDataSample(sampleData.sample, sampleData.uInLambda, sampleData.vInLambda, sampleData.wInLambda);
	}
}

void WSMSGridder::predictMeasurementSet(MSData &msData)
{
	msData.msProvider->ReopenRW();
	const MultiBandData selectedBandData(msData.SelectedBand());
	_gridder->PrepareBand(selectedBandData);
	
	size_t rowsProcessed = 0;
	
	aocommon::Lane<PredictionWorkItem>
		calcLane(_laneBufferSize+_cpuCount),
		writeLane(_laneBufferSize);
	set_lane_debug_name(calcLane, "Prediction calculation lane (buffered) containing full row data");
	set_lane_debug_name(writeLane, "Prediction write lane containing full row data");
	lane_write_buffer<PredictionWorkItem> bufferedCalcLane(&calcLane, _laneBufferSize);
	std::thread writeThread(&WSMSGridder::predictWriteThread, this, &writeLane, &msData);
	std::vector<std::thread> calcThreads;
	for(size_t i=0; i!=_cpuCount; ++i)
		calcThreads.emplace_back(&WSMSGridder::predictCalcThread, this, &calcLane, &writeLane);
		
	/* Start by reading the u,v,ws in, so we don't need IO access
	 * from this thread during further processing */
	std::vector<double> us, vs, ws;
	std::vector<size_t> rowIds, dataIds;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBandData[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_gridder->IsInLayerRange(w1, w2))
		{
			us.push_back(uInMeters);
			vs.push_back(vInMeters);
			ws.push_back(wInMeters);
			dataIds.push_back(dataDescId);
			rowIds.push_back(msData.msProvider->RowId());
			++rowsProcessed;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(size_t i=0; i!=us.size(); ++i)
	{
		PredictionWorkItem newItem;
		newItem.u = us[i];
		newItem.v = vs[i];
		newItem.w = ws[i];
		newItem.dataDescId = dataIds[i];
		newItem.data.reset(new std::complex<float>[selectedBandData[dataIds[i]].ChannelCount()]);
		newItem.rowId = rowIds[i];
				
		bufferedCalcLane.write(std::move(newItem));
	}
	if(Verbose())
		Logger::Info << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsProcessed;
	
	bufferedCalcLane.write_end();
	for(std::thread& thr : calcThreads)
		thr.join();
	writeLane.write_end();
	writeThread.join();
}

void WSMSGridder::predictCalcThread(aocommon::Lane<PredictionWorkItem>* inputLane, aocommon::Lane<PredictionWorkItem>* outputLane)
{
	lane_write_buffer<PredictionWorkItem> writeBuffer(outputLane, _laneBufferSize);
	
	PredictionWorkItem item;
	while(inputLane->read(item))
	{
		_gridder->SampleData(item.data.get(), item.dataDescId, item.u, item.v, item.w);
		
		writeBuffer.write(std::move(item));
	}
}

void WSMSGridder::predictWriteThread(aocommon::Lane<PredictionWorkItem>* predictionWorkLane, const MSData* msData)
{
	lane_read_buffer<PredictionWorkItem> buffer(predictionWorkLane, std::min(_laneBufferSize, predictionWorkLane->capacity()));
	PredictionWorkItem workItem;
	while(buffer.read(workItem))
	{
		msData->msProvider->WriteModel(workItem.rowId, workItem.data.get());
	}
}

void WSMSGridder::Invert()
{
	std::vector<MSData> msDataVector;
	initializeMSDataVector(msDataVector);
	
	_gridder.reset(new GridderType(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(GridMode());
	if(HasDenormalPhaseCentre())
		_gridder->SetDenormalPhaseCentre(PhaseCentreDL(), PhaseCentreDM());
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == aocommon::Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(ActualWGridSize(), double(_memSize)*(6.0/10.0), _minW, _maxW);
	
	if(Verbose() && Logger::IsVerbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	resetVisibilityCounters();
	for(size_t pass=0; pass!=_gridder->NPasses(); ++pass)
	{
		Logger::Info << "Gridding pass " << pass << "... ";
		if(Verbose()) Logger::Info << '\n';
		else Logger::Info.Flush();
		
		_gridder->StartInversionPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			MSData& msData = msDataVector[i];
			
			const MultiBandData selectedBand(msData.SelectedBand());
			
			startInversionWorkThreads(selectedBand.MaxChannels());
		
			gridMeasurementSet(msData);
			
			finishInversionWorkThreads();
		}
		
		Logger::Info << "Fourier transforms...\n";
		_gridder->FinishInversionPass();
	}
	
	if(Verbose())
	{
		size_t totalRowsRead = 0, totalMatchingRows = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			totalRowsRead += msDataVector[i].totalRowsProcessed;
			totalMatchingRows += msDataVector[i].matchingRows;
		}
		
		Logger::Debug << "Total rows read: " << totalRowsRead;
		if(totalMatchingRows != 0)
			Logger::Debug << " (overhead: " << std::max(0.0, round(totalRowsRead * 100.0 / totalMatchingRows - 100.0)) << "%)";
		Logger::Debug << '\n';
	}
	
	_gridder->FinalizeImage(1.0/totalWeight(), false);
	Logger::Info << "Gridded visibility count: " << double(GriddedVisibilityCount());
	if(Weighting().IsNatural())
		Logger::Info << ", effective count after weighting: " << EffectiveGriddedVisibilityCount();
	Logger::Info << '\n';
	
	_realImage = _gridder->RealImageDouble();
	if(IsComplex())
		_imaginaryImage = _gridder->ImaginaryImageDouble();
	else
		_imaginaryImage = Image();
	
	if(ImageWidth()!=_actualInversionWidth || ImageHeight()!=_actualInversionHeight)
	{
		// Interpolate the image
		// The input is of size _actualInversionWidth x _actualInversionHeight
		FFTResampler resampler(_actualInversionWidth, _actualInversionHeight, ImageWidth(), ImageHeight(), _cpuCount);
		
		if(IsComplex())
		{
			Image resizedReal(ImageWidth(), ImageHeight());
			Image resizedImag(ImageWidth(), ImageHeight());
			resampler.Start();
			resampler.AddTask(_realImage.data(), resizedReal.data());
			resampler.AddTask(_imaginaryImage.data(), resizedImag.data());
			resampler.Finish();
			_realImage = std::move(resizedReal);
			_imaginaryImage = std::move(resizedImag);
		}
		else {
			Image resized = Image(ImageWidth(), ImageHeight());
			resampler.Resample(_realImage.data(), resized.data());
			_realImage = std::move(resized);
		}
	}
	
	if(TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight())
	{
		Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight() << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';
		// Perform trimming
		
		Image trimmed(TrimWidth(), TrimHeight());
		Image::Trim(trimmed.data(), TrimWidth(), TrimHeight(), _realImage.data(), ImageWidth(), ImageHeight());
		_realImage = std::move(trimmed);
		
		if(IsComplex())
		{
			Image trimmedImag = Image(TrimWidth(), TrimHeight());
			Image::Trim(trimmedImag.data(), TrimWidth(), TrimHeight(), _imaginaryImage.data(), ImageWidth(), ImageHeight());
			_imaginaryImage = std::move(trimmedImag);
		}
	}
	Logger::Debug << "Inversion finished.\n";
}

void WSMSGridder::Predict(Image real, Image imaginary)
{
	if(imaginary.empty() && IsComplex())
		throw std::runtime_error("Missing imaginary in complex prediction");
	if(!imaginary.empty() && !IsComplex())
		throw std::runtime_error("Imaginary specified in non-complex prediction");
	
	std::vector<MSData> msDataVector;
	initializeMSDataVector(msDataVector);
	
	_gridder = std::unique_ptr<GridderType>(new GridderType(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(GridMode());
	if(HasDenormalPhaseCentre())
		_gridder->SetDenormalPhaseCentre(PhaseCentreDL(), PhaseCentreDM());
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == aocommon::Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(ActualWGridSize(), double(_memSize)*(6.0/10.0), _minW, _maxW);
	
	if(Verbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	Image untrimmedReal, untrimmedImag;
	if(TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight())
	{
		Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight() << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
		// Undo trimming (i.e., extend with zeros)
		// The input is of size TrimWidth() x TrimHeight()
		// This will make the model image of size ImageWidth() x ImageHeight()
		untrimmedReal = Image(ImageWidth(), ImageHeight());
		Image::Untrim(untrimmedReal.data(), ImageWidth(), ImageHeight(), real.data(), TrimWidth(), TrimHeight());
		real = std::move(untrimmedReal);
		
		if(IsComplex())
		{
			untrimmedImag = Image(ImageWidth(), ImageHeight());
			Image::Untrim(untrimmedImag.data(), ImageWidth(), ImageHeight(), imaginary.data(), TrimWidth(), TrimHeight());
			imaginary = std::move(untrimmedImag);
		}
	}
	
	Image resampledReal, resampledImag;
	if(ImageWidth()!=_actualInversionWidth || ImageHeight()!=_actualInversionHeight)
	{
		// Decimate the image
		// Input is ImageWidth() x ImageHeight()
		FFTResampler resampler(ImageWidth(), ImageHeight(), _actualInversionWidth, _actualInversionHeight, _cpuCount);
		
		resampledReal = Image(ImageWidth(), ImageHeight());
		if(imaginary.empty())
		{
			resampler.Resample(real.data(), resampledReal.data());
		}
		else {
			resampledImag = Image(ImageWidth(), ImageHeight());
			resampler.Start();
			resampler.AddTask(real.data(), resampledReal.data());
			resampler.AddTask(imaginary.data(), resampledImag.data());
			resampler.Finish();
			imaginary = std::move(resampledImag);
		}
		real = std::move(resampledReal);
	}
	
	if(imaginary.empty())
		_gridder->InitializePrediction(std::move(real));
	else
		_gridder->InitializePrediction(std::move(real), std::move(imaginary));
	for(size_t pass=0; pass!=_gridder->NPasses(); ++pass)
	{
		Logger::Info << "Fourier transforms for pass " << pass << "... ";
		if(Verbose()) Logger::Info << '\n';
		else Logger::Info.Flush();
		
		_gridder->StartPredictionPass(pass);
		
		Logger::Info << "Predicting...\n";
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			predictMeasurementSet(msDataVector[i]);
	}
	
	size_t totalRowsWritten = 0, totalMatchingRows = 0;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		totalRowsWritten += msDataVector[i].totalRowsProcessed;
		totalMatchingRows += msDataVector[i].matchingRows;
	}
	
	Logger::Debug << "Total rows written: " << totalRowsWritten;
	if(totalMatchingRows != 0)
		Logger::Debug << " (overhead: " << std::max(0.0, round(totalRowsWritten * 100.0 / totalMatchingRows - 100.0)) << "%)";
	Logger::Debug << '\n';
}
