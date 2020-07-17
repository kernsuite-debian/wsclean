// kate: space-indent off; tab-width 2; indent-width 2; replace-tabs off; eol unix;

#include "idgmsgridder.h"

#include <cmath>
#include <thread>

#include <idg-api.h>

#include "../msproviders/msprovider.h"
#include "../msproviders/timestepbuffer.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/logger.h"
#include "../wsclean/wscleansettings.h"

#include "../fitsreader.h"
#include "../fitswriter.h"

#include "../aterms/atermconfig.h"
#include "../aterms/dishaterm.h"
#include "../aterms/lofarbeamterm.h"
#include "../aterms/mwabeamterm.h"
#include "../aterms/telescope.h"

#include "idgconfiguration.h"

IdgMsGridder::IdgMsGridder(const WSCleanSettings& settings) :
	_averageBeam(nullptr),
	_outputProvider(nullptr),
	_settings(settings),
	_proxyType(idg::api::Type::CPU_OPTIMIZED),
	_buffersize(0)
{
	IdgConfiguration::Read(_proxyType, _buffersize, _options);
	setIdgType();
	_bufferset = std::unique_ptr<idg::api::BufferSet>(
		idg::api::BufferSet::create(_proxyType));
	if(settings.gridWithBeam || !settings.atermConfigFilename.empty())
		_options["a_term_kernel_size"] = float(_settings.atermKernelSize);
	_options["max_threads"] = int(settings.threadCount);
	if(settings.gridMode == BlackmanHarrisKernel)
		_options["taper"] = std::string("blackman-harris");
}

void IdgMsGridder::Invert()
{
	const size_t untrimmedWidth = ImageWidth();
	const size_t width = TrimWidth(), height = TrimHeight();

	assert(width == height);
	assert(untrimmedWidth == ImageHeight());

	_options["padded_size"] = untrimmedWidth;

	// Stokes I is always the first requested pol. So, only when Stokes I
	// is requested, do the actual inversion. Since all pols are produced at once by IDG,
	// when Stokes Q/U/V is requested, the result of earlier gridding is returned.
	// For the first inversion, WSClean will ask for the PSF of Stokes I. The next run
	// is the dirty of Stokes I, which should thus overwrite all images.
	if(Polarization() == aocommon::Polarization::StokesI)
	{
		if (!_metaDataCache->averageBeam) _metaDataCache->averageBeam.reset(new AverageBeam());
		_averageBeam = static_cast<AverageBeam*>(_metaDataCache->averageBeam.get());

		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector);
		
		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxWWithFlags);
		}

		double shiftl = 0.0, shiftm = 0.0, shiftp = 0.0; // TODO
		_bufferset->init(width, _actualPixelSizeX, max_w+1.0, shiftl, shiftm, shiftp, _options);
		Logger::Debug << "IDG subgrid size: " << _bufferset->get_subgridsize() << '\n';

		if (DoImagePSF())
		{
			// Computing the PSF
			// For the PSF the aterm is computed but not applied
			// The aterm is computed so that the average beam can be computed
			_bufferset->set_apply_aterm(false);
			_bufferset->unset_matrix_inverse_beam();
			_bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_and_grid);
			resetVisibilityCounters();
			for(size_t i=0; i!=MeasurementSetCount(); ++i)
			{
				// Adds the gridding result to _image member
				gridMeasurementSet(msDataVector[i]);
			}
			_bufferset->finalize_compute_avg_beam();
			_averageBeam->SetScalarBeam(_bufferset->get_scalar_beam());
			_averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam());
			_image.assign(4 * width * height, 0.0);
			_bufferset->get_image(_image.data());

			Logger::Debug << "Total weight: " << totalWeight() << '\n';

			double center_pixel_value = _image[height/2 * width + width/2]; // TODO check memory layout, is this correct? for now it does not matter, because width == height

			if (center_pixel_value)
			{
				for(size_t ii=0; ii != 4 * width * height; ++ii)
				{
					_image[ii] /= center_pixel_value;
				}
			}

		}
		else {
			// Compute a dirty/residual image
			// with application of the a term
			_bufferset->set_apply_aterm(true);

			// Because compensation for the average beam happens at subgrid level
			// it needs to be known in advance.
			// If it is not in the cache it needs to be computed first
			if (!_averageBeam->Empty())
			{
				// Set avg beam from cache
				Logger::Debug << "Using average beam from cache.\n";
				_bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
				_bufferset->set_matrix_inverse_beam(_averageBeam->MatrixInverseBeam());
			}
			else {
				// Compute avg beam
				Logger::Debug << "Computing average beam.\n";
				_bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_only);
				for(size_t i=0; i!=MeasurementSetCount(); ++i)
				{
					gridMeasurementSet(msDataVector[i]);
				}
				_bufferset->finalize_compute_avg_beam();
				Logger::Debug << "Finished computing average beam.\n";
				_averageBeam->SetScalarBeam(_bufferset->get_scalar_beam());
				_averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam());
			}

			resetVisibilityCounters();
			for(size_t i=0; i!=MeasurementSetCount(); ++i)
			{
				// Adds the gridding result to _image member
				gridMeasurementSet(msDataVector[i]);
			}
			_image.assign(4 * width * height, 0.0);
			_bufferset->get_image(_image.data());
		}

		// result is now in _image member
		// Can be accessed by subsequent calls to ImageRealResult()
		
	}
	else if(_image.empty()) {
		throw std::runtime_error("IdgMsGridder::Invert() was called out of sequence");
	}
}

std::unique_ptr<class ATermBase> IdgMsGridder::getATermMaker(MSGridderBase::MSData& msData)
{
	SynchronizedMS ms = msData.msProvider->MS();
	size_t nr_stations = ms->antenna().nrow();
	std::unique_ptr<ATermBase> aTermMaker;
	aocommon::UVector<std::complex<float>> aTermBuffer;
	if(!_settings.atermConfigFilename.empty() || _settings.gridWithBeam)
	{
		ATermBase::CoordinateSystem system;
		// IDG uses a flipped coordinate system which is moved by half a pixel:
		system.dl = -_bufferset->get_subgrid_pixelsize();
		system.dm = -_bufferset->get_subgrid_pixelsize();
		system.phaseCentreDL = PhaseCentreDL() - 0.5*system.dl;
		system.phaseCentreDM = PhaseCentreDM() + 0.5*system.dm;
		system.width = _bufferset->get_subgridsize();
		system.height = system.width;
		system.maxSupport = _settings.atermKernelSize;
		system.ra = PhaseCentreRA();
		system.dec = PhaseCentreDec();
		if(!_settings.atermConfigFilename.empty())
		{
			std::unique_ptr<ATermConfig> config(new ATermConfig(nr_stations, system, _settings));
			config->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
			config->Read(ms.Filename(), *ms, _settings.atermConfigFilename);
			return std::move(config);
		}
		else {
			switch(Telescope::GetType(*ms))
			{
				case Telescope::AARTFAAC:
				case Telescope::LOFAR: {
					std::unique_ptr<LofarBeamTerm> beam(new LofarBeamTerm(*ms, system, _settings.dataColumnName));
					beam->SetUseDifferentialBeam(_settings.useDifferentialLofarBeam);
					beam->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
					beam->SetUpdateInterval(_settings.beamAtermUpdateTime);
					return std::move(beam);
				}
				case Telescope::MWA: {
					std::unique_ptr<MWABeamTerm> beam(new MWABeamTerm(*ms, system));
					beam->SetUpdateInterval(_settings.beamAtermUpdateTime);
					beam->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
					beam->SetSearchPath(_settings.mwaPath);
					return std::move(beam);
				}
				case Telescope::VLA: {
					std::unique_ptr<DishATerm> beam(new DishATerm(*ms, system));
					beam->SetUpdateInterval(_settings.beamAtermUpdateTime);
					beam->SetSaveATerms(_settings.saveATerms, _settings.prefixName);
					return std::move(beam);
				}
				default:
					throw std::runtime_error("Can't make beam for this telescope");
			}
		}
	}
	else {
		return std::unique_ptr<ATermBase>();
	}
}

void IdgMsGridder::gridMeasurementSet(MSGridderBase::MSData& msData)
{
	std::unique_ptr<ATermBase> aTermMaker;
	aocommon::UVector<std::complex<float>> aTermBuffer;
	if(!prepareForMeasurementSet(msData, aTermMaker, aTermBuffer, idg::api::BufferSetType::gridding))
		return;
	
	aocommon::UVector<float> weightBuffer(_selectedBands.MaxChannels()*4);
	aocommon::UVector<std::complex<float>> modelBuffer(_selectedBands.MaxChannels()*4);
	aocommon::UVector<bool> isSelected(_selectedBands.MaxChannels()*4, true);
	aocommon::UVector<std::complex<float>> dataBuffer(_selectedBands.MaxChannels()*4);
	// The gridder doesn't need to know the absolute time index; this value indexes relatively to where we
	// start in the measurement set, and only increases when the time changes.
	int timeIndex = -1;
	double currentTime = -1.0;
	aocommon::UVector<double> uvws(msData.msProvider->NAntennas()*3, 0.0);
	for(TimestepBuffer timestepBuffer(msData.msProvider, DoSubtractModel()) ; timestepBuffer.CurrentRowAvailable() ; timestepBuffer.NextRow())
	{
		MSProvider::MetaData metaData;
		timestepBuffer.ReadMeta(metaData);

		if(currentTime != metaData.time)
		{
			currentTime = metaData.time;
			timeIndex++;
			
			if(aTermMaker)
			{
				timestepBuffer.GetUVWsForTimestep(uvws);
				if(aTermMaker->Calculate(aTermBuffer.data(), currentTime, _selectedBands[metaData.dataDescId].CentreFrequency(), metaData.fieldId, uvws.data()))
				{
					_bufferset->get_gridder(metaData.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
					Logger::Debug << "Calculated a-terms for timestep " << timeIndex << "\n";
				}
			}
		}
		const BandData& curBand(_selectedBands[metaData.dataDescId]);
		IDGInversionRow rowData;
		
		rowData.data = dataBuffer.data();
		rowData.uvw[0] = metaData.uInM;
		rowData.uvw[1] = metaData.vInM;
		rowData.uvw[2] = metaData.wInM;
		
		rowData.antenna1 = metaData.antenna1;
		rowData.antenna2 = metaData.antenna2;
		rowData.timeIndex = timeIndex;
		rowData.dataDescId = metaData.dataDescId;
		readAndWeightVisibilities<4>(timestepBuffer, rowData, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());

		rowData.uvw[1] = -metaData.vInM;  // DEBUG vdtol, flip axis
		rowData.uvw[2] = -metaData.wInM;  //

		_bufferset->get_gridder(rowData.dataDescId)->grid_visibilities(timeIndex, metaData.antenna1, metaData.antenna2, rowData.uvw, rowData.data, weightBuffer.data());
	}
	
	_bufferset->finished();
}

void IdgMsGridder::Predict(Image image)
{
	const size_t untrimmedWidth = ImageWidth();
	const size_t width = TrimWidth(), height = TrimHeight();

	assert(width == height);
	assert(untrimmedWidth == ImageHeight());

	_options["padded_size"] = untrimmedWidth;

	if (Polarization() == aocommon::Polarization::StokesI)
	{
		_image.assign(4 * width * height, 0.0);
		if (!_metaDataCache->averageBeam)
		{
			Logger::Info << "no average_beam in cache, creating an empty one.\n";
			_metaDataCache->averageBeam.reset(new AverageBeam());
		}
		_averageBeam = static_cast<AverageBeam*>(_metaDataCache->averageBeam.get());
	}

	size_t polIndex = aocommon::Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_image[i + polIndex*width*height] = image[i];

	// Stokes V is the last requested pol, unless only Stokes I is imaged. Only when the last
	// polarization is given, do the actual prediction.
	if (Polarization() == aocommon::Polarization::StokesV || (Polarization() == aocommon::Polarization::StokesI && _settings.polarizations.size()==1))
	{
		// Do actual predict

		bool do_scale = false;
		if (!_averageBeam->Empty())
		{
			// Set avg beam from cache
			Logger::Debug << "Average beam is already in cache.\n";
			_bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
			do_scale = true;
		}

		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector);

		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxWWithFlags);
		}

		float shiftl = 0.0, shiftm = 0.0, shiftp = 0.0; // TODO
		_bufferset->init(width, _actualPixelSizeX, max_w+1.0, shiftl, shiftm, shiftp, _options);
		_bufferset->set_image(_image.data(), do_scale);

		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			predictMeasurementSet(msDataVector[i]);
		}
	}
}

void IdgMsGridder::setIdgType()
{
	switch(_settings.idgMode)
	{
		default:
			return;
		case WSCleanSettings::IDG_CPU:
			_proxyType = idg::api::Type::CPU_OPTIMIZED;
			return;
		case WSCleanSettings::IDG_GPU:
			_proxyType = idg::api::Type::CUDA_GENERIC;
			return;
		case WSCleanSettings::IDG_HYBRID:
			_proxyType = idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
			return;
	}
}

bool IdgMsGridder::prepareForMeasurementSet(MSGridderBase::MSData& msData, std::unique_ptr<ATermBase>& aTermMaker, aocommon::UVector<std::complex<float>>& aTermBuffer, idg::api::BufferSetType bufferSetType)
{
	const float max_baseline = msData.maxBaselineInM;
	// Skip this ms if there is no data in it
	if (!max_baseline) return false;
	
	_selectedBands = msData.SelectedBand();

	// TODO for now we map the ms antennas directly to the gridder's antenna,
	// including non-selected antennas. Later this can be made more efficient.
	size_t nStations = msData.msProvider->MS()->antenna().nrow();
	
	std::vector<std::vector<double>> bands;
	for(size_t i=0; i!=_selectedBands.BandCount(); ++i)
	{
		bands.push_back(std::vector<double>(_selectedBands[i].begin(), _selectedBands[i].end()));
	}
	size_t nChannels = _selectedBands.MaxChannels();

	aTermMaker = getATermMaker(msData);
	size_t subgridsize = _bufferset->get_subgridsize();
	if(aTermMaker)
		aTermBuffer.resize(subgridsize*subgridsize*4*nStations);
	
	uint64_t memSize = getAvailableMemory(_settings.memFraction, _settings.absMemLimit);
	uint64_t memPerTimestep = idg::api::BufferSet::get_memory_per_timestep(nStations, nChannels);
	if(aTermMaker)
	{
		// When a-terms are used, they will also take memory. Here we calculate their approx contribution.
		double avgUpdate = aTermMaker->AverageUpdateTime();
		Logger::Debug << "A-terms change on average every " << avgUpdate << " s, once every " << (avgUpdate / msData.integrationTime) << " timesteps.\n";
		uint64_t atermMemPerTimestep =
			subgridsize*subgridsize*nStations *   // size of grid x nr of grids
			(4*8) *                               // 4 pol, 8 bytes per complex value
			(msData.integrationTime / avgUpdate); // Average number of aterms per timestep
		Logger::Debug << "A-terms increase mem per timestep from " << memPerTimestep << " bytes to " << (memPerTimestep + atermMemPerTimestep) << " bytes.\n";
		memPerTimestep += atermMemPerTimestep;
	}
	
	// IDG can allocate two visibility buffers: (for parallel processing)
	memPerTimestep *= 2;
	
	// Only one-third of the mem is allocated to the buffers, so that memory remains available for the images
	// and other things done by IDG.
	_buffersize = std::max<size_t>(1, memSize/3 / memPerTimestep);
	
	Logger::Debug << "Allocatable timesteps (" << nStations << " stations, " << nChannels << " channels, " <<
		memSize/(1024*1024*1024) << " GB mem): " << _buffersize << '\n';
	_bufferset->init_buffers(_buffersize, bands, nStations, max_baseline, _options, bufferSetType);
	
	return true;
}

void IdgMsGridder::predictMeasurementSet(MSGridderBase::MSData& msData)
{
	std::unique_ptr<ATermBase> aTermMaker;
	aocommon::UVector<std::complex<float>> aTermBuffer;
	if(!prepareForMeasurementSet(msData, aTermMaker, aTermBuffer, idg::api::BufferSetType::degridding))
		return;
	
	msData.msProvider->ReopenRW();

	_outputProvider = msData.msProvider;
	
	aocommon::UVector<std::complex<float>> buffer(_selectedBands.MaxChannels()*4);
	int timeIndex = -1;
	double currentTime = -1.0;
	aocommon::UVector<double> uvws(msData.msProvider->NAntennas()*3, 0.0);
	for(TimestepBuffer timestepBuffer(msData.msProvider, false) ; timestepBuffer.CurrentRowAvailable() ; timestepBuffer.NextRow())
	{
		MSProvider::MetaData metaData;
		timestepBuffer.ReadMeta(metaData);

		size_t provRowId = timestepBuffer.RowId();
		if(currentTime != metaData.time)
		{
			currentTime = metaData.time;
			timeIndex++;
			
			if(aTermMaker)
			{
				timestepBuffer.GetUVWsForTimestep(uvws);
				if(aTermMaker->Calculate(aTermBuffer.data(), currentTime, _selectedBands[metaData.dataDescId].CentreFrequency(), metaData.fieldId, uvws.data()))
				{
					_bufferset->get_degridder(metaData.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
					Logger::Debug << "Calculated new a-terms for timestep " << timeIndex << "\n";
				}
			}
		}
		
		IDGPredictionRow row;
		row.uvw[0] = metaData.uInM;
		row.uvw[1] = -metaData.vInM;
		row.uvw[2] = -metaData.wInM;
		row.antenna1 = metaData.antenna1;
		row.antenna2 = metaData.antenna2;
		row.timeIndex = timeIndex;
		row.dataDescId = metaData.dataDescId;
		row.rowId = provRowId;
		predictRow(row);
  }
	
	for(size_t d=0; d!=_selectedBands.DataDescCount(); ++d)
		computePredictionBuffer(d);
}

void IdgMsGridder::predictRow(IDGPredictionRow& row)
{
	while (_bufferset->get_degridder(row.dataDescId)->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw))
	{
		computePredictionBuffer(row.dataDescId);
	}
}

void IdgMsGridder::computePredictionBuffer(size_t dataDescId)
{
	auto available_row_ids = _bufferset->get_degridder(dataDescId)->compute();
	Logger::Debug << "Computed " << available_row_ids.size() << " rows.\n";
	for(auto i : available_row_ids)
	{
		_outputProvider->WriteModel(i.first, i.second);
	}
	_bufferset->get_degridder(dataDescId)->finished_reading();
}

void IdgMsGridder::Predict(Image /*real*/, Image /*imaginary*/)
{
	throw std::runtime_error("IDG gridder cannot make complex images");
}

Image IdgMsGridder::ImageRealResult()
{
	const size_t width = TrimWidth(), height = TrimHeight();
	size_t polIndex = aocommon::Polarization::StokesToIndex(Polarization());
	Image image(height, width);
	std::copy_n(_image.data()+polIndex*width*height, width*height, image.data());
	return std::move(image);
}

Image IdgMsGridder::ImageImaginaryResult()
{
	throw std::runtime_error("IDG gridder cannot make complex images");
}

void IdgMsGridder::SaveBeamImage(const ImagingTableEntry& entry, ImageFilename& filename, const WSCleanSettings& settings, double ra, double dec, double pdl, double pdm, const MetaDataCache& cache)
{
	if (!cache.averageBeam || cache.averageBeam->Empty())
	{
		throw std::runtime_error("IDG gridder can not save the beam image. Beam has not been computed yet.");
	}
	FitsWriter writer;
	writer.SetImageDimensions(settings.trimmedImageWidth, settings.trimmedImageHeight, ra, dec, settings.pixelScaleX, settings.pixelScaleY);
	writer.SetPhaseCentreShift(pdl, pdm);
	ImageFilename polName(filename);
	polName.SetPolarization(aocommon::Polarization::StokesI);
	writer.SetPolarization(aocommon::Polarization::StokesI);
	writer.SetFrequency(entry.CentralFrequency(), entry.bandEndFrequency - entry.bandStartFrequency);
	writer.Write(polName.GetBeamPrefix(settings) + ".fits", cache.averageBeam->ScalarBeam()->data());
}

void IdgMsGridder::SavePBCorrectedImages(FitsWriter& writer, const ImageFilename& filename, const std::string& filenameKind, const WSCleanSettings& settings)
{
	ImageFilename beamName(filename);
	beamName.SetPolarization(aocommon::Polarization::StokesI);
	FitsReader reader(beamName.GetBeamPrefix(settings) + ".fits");
	
	Image beam(reader.ImageWidth(), reader.ImageHeight());
	reader.Read(beam.data());
	
	Image image;
	for(size_t polIndex = 0; polIndex != 4; ++polIndex)
	{
		aocommon::PolarizationEnum pol = aocommon::Polarization::IndexToStokes(polIndex);
		ImageFilename name(filename);
		name.SetPolarization(pol);
		FitsReader reader(name.GetPrefix(settings) + "-" + filenameKind + ".fits");
		if(image.empty())
			image = Image(reader.ImageWidth(), reader.ImageHeight());
		reader.Read(image.data());
		
		for(size_t i=0; i!=reader.ImageWidth() * reader.ImageHeight(); ++i)
		{
			if(beam[i] > 1e-6)
				image[i] /= beam[i];
			else
				image[i] = std::numeric_limits<double>::quiet_NaN();
		}
		
		writer.SetPolarization(pol);
		writer.Write(name.GetPrefix(settings) + "-" + filenameKind + "-pb.fits", image.data());
	}
}
