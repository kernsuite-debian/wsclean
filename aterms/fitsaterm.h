#ifndef FITS_ATERM_H
#define FITS_ATERM_H

#include "atermbase.h"
#include "cache.h"

#include "../fitsreader.h"

#include "../uvector.h"
#include "../windowfunction.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

/**
 * Class that reads in FITS images and resamples them onto aterm grids.
 * The fits file is supposed to have a TIME, FREQ and ANTENNA axis. 
 */
class FitsATerm : public ATermBase
{
public:
	FitsATerm(size_t nAntenna, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM, size_t atermSize);
	~FitsATerm();
	
	void OpenTECFiles(const std::vector<std::string>& filenames);
	void OpenDiagGainFiles(const std::vector<std::string>& filenames);
	
	virtual bool Calculate(std::complex<float>* buffer, double time, double frequency);
	
	virtual double AverageUpdateTime() const override;
	
	void SetTukeyWindow(double padding)
	{
		_window = WindowFunction::Tukey;
		_padding = padding;
	}
	
	void SetWindow(WindowFunction::Type window)
	{
		_window = window;
	}
	
private:
	void initializeFromFile();
	
	enum Mode { TECMode, DiagonalMode } _mode;
	
	void readImages(std::complex<float>* buffer, size_t timeIndex, double frequency);
	
	void resample(const FitsReader& reader, double* dest, const double* source);
	
	void evaluateTEC(std::complex<float>* dest, const double* source, double frequency);
	
	void copyToRealPolarization(std::complex<float>* dest, const double* source, size_t polIndex);
	void copyToImaginaryPolarization(std::complex<float>* dest, const double* source, size_t polIndex);
	void setPolarization(std::complex<float>* dest, size_t polIndex, std::complex<float> value);
	
	size_t _nAntenna, _nFrequencies, _width, _height;
	double _ra, _dec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	size_t _atermWidth, _atermHeight;
	WindowFunction::Type _window;
	double _padding;
	struct Timestep {
		double time;
		size_t readerIndex;
		size_t imgIndex;
	};
	std::vector<Timestep> _timesteps;
	Cache _cache;
	ao::uvector<double> _scratchA, _scratchB;
	size_t _curTimeindex;
	double _curFrequency;
	std::vector<FitsReader> _readers;
	std::unique_ptr<class FFTResampler> _resampler;
};

#endif
