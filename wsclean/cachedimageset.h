#ifndef CACHED_IMAGE_SET_H
#define CACHED_IMAGE_SET_H

#include "../fitswriter.h"
#include "../fitsreader.h"

#include "../wsclean/logger.h"

#include <string.h>
#include <set>

class CachedImageSet
{
public:
	CachedImageSet() : _polCount(0), _freqCount(0), _image()
	{ }
	
	~CachedImageSet()
	{
		for(const std::string& filename : _storedNames)
			std::remove(filename.c_str());
	}
	
	CachedImageSet(const CachedImageSet& source) = delete;
	CachedImageSet& operator=(const CachedImageSet& source) = delete;
	
	void Initialize(const FitsWriter& writer, size_t polCount, size_t freqCount, const std::string& prefix)
	{
		_writer = writer;
		_polCount = polCount;
		_freqCount = freqCount;
		_prefix = prefix;
		_image.reset();
	}
	
	void SetFitsWriter(const FitsWriter& writer)
	{
		_writer = writer;
	}
	
	void Load(double* image, aocommon::PolarizationEnum polarization, size_t freqIndex, bool isImaginary) const
	{
		if(_writer.Width() == 0 || _writer.Height() == 0)
			throw std::runtime_error("Writer is not set.");
		Logger::Debug << "Loading " << name(polarization, freqIndex, isImaginary) << '\n';
		if(_polCount == 1 && _freqCount == 1)
			if(_image.empty())
				throw std::runtime_error("Loading image before store");
			else
				std::copy(_image.data(), _image.data() + _writer.Width()*_writer.Height(), image);
		else {
			FitsReader reader(name(polarization, freqIndex, isImaginary));
			reader.Read(image);
		}
	}
	
	void Store(const double* image, aocommon::PolarizationEnum polarization, size_t freqIndex, bool isImaginary)
	{
		if(_writer.Width() == 0 || _writer.Height() == 0)
			throw std::runtime_error("Writer is not set.");
		Logger::Debug << "Storing " << name(polarization, freqIndex, isImaginary) << '\n';
		if(_polCount == 1 && _freqCount == 1)
		{
			if(_image.empty())
			{
				_image = Image(_writer.Width(), _writer.Height());
			}
			std::copy(image, image + _writer.Width() * _writer.Height(), _image.data());
			//memcpy(_image, image, _writer.Width() * _writer.Height() * sizeof(double));
		}
		else {
			std::string n = name(polarization, freqIndex, isImaginary);
			_writer.Write(n, image);
			_storedNames.insert(n);
		}
	}
	
private:
	std::string name(aocommon::PolarizationEnum polarization, size_t freqIndex, bool isImaginary) const
	{
		if(_freqCount == 1)
		{
			if(isImaginary)
				return _prefix + '-' + aocommon::Polarization::TypeToShortString(polarization) + "i-tmp.fits";
			else
				return _prefix + '-' + aocommon::Polarization::TypeToShortString(polarization) + "-tmp.fits";
		}
		else {
			std::ostringstream str;
			str <<  _prefix + '-' + aocommon::Polarization::TypeToShortString(polarization);
			if(isImaginary)
				str << 'i';
			str << '-';
			if(freqIndex < 10) str << '0';
			if(freqIndex < 100) str << '0';
			if(freqIndex < 1000) str << '0';
			str << freqIndex << "-tmp.fits";
			return str.str();
		}
	}
	FitsWriter _writer;
	size_t _polCount, _freqCount;
	std::string _prefix;
	
	Image _image;
	std::set<std::string> _storedNames;
};

#endif
