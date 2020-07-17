#ifndef ATERM_STUB
#define ATERM_STUB

#include "atermbeam.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class ATermStub : public ATermBeam
{
public:
	ATermStub(casacore::MeasurementSet&, const CoordinateSystem&, const std::string& /*dataColumnName*/)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	virtual bool calculateBeam(std::complex<float>* /*buffer*/, double /*time*/, double /*frequency*/, size_t /*fieldId*/)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	
	void SetUseDifferentialBeam(bool) { }
	void SetUseChannelFrequency(bool) { }
};

#endif

