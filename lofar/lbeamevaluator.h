#ifndef LBEAM_EVALUATOR_H
#define LBEAM_EVALUATOR_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../matrix2x2.h"

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/Station.h>
#include <StationResponse/ITRFConverter.h>
#endif

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>

#include <memory>

class LBeamEvaluator
{
public:
	class PrecalcPosInfo
	{
	public:
		friend class LBeamEvaluator;
	private:
#ifdef HAVE_LOFAR_BEAM
		LOFAR::StationResponse::vector3r_t itrfDirection;
#endif
	};
	
	explicit LBeamEvaluator(casacore::MeasurementSet& ms);
	~LBeamEvaluator();

	void Evaluate(double ra, double dec, double frequency, size_t antennaIndex, MC2x2& beamValues);
		
	void Evaluate(const PrecalcPosInfo& posInfo, double frequency, size_t antennaIndex, MC2x2& beamValues);
	
	void EvaluateFullArray(const PrecalcPosInfo& posInfo, double frequency, MC2x2& beamValues);
		
	void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, double raRad, double decRad);
	
	void SetTime(const casacore::MEpoch& time);
	
	const casacore::MEpoch& Time() const { return _time; }
	
	static void EvaluateFullCorrection(casacore::MeasurementSet& ms, double ra, double dec, const class BandData& band, MC2x2* beamValues);
	
private:
	casacore::MeasurementSet _ms;
	casacore::MEpoch _time;
	
#ifdef HAVE_LOFAR_BEAM
	std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	double _subbandFrequency, _timeAsDouble;
	casacore::MDirection _delayDir, _tileBeamDir;
	std::unique_ptr<LOFAR::StationResponse::ITRFConverter> _itrfConverter;
	LOFAR::StationResponse::vector3r_t _station0, _tile0;

	void dirToITRF(const casacore::MDirection& dir, LOFAR::StationResponse::vector3r_t& itrf)
	{
		casacore::MDirection itrfDir = _itrfConverter->toDirection(dir);
		casacore::Vector<double> itrfVal = itrfDir.getValue().getValue();
		itrf[0] = itrfVal[0];
		itrf[1] = itrfVal[1];
		itrf[2] = itrfVal[2];
	}
#endif

};

#ifndef HAVE_LOFAR_BEAM

#include <stdexcept>

// If the LOFAR beam library is not available, replace methods by dummies
inline LBeamEvaluator::LBeamEvaluator(casacore::MeasurementSet&) { }

inline LBeamEvaluator::~LBeamEvaluator() { }

inline void LBeamEvaluator::Evaluate(double, double, double, size_t, MC2x2&)
{
	throw std::runtime_error("LBeamEvaluator::Evaluate() was called, but the LOFAR beam library was not available/specified during compilation");
}
		
inline void LBeamEvaluator::Evaluate(const PrecalcPosInfo&, double, size_t, MC2x2&)
{
	throw std::runtime_error("LBeamEvaluator::Evaluate() was called, but the LOFAR beam library was not available/specified during compilation");
}
	
inline void LBeamEvaluator::EvaluateFullArray(const PrecalcPosInfo&, double, MC2x2&)
{
	throw std::runtime_error("LBeamEvaluator::Evaluate() was called, but the LOFAR beam library was not available/specified during compilation");
}
	
inline void LBeamEvaluator::PrecalculatePositionInfo(PrecalcPosInfo&, double, double) { }
	
inline void LBeamEvaluator::SetTime(const casacore::MEpoch& time) { _time = time; }

inline void LBeamEvaluator::EvaluateFullCorrection(casacore::MeasurementSet&, double, double, const class BandData&, MC2x2*)
{
	throw std::runtime_error("LBeamEvaluator::Evaluate() was called, but the LOFAR beam library was not available/specified during compilation");
}
	
#endif // !HAVE_LOFAR_BEAM

#endif
