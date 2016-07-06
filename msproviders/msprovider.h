#ifndef MSPROVIDER_H
#define MSPROVIDER_H

#include "../polarizationenum.h"

#include <casacore/casa/Arrays/Array.h>

#include <complex>

namespace casacore {
	class MeasurementSet;
}
class MSSelection;

/**
 * The abstract MSProvider class is the base class for classes that read and write the visibilities.
 * An MSProvider knows which rows are selected and doesn't read or write to unselected rows. 
 * It provides the visibilities weighted with the visibility weight and converts the visibilities
 * to a requested polarization.
 * Currently, the @ref ContiguousMS and @ref PartitionedMS classes implement the MSProvider interface.
 */
class MSProvider
{
public:
	virtual ~MSProvider() { }
	
	virtual casacore::MeasurementSet &MS() = 0;
	
	virtual size_t RowId() const = 0;
	
	virtual bool CurrentRowAvailable() = 0;
	
	virtual void NextRow() = 0;
	
	virtual void Reset() = 0;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) = 0;
	
	virtual void ReadData(std::complex<float>* buffer) = 0;
	
	virtual void ReadModel(std::complex<float>* buffer) = 0;
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer) = 0;
	
	virtual void ReadWeights(float* buffer) = 0;
	
	virtual void ReadWeights(std::complex<float>* buffer) = 0;
	
	virtual void ReopenRW() = 0;
	
	virtual double StartTime() = 0;
	
	/**
	 * This function should become deprecated.
	 */
	virtual void MakeMSRowToRowIdMapping(std::vector<size_t>& msToId, const MSSelection& selection) = 0;
	
	virtual void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow, const MSSelection& selection) = 0;
	
	static std::vector<PolarizationEnum> GetMSPolarizations(casacore::MeasurementSet& ms);
protected:
	static void copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);
	
	template<typename NumType>
	static void copyWeights(NumType* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);
	
	static void reverseCopyData(casacore::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsDest, const std::complex<float>* source, PolarizationEnum polSource);
	
	static void getRowRange(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow);
	
	static void copyRealToComplex(std::complex<float>* dest, const float* source, size_t n)
	{
		const float* end = source + n;
		while(source != end)
		{
			*dest = *source;
			++dest;
			++source;
		}
	}
	
	static void initializeModelColumn(casacore::MeasurementSet& ms);
	
	MSProvider() { }
private:
	MSProvider(const MSProvider&) { }
	void operator=(const MSProvider&) { }
};

#endif
