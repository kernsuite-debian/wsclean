#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../msselection.h"
#include "../multibanddata.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <memory>

class ContiguousMS final : public MSProvider
{
public:
	ContiguousMS(const string& msPath, const std::string& dataColumnName, MSSelection selection, PolarizationEnum polOut, size_t dataDescIndex, bool includeModel);
	
	ContiguousMS(const ContiguousMS&) = delete;
	
	ContiguousMS& operator=(const ContiguousMS&) = delete;
	
	virtual casacore::MeasurementSet &MS() override { return _ms; }
	
	virtual size_t RowId() const override { return _rowId; }
	
	virtual bool CurrentRowAvailable() override;
	
	virtual void NextRow() override;
	
	virtual void Reset() override;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) override;
	
	virtual void ReadData(std::complex<float>* buffer) override;
	
	virtual void ReadModel(std::complex<float>* buffer) override;
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer) override;
	
	virtual void ReadWeights(float* buffer) override;
	
	virtual void ReadWeights(std::complex<float>* buffer) override;
	
	virtual void ReopenRW() override
	{
		_ms.reopenRW();
	}
	
	virtual double StartTime() override;
	
	virtual void MakeMSRowToRowIdMapping(std::vector<size_t>& msToId, const MSSelection& selection) override;
	
	virtual void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow, const MSSelection& selection) override;
private:
	size_t _row, _rowId;
	size_t _timestep;
	double _time;
	int _dataDescId;
	bool _isMetaRead, _isDataRead, _isModelRead, _isWeightRead;
	bool _isModelColumnPrepared;
	size_t _startRow, _endRow;
	std::vector<PolarizationEnum> _inputPolarizations;
	MSSelection _selection;
	PolarizationEnum _polOut;
	std::string _msPath;
	casacore::MeasurementSet _ms;
	MultiBandData _bandData;
	bool _msHasWeights;

	casacore::ROScalarColumn<int> _antenna1Column, _antenna2Column, _fieldIdColumn, _dataDescIdColumn;
	casacore::ROScalarColumn<double> _timeColumn;
	casacore::ROArrayColumn<double> _uvwColumn;
	std::unique_ptr<casacore::ROArrayColumn<float>> _weightColumn;
	std::string _dataColumnName;
	casacore::ROArrayColumn<casacore::Complex> _dataColumn;
	casacore::ROArrayColumn<bool> _flagColumn;
	std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _modelColumn;
	
	casacore::Array<std::complex<float>> _dataArray, _modelArray;
	casacore::Array<float> _weightArray;
	casacore::Array<bool> _flagArray;
	
	void prepareModelColumn();
	void readMeta()
	{
		if(!_isMetaRead)
		{
			_dataDescId = _dataDescIdColumn(_row);
			_isMetaRead = true;
		}
	}
	void readData()
	{
		if(!_isDataRead)
		{
			_dataColumn.get(_row, _dataArray);
			_isDataRead = true;
		}
	}
	void readWeights()
	{
		if(!_isWeightRead)
		{
			_flagColumn.get(_row, _flagArray);
			if(_msHasWeights)
				_weightColumn->get(_row, _weightArray);
			_isWeightRead = true;
		}
	}
	void readModel()
	{
		if(!_isModelRead)
		{
			_modelColumn->get(_row, _modelArray);
			_isModelRead = true;
		}
	}
};

#endif
