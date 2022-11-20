#pragma once

#include "gomea/src/fitness/fitness_custom.hpp"
#include "gomea/fitness.h"
#include "numpy/arrayobject.h"

namespace gomea{
namespace fitness{

template<class T>
class pyFitnessFunction_t : public customFitnessFunction_t<T>
{
	public:
		pyFitnessFunction_t( int number_of_parameters, PyObject *obj );
		pyFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj );
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );

		double mappingFunction( int objective_index, vec_t<double> &fitness_buffers );
		double mappingFunctionConstraintValue( int objective_index, vec_t<double> &fitness_buffers );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
	
	protected:
		PyObject *py_class;
		double subfunction( int subfunction_index, vec_t<T> &variables );
};

template class pyFitnessFunction_t<char>;
template class pyFitnessFunction_t<double>;

}}
