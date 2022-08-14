#pragma once

#include "fitness/fitness_custom.hpp"
#include "fitness/cython/Fitness.h"
#include "numpy/arrayobject.h"

namespace gomea{
namespace fitness{

class pyFitnessFunction_t : public customFitnessFunction_t 
{
	public:
		pyFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj );
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		PyObject *py_class;
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

}}
