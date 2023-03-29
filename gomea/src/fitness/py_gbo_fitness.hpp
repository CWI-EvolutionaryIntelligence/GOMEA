#pragma once

#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/fitness.h"
#include "numpy/arrayobject.h"

namespace gomea{
namespace fitness{

template<class T>
class pyGBOFitnessFunction_t : public GBOFitnessFunction_t<T>
{
	public:
		pyGBOFitnessFunction_t( int number_of_parameters, PyObject *obj );
		pyGBOFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj );
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		double getSimilarityMeasure( size_t var_a, size_t var_b );

		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( int objective_index, vec_t<double> &fitness_buffers );
		int getNumberOfFitnessBuffers();
		int getIndexOfFitnessBuffer( int subfunction_index );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
	
	protected:
		PyObject *py_class;
		double subfunction( int subfunction_index, vec_t<T> &variables );
};

template class pyGBOFitnessFunction_t<char>;
template class pyGBOFitnessFunction_t<double>;

}}
