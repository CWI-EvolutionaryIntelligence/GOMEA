#pragma once

#ifndef CPP_STANDALONE

#include "gomea/src/fitness/bbo_fitness.hpp"
#include "gomea/fitness.h"
#include "numpy/arrayobject.h"

namespace gomea{
namespace fitness{

template<class T>
class pyBBOFitnessFunction_t : public BBOFitnessFunction_t<T>
{
	public:
		pyBBOFitnessFunction_t( int number_of_parameters, PyObject *obj );
		pyBBOFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj );

		double objectiveFunction( int objective_index, vec_t<T> &variables );
		double constraintFunction( int objective_index, vec_t<T> &variables );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
	
	protected:
		PyObject *py_class;
};

}}

#endif