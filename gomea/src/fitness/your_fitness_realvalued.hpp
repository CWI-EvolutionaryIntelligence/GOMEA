#pragma once

#include "gomea/src/fitness/fitness_custom.hpp"

namespace gomea{
namespace fitness{

class yourFitnessFunctionRealValued : public customFitnessFunction_t<double> 
{
	public:
		yourFitnessFunctionRealValued( int number_of_variables, double vtr );
		
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

}}
