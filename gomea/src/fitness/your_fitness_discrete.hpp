#pragma once

#include "gomea/src/fitness/gbo_fitness.hpp"

namespace gomea{
namespace fitness{

class yourFitnessFunctionDiscrete : public GBOFitnessFunction_t<char> 
{
	public:
		yourFitnessFunctionDiscrete( int number_of_variables, double vtr );
		
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );

	private:
		int trap_size = 5;
		double subfunction( int subfunction_index, vec_t<char> &variables );
};

}}
