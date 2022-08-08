#pragma once

#include "fitness/fitness_custom.hpp"

namespace gomea{
namespace fitness{

class yourFitnessFunction_t : public customFitnessFunction_t 
{
	public:
		yourFitnessFunction_t( int number_of_parameters, int number_of_subfunctions, double vtr );
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		virtual double subfunction( int subfunction_index, vec_t<double> &variables );
};

}}
