#pragma once

#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class sphereFunction_t : public GBOFitnessFunction_t<double>
{
	public:
		sphereFunction_t( int number_of_variables, double vtr );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

class rosenbrockFunction_t : public GBOFitnessFunction_t<double>
{
	public:
		rosenbrockFunction_t( int number_of_variables, double vtr );
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

}}
