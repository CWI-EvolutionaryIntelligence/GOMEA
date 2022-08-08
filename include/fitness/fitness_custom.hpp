#pragma once

#include "fitness/fitness_basic.hpp"
#include "common/solution.hpp"
#include "common/partial_solution.hpp"
#include "common/gomea_defs.hpp"
#include "utils/embed.hpp"
#include "utils/tools.hpp"

namespace gomea{
namespace fitness{

class customFitnessFunction_t : public fitness_t 
{
	public:
		customFitnessFunction_t( int number_of_parameters, double vtr );

	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		virtual double subfunction( int subfunction_index, vec_t<double> &variables ) = 0;
};

}}
