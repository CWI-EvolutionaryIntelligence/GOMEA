#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class sphereFunction_t : public fitness_t 
{
	public:
		sphereFunction_t( int number_of_parameters, double vtr );
		
	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double x );
};

class rosenbrockFunction_t : public fitness_t 
{
	public:
		rosenbrockFunction_t( int number_of_parameters, double vtr );

		int getNumberOfSubfunctions();

		void initializeVariableInteractionGraph();
		
	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		void univariatePartialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double x, double y );
};

}}
