/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

oneMax_t::oneMax_t( int number_of_variables ) : fitness_t(number_of_variables)
{
	this->name = "OneMax function";
	this->vtr = number_of_variables;
	this->use_vtr = true;
}
		
void oneMax_t::evaluationFunction( solution_t<char> *solution )
{
	double result = 0.0;
	for( int i = 0; i < getNumberOfSubfunctions(); i++ )
		result += subfunction( solution->variables[i] );

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void oneMax_t::partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution )
{
	double result = 0.0;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		result += subfunction( solution->touched_variables[i] );
		result -= subfunction( parent->variables[ind] );
	}
	
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	
	full_number_of_evaluations++;
	number_of_evaluations += solution->getNumberOfTouchedVariables() / (double) getNumberOfSubfunctions();
}

double oneMax_t::subfunction( char x )
{
	return( x );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
