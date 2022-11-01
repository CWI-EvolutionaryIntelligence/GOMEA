#include "gomea/src/fitness/fitness_custom.hpp"

namespace gomea{
namespace fitness{

customFitnessFunction_t::customFitnessFunction_t( int number_of_parameters, double vtr ) : fitness_t(number_of_parameters,vtr)
{
	this->name = "Custom fitness function (C++)";
}
		
void customFitnessFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < getNumberOfSubfunctions(); i++ )
	{
		double fsub = subfunction( i, solution->variables );
		//solution->setPartialObjectiveValue(i,fsub);
		result += fsub;
	}

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void customFitnessFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	std::set<int> touched_subfunctions;
	for( int ind : solution->touched_indices )
	{
		assert( subfunction_dependency_map[ind].size() > 0 );
		touched_subfunctions.insert(subfunction_dependency_map[ind].begin(), subfunction_dependency_map[ind].end());
	}


	double objective_value_delta = 0.0;
	// Calculate sum of touched subfunctions for parent
	for( int subfunction_index : touched_subfunctions )
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		objective_value_delta -= subf_result; 
	}
	
	// Create backup of parent variables before modification
	vec_t<double> partial_backup = parent->createPartialBackup( solution->touched_indices );

	// Insert variables of partial solution and then calculate sum of touched subfunctions for offspring
	parent->insertVariables( solution->touched_variables, solution->touched_indices );
	for( int subfunction_index : touched_subfunctions ) 
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		//solution->partial_objective_values[subfunction_index] = subf_result;
		objective_value_delta += subf_result;
	}

	// Return parent to original state
	parent->insertVariables(partial_backup, solution->touched_indices);
	
	solution->setObjectiveValue(parent->getObjectiveValue() + objective_value_delta);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += touched_subfunctions.size() / (double) getNumberOfSubfunctions();
}

}}