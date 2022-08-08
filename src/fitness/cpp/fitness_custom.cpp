#include "fitness/fitness_custom.hpp"

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
		solution->setPartialObjectiveValue(i,fsub);
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
	for( int ind : solution->touched_variables )
	{
		touched_subfunctions.insert(subfunction_dependency_graph[ind].begin(), subfunction_dependency_graph[ind].end());
	}

	vec_t<double> partial_backup = parent->createPartialBackup( solution->touched_indices );
	
	double objective_value_delta = 0.0;
	for( int subfunction_index : touched_subfunctions ) 
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		solution->partial_objective_values[subfunction_index] = subf_result;
		objective_value_delta += subf_result - parent->getPartialObjectiveValue(subfunction_index);
	}

	parent->insertPartialBackup(partial_backup,solution->touched_indices);
	
	solution->setObjectiveValue(parent->getObjectiveValue() + objective_value_delta);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += solution->getNumberOfTouchedVariables() / (double) getNumberOfSubfunctions();
}

double customFitnessFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	return( variables[subfunction_index] * variables[subfunction_index] );
}

}}