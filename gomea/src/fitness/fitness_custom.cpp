#include "gomea/src/fitness/fitness_custom.hpp"

namespace gomea{
namespace fitness{

template<class T>
customFitnessFunction_t<T>::customFitnessFunction_t( int number_of_variables ) : fitness_t<T>(number_of_variables)
{
	this->name = "Custom fitness function (C++)";
}

template<class T>
customFitnessFunction_t<T>::customFitnessFunction_t( int number_of_variables, double vtr ) : fitness_t<T>(number_of_variables,vtr)
{
	this->name = "Custom fitness function (C++)";
}
		
template<class T>
void customFitnessFunction_t<T>::evaluationFunction( solution_t<T> *solution )
{
	double result = 0.0;
	for( int i = 0; i < this->getNumberOfSubfunctions(); i++ )
	{
		double fsub = subfunction( i, solution->variables );
		//solution->setPartialObjectiveValue(i,fsub);
		result += fsub;
	}

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	this->full_number_of_evaluations++;
	this->number_of_evaluations++;
}

template<class T>
void customFitnessFunction_t<T>::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution )
{
	std::set<int> touched_subfunctions;
	for( int ind : solution->touched_indices )
	{
		assert( this->subfunction_dependency_map[ind].size() > 0 );
		touched_subfunctions.insert(this->subfunction_dependency_map[ind].begin(), this->subfunction_dependency_map[ind].end());
	}

	double objective_value_delta = 0.0;
	// Calculate sum of touched subfunctions for parent
	for( int subfunction_index : touched_subfunctions )
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		objective_value_delta -= subf_result; 
	}
	
	// Create backup of parent variables before modification
	vec_t<T> partial_backup = parent->createPartialBackup( solution->touched_indices );

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
	this->full_number_of_evaluations++;
	this->number_of_evaluations += touched_subfunctions.size() / (double) this->getNumberOfSubfunctions();
}

}}