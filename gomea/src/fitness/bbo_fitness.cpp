#include "gomea/src/fitness/bbo_fitness.hpp"

namespace gomea{
namespace fitness{

template<class T>
BBOFitnessFunction_t<T>::BBOFitnessFunction_t( int number_of_variables ) : fitness_t<T>(number_of_variables)
{
	this->name = "Custom fitness function (C++)";
}

template<class T>
BBOFitnessFunction_t<T>::BBOFitnessFunction_t( int number_of_variables, double vtr ) : fitness_t<T>(number_of_variables,vtr)
{
	this->name = "Custom fitness function (C++)";
}

template<class T>
void BBOFitnessFunction_t<T>::initialize() 
{
	return;
}
		
template<class T>
void BBOFitnessFunction_t<T>::evaluationFunction( solution_t<T> *solution )
{
	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		double ffitness = objectiveFunction( i, solution );
		solution->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(solution);
	solution->setConstraintValue(fcons);

	this->full_number_of_evaluations++;
	this->number_of_evaluations++;
}

template<class T>
void BBOFitnessFunction_t<T>::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution )
{
	// Create backup of parent variables before modification
	vec_t<T> partial_backup = parent->getCopyOfVariables( solution->touched_indices );

	// Insert variables of partial solution and then calculate fitness
	parent->insertVariables( solution->touched_variables, solution->touched_indices );
	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		double ffitness = objectiveFunction( i, parent );
		solution->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(parent);
	solution->setConstraintValue(fcons);

	// Return parent variables to original state
	parent->insertVariables(partial_backup, solution->touched_indices);

	this->full_number_of_evaluations++;
	this->number_of_evaluations++;
}

template<class T>
double BBOFitnessFunction_t<T>::objectiveFunction( int objective_index, solution_t<T> *solution )
{
	return objectiveFunction(objective_index,solution->variables);
}

template<class T>
double BBOFitnessFunction_t<T>::constraintFunction( solution_t<T> *solution )
{
	return constraintFunction(solution->variables);
}

template<class T>
double BBOFitnessFunction_t<T>::constraintFunction( vec_t<T> &variables )
{
	return 0;
}

template class BBOFitnessFunction_t<char>;
template class BBOFitnessFunction_t<double>;

}}
