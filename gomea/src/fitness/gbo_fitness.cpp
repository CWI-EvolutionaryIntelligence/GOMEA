#include "gomea/src/fitness/gbo_fitness.hpp"

namespace gomea{
namespace fitness{

template<class T>
GBOFitnessFunction_t<T>::GBOFitnessFunction_t( int number_of_variables ) : fitness_t<T>(number_of_variables)
{
	this->name = "Custom fitness function (C++)";
}

template<class T>
GBOFitnessFunction_t<T>::GBOFitnessFunction_t( int number_of_variables, double vtr ) : fitness_t<T>(number_of_variables,vtr)
{
	this->name = "Custom fitness function (C++)";
}
		
template<class T>
void GBOFitnessFunction_t<T>::evaluationFunction( solution_t<T> *solution )
{
	solution->initFitnessBuffers(getNumberOfFitnessBuffers());
	solution->clearFitnessBuffers();
	for( int i = 0; i < this->getNumberOfSubfunctions(); i++ )
	{
		int buffer_index = this->getIndexOfFitnessBuffer(i);
		double fsub = subfunction( i, solution->variables );
		//solution->setPartialObjectiveValue(i,fsub);
		solution->addToFitnessBuffer(buffer_index, fsub);
	}

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
void GBOFitnessFunction_t<T>::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution ) //, const std::set<int> &dependent_subfunctions )
{
	solution->initFitnessBuffers(getNumberOfFitnessBuffers());
	solution->resetFitnessBuffers();

	std::set<int> touched_subfunctions;
	//if( dependent_subfunctions.size() == 0 )
	{
		for( int ind : solution->touched_indices )
		{
			assert( this->subfunction_dependency_map[ind].size() > 0 );
			touched_subfunctions.insert(this->subfunction_dependency_map[ind].begin(), this->subfunction_dependency_map[ind].end());
		}
	}
	//else
		//touched_subfunctions = dependent_subfunctions;
	
	double objective_value_delta = 0.0;
	// Calculate sum of touched subfunctions for parent
	for( int subfunction_index : touched_subfunctions )
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		objective_value_delta -= subf_result; 
		int buffer_index = this->getIndexOfFitnessBuffer(subfunction_index);
		solution->subtractFromFitnessBuffer( buffer_index, subf_result );
	}
	
	// Create backup of parent variables before modification
	vec_t<T> partial_backup = parent->getCopyOfVariables( solution->touched_indices );

	// Insert variables of partial solution and then calculate sum of touched subfunctions for offspring
	parent->insertVariables( solution->touched_variables, solution->touched_indices );
	for( int subfunction_index : touched_subfunctions ) 
	{
		double subf_result = subfunction( subfunction_index, parent->variables );
		//solution->partial_objective_values[subfunction_index] = subf_result;
		objective_value_delta += subf_result;
		int buffer_index = this->getIndexOfFitnessBuffer(subfunction_index);
		solution->addToFitnessBuffer( buffer_index, subf_result );
	}

	// Return parent variables to original state
	parent->insertVariables(partial_backup, solution->touched_indices);

	// Update fitness of partial solution
	//solution->setObjectiveValue(parent->getObjectiveValue() + objective_value_delta);
	//solution->setConstraintValue(parent->getConstraintValue());

	// Add parent buffer for final result of buffer	
	vec_t<double> parent_buffers = parent->fitness_buffers;
	for( size_t i = 0; i < parent_buffers.size(); i++ )
		solution->addToFitnessBuffer(i,parent_buffers[i]);

	// Apply function to calculate objective value from fitness buffers
	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		double ffitness = objectiveFunction( i, solution );
		solution->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(solution);
	solution->setConstraintValue(fcons);

	this->full_number_of_evaluations++;
	this->number_of_evaluations += touched_subfunctions.size() / (double) this->getNumberOfSubfunctions();
}

template<class T>
void GBOFitnessFunction_t<T>::initialize()
{
	initializeSubfunctionDependencyMap();
	initializeVariableInteractionGraph();
}

template<class T>
void GBOFitnessFunction_t<T>::initializeSubfunctionDependencyMap()
{
	for( int i = 0; i < this->number_of_variables; i++ )
	{
		this->subfunction_dependency_map[i] = std::set<int>();
	}
	for( int i = 0; i < getNumberOfSubfunctions(); i++ )
	{
		vec_t<int> dependent_variables = inputsToSubfunction(i);
		for( int j : dependent_variables )
		{
			this->subfunction_dependency_map[j].insert(i);
		}
	}
}

template<class T>
void GBOFitnessFunction_t<T>::initializeVariableInteractionGraph() 
{
	for( int i = 0; i < this->number_of_variables; i++ )
	{
		this->variable_interaction_graph[i] = std::set<int>();
	}
	for( int i = 0; i < this->getNumberOfSubfunctions(); i++ )
	{
		vec_t<int> dependent_variables = this->inputsToSubfunction(i);
		for( int vi : dependent_variables )
		{
			for (int vj : dependent_variables)
			{
				if( vi != vj )
					this->variable_interaction_graph[vi].insert(vj);
			}
		}
	}
	//this->printVariableInteractionGraph();
}

template<class T>
int GBOFitnessFunction_t<T>::getNumberOfFitnessBuffers()
{
	return 1;
}

template<class T>
int GBOFitnessFunction_t<T>::getIndexOfFitnessBuffer( int subfunction_index )
{
	return 0;
}

template<class T>
double GBOFitnessFunction_t<T>::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	if( this->variable_interaction_graph.size() == 0 )
		initializeVariableInteractionGraph();
	std::set<int> dependencies = this->variable_interaction_graph[var_a];
	if( dependencies.size() == this->number_of_variables )
		return 1.0;
	for( int v : dependencies )
	{
		if( (int) var_b == v )
			return 1.0;
	}
	return 0.0;
}

template<class T>
double GBOFitnessFunction_t<T>::objectiveFunction( int objective_index, solution_t<T> *solution )
{
	return objectiveFunction(objective_index,solution->fitness_buffers);
}

template<class T>
double GBOFitnessFunction_t<T>::objectiveFunction( int objective_index, partial_solution_t<T> *solution )
{
	return objectiveFunction(objective_index,solution->fitness_buffers);
}

template<class T>
double GBOFitnessFunction_t<T>::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
	return fitness_buffers[objective_index];
}

template<class T>
double GBOFitnessFunction_t<T>::constraintFunction( solution_t<T> *solution )
{
	return constraintFunction(solution->fitness_buffers);
}

template<class T>
double GBOFitnessFunction_t<T>::constraintFunction( partial_solution_t<T> *solution )
{
	return constraintFunction(solution->fitness_buffers);
}

template<class T>
double GBOFitnessFunction_t<T>::constraintFunction( vec_t<double> &fitness_buffers )
{
	return 0;
}

template class GBOFitnessFunction_t<char>;
template class GBOFitnessFunction_t<double>;

}}