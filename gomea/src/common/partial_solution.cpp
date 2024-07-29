#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/solution.hpp"

namespace gomea{

template<class T>
partial_solution_t<T>::partial_solution_t( int num_touched_variables )
{
	this->touched_variables = std::vector<T>(num_touched_variables); 
}

template<class T>
partial_solution_t<T>::partial_solution_t( const vec_t<T> &touched_variables, const vec_t<int> &touched_indices )
{
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
}

template<class T>
void partial_solution_t<T>::initMemory( int number_of_objectives, int number_of_fitness_buffers )
{
	initObjectiveValues(number_of_objectives);
	initFitnessBuffers(number_of_fitness_buffers);
}

template<class T>
void partial_solution_t<T>::initObjectiveValues( int number_of_objectives )
{
	if( objective_values.size() != number_of_objectives )
	{
		this->objective_values.resize(number_of_objectives);
		for (int i = 0; i < number_of_objectives; i++)
			this->objective_values[i] = INFINITY;
	}
}

template<class T>
void partial_solution_t<T>::initFitnessBuffers( int number_of_fitness_buffers )
{
	if( fitness_buffers.size() != number_of_fitness_buffers )
	{
		this->fitness_buffers.resize(number_of_fitness_buffers);
		for (int i = 0; i < number_of_fitness_buffers; i++)
			this->fitness_buffers[i] = 0.0;
	}
}

template<class T>
int partial_solution_t<T>::getNumberOfTouchedVariables()
{
	return touched_variables.size();
}

template<class T>
double partial_solution_t<T>::getObjectiveValue( int objective_value_index ) const
{
	return objective_values[objective_value_index];
}

template<class T>
const vec_t<double> partial_solution_t<T>::getObjectiveValues() const
{
	return objective_values;
}

template<class T>
double partial_solution_t<T>::getConstraintValue() const
{
	return constraint_value;
}
		
template<class T>
void partial_solution_t<T>::setObjectiveValue( double v )
{
	objective_values[0] = v;
}

template<class T>
void partial_solution_t<T>::setObjectiveValue( int objective_value_index, double v )
{
	objective_values[objective_value_index] = v;
}

template<class T>
void partial_solution_t<T>::setObjectiveValues( vec_t<double> objective_values )
{
	assert( objective_values.size() == this->objective_values.size() );
	for( size_t i = 0; i < objective_values.size(); i++ )
	{
		this->objective_values[i] = objective_values[i];
	}
}

template<class T>
void partial_solution_t<T>::setConstraintValue( double v )
{
	constraint_value = v;
}

template<class T>
double partial_solution_t<T>::getFitnessBuffer( int buffer_index ) const 
{
	return( fitness_buffers[buffer_index] );
}

template<class T>
const vec_t<double> partial_solution_t<T>::getFitnessBuffers() const
{
	return( fitness_buffers );
}
		
template<class T>
void partial_solution_t<T>::setFitnessBuffers( vec_t<double> buffers )
{
	assert( this->fitness_buffers.size() == buffers.size() );
	for( size_t i = 0; i < fitness_buffers.size(); i++ )
		fitness_buffers[i] = buffers[i];
}

template<class T>
void partial_solution_t<T>::resetFitnessBuffers()
{
	for( size_t i = 0; i < fitness_buffers.size(); i++ )
		fitness_buffers[i] = 0.0;
}
		
template<class T>
void partial_solution_t<T>::addToFitnessBuffer( int buffer_index, double partial_fitness )
{
	fitness_buffers[buffer_index] += partial_fitness;
}

template<class T>
void partial_solution_t<T>::subtractFromFitnessBuffer( int buffer_index, double partial_fitness )
{
	fitness_buffers[buffer_index] -= partial_fitness;
}

template<class T>
void partial_solution_t<T>::insertSolution( solution_t<T> *solution )
{
	for (int i = 0; i < touched_indices.size(); i++)
	{
		int ind = touched_indices[i];
		touched_variables[i] = solution->variables[ind];
	}
	setObjectiveValues( solution->getObjectiveValues() );
	setConstraintValue( solution->getConstraintValue() );
	setFitnessBuffers( solution->getFitnessBuffers() );
}

template<class T>
int partial_solution_t<T>::getTouchedIndex( int ind )
{
	if( touched_index_map.empty() )
	{
		for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
			touched_index_map[touched_indices[i]] = i;
	}

	auto map_ind = touched_index_map.find(ind);
	if( map_ind == touched_index_map.end() )
		return( -1 );
	else
		return map_ind->second;
}

template<class T>
void partial_solution_t<T>::print()
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		printf("[%d][%6.3e]",touched_indices[i],(double)touched_variables[i]);
	//printf("\n");
}

template<>
void partial_solution_t<char>::print()
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		printf("[%d][%c]",touched_indices[i],touched_variables[i]);
	//printf("\n");
}

template<>
void partial_solution_t<int>::print()
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		printf("[%d][%d]",touched_indices[i],touched_variables[i]);
	//printf("\n");
}

template<>
partial_solution_t<double>::partial_solution_t( vec_t<double> &touched_variables, vec_t<double> &sample_zs, vec_t<int> &touched_indices ) : partial_solution_t(touched_variables,touched_indices)
{
	this->sample_zs = sample_zs;
}

template<>
void partial_solution_t<double>::setSampleMean( vec_t<double> &means )
{
	this->sample_means = means;
}


template class partial_solution_t<char>;
template class partial_solution_t<double>;

}
