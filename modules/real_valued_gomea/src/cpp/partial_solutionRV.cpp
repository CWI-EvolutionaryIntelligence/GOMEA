#include "partial_solutionRV.hpp"

namespace gomea{
namespace realvalued{

template<class T>
partial_solution_t<T>::partial_solution_t( int num_touched_variables )
{
	this->touched_variables = vec_t<T>(num_touched_variables, 0.0); 
	this->sample_zs = vec_t<T>(num_touched_variables, 0.0);
}

template<class T>
partial_solution_t<T>::partial_solution_t( vec_t<T> &touched_variables, vec_t<int> &touched_indices )
{
	assert( touched_variables.size() == touched_indices.size() );
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
	this->sample_zs = vec_t<T>(getNumberOfTouchedVariables(), 0.0);
}

template<class T>
partial_solution_t<T>::partial_solution_t( vec_t<T> &touched_variables, vec_t<T> &sample_zs, vec_t<int> &touched_indices )
{
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
	this->sample_zs = sample_zs;
}
		
template<class T>
partial_solution_t<T>::partial_solution_t( partial_solution_t<T> &other )
{
	this->touched_indices = other.touched_indices;
	this->touched_variables = other.touched_variables;
	this->sample_zs = other.sample_zs;
	this->objective_value = other.objective_value;
	this->constraint_value = other.constraint_value;
}

template<class T>
int partial_solution_t<T>::getNumberOfTouchedVariables()
{
	return touched_indices.size();
}

template<class T>
double partial_solution_t<T>::getObjectiveValue()
{
	return objective_value;
}

template<class T>
double partial_solution_t<T>::getObjectiveValue( int objective_value_index )
{
	return objective_value;
}

template<class T>
double partial_solution_t<T>::getConstraintValue()
{
	return constraint_value;
}
		
template<class T>
void partial_solution_t<T>::setObjectiveValue( double v )
{
	objective_value = v;
}

template<class T>
void partial_solution_t<T>::setObjectiveValue( double v, int objective_value_index )
{
	objective_value = v;
}

template<class T>
void partial_solution_t<T>::setConstraintValue( double v )
{
	constraint_value = v;
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
void partial_solution_t<T>::setSampleMean( vec_t<T> &means )
{
	this->sample_means = means;
}

template<class T>
void partial_solution_t<T>::print()
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		printf("[%d][%6.3e]",touched_indices[i],touched_variables[i]);
	//printf("\n");
}

template class partial_solution_t<float>;
template class partial_solution_t<double>;

}}
