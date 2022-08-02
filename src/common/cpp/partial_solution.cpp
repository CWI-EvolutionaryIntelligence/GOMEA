#include "common/partial_solution.hpp"

namespace gomea{

template<class T>
partial_solution_t<T>::partial_solution_t() : objective_values(vec_t<double>(1))
{
	objective_values[0] = 1e308;
}

template<class T>
partial_solution_t<T>::partial_solution_t( int num_touched_variables ) : partial_solution_t()
{
	this->touched_variables = std::vector<T>(num_touched_variables); 
}

template<class T>
partial_solution_t<T>::partial_solution_t( std::vector<T> &touched_variables, std::vector<int> &touched_indices ) : partial_solution_t()
{
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
}

template<class T>
int partial_solution_t<T>::getNumberOfTouchedVariables()
{
	return touched_variables.size();
}

template<class T>
double partial_solution_t<T>::getObjectiveValue()
{
	return objective_values[0];
}

template<class T>
double partial_solution_t<T>::getObjectiveValue( int objective_value_index )
{
	return objective_values[objective_value_index];
}

template<class T>
double partial_solution_t<T>::getConstraintValue()
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
void partial_solution_t<T>::print()
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		printf("[%d][%6.3e]",touched_indices[i],touched_variables[i]);
	//printf("\n");
}

}
