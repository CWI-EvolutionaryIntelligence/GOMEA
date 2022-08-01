#include "partial_solution.hpp"

namespace gomea{
namespace common{

template<class T>
partial_solution_t<T>::partial_solution_t( int num_touched_variables )
{
	this->touched_variables = std::vector<T>(num_touched_variables); 
}

template<class T>
partial_solution_t<T>::partial_solution_t( std::vector<T> &touched_variables, std::vector<int> &touched_indices )
{
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
}

template<class T>
partial_solution_t<T>::partial_solution_t( partial_solution_t<T> &other )
{
	for( int i = 0; i < getNumberOfTouchedVariables(); i++ )
		this->touched_indices.push_back(other.touched_indices[i]);
	this->touched_variables = other.touched_variables;
	this->objective_values = other.objective_values;
	this->constraint_value = other.constraint_value;
}

template<class T>
int partial_solution_t<T>::getNumberOfTouchedVariables()
{
	return touched_variables.size();
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

}}
