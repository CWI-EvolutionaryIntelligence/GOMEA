#include "gomea/src/common/solution.hpp"

namespace gomea{

template<class T>
solution_t<T>::solution_t( int number_of_variables )
{
	this->variables = vec_t<T>(number_of_variables);
 	this->objective_values = vec_t<double>(1);
	objective_values[0] = INFINITY;
}

template<class T>
solution_t<T>::solution_t( vec_t<T> &variables )
{
	this->variables = vec_t<T>(variables);
 	this->objective_values = vec_t<double>(1);
	objective_values[0] = INFINITY;
}

template<class T>
int solution_t<T>::getNumberOfVariables()
{
	return variables.size();
}

template<class T>
int solution_t<T>::getNumberOfObjectives()
{
	return objective_values.size();
}

template<class T>
double solution_t<T>::getObjectiveValue()
{
	return objective_values[0];
}

template<class T>
double solution_t<T>::getObjectiveValue( int objective_value_index )
{
	return objective_values[objective_value_index];
}

template<class T>
double solution_t<T>::getConstraintValue()
{
	return constraint_value;
}

template<class T>
double solution_t<T>::getPartialObjectiveValue( int subfunction_index )
{
	return partial_objective_values[subfunction_index];
}

template<class T>
double solution_t<T>::getPartialConstraintValue( int subfunction_index )
{
	return partial_constraint_values[subfunction_index];
}
		
template<class T>
void solution_t<T>::setObjectiveValue( double v )
{
	objective_values[0] = v;
}

template<class T>
void solution_t<T>::setObjectiveValue( int objective_value_index, double v )
{
	objective_values[objective_value_index] = v;
}

template<class T>
void solution_t<T>::setConstraintValue( double v )
{
	constraint_value = v;
}

template<class T>
void solution_t<T>::setPartialObjectiveValue( int subfunction_index, double v )
{
	partial_objective_values[subfunction_index] = v;
}

template<class T>
void solution_t<T>::setPartialConstraintValue( int subfunction_index, double v )
{
	partial_constraint_values[subfunction_index] = v;
}

template<class T>
vec_t<T> solution_t<T>::createPartialBackup( vec_t<int> variable_indices )
{
	vec_t<T> backup = vec_t<T>(variable_indices.size());
	for( int i = 0; i < variable_indices.size(); i++ )
	{
		int ind = variable_indices[i];
		backup[i] = variables[ind];
	}
	return backup;
}

template<class T>
void solution_t<T>::insertVariables( vec_t<T> vars_to_insert, vec_t<int> indices_to_insert )
{
	for( int i = 0; i < indices_to_insert.size(); i++ )
	{
		int ind = indices_to_insert[i];
		variables[ind] = vars_to_insert[i];
	}
}

template<class T>
void solution_t<T>::print()
{
	for( int i = 0; i < variables.size(); i++ )
		printf("%6.3e ",variables[i]);
	printf("\n");
}

}
