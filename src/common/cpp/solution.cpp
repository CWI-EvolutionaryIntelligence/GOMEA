#include "common/solution.hpp"

namespace gomea{

template<class T>
solution_t<T>::solution_t() : objective_values(vec_t<double>(1))
{
}

template<class T>
solution_t<T>::solution_t( int number_of_variables ) : solution_t()
{
	this->variables = vec_t<T>(number_of_variables);
}

template<class T>
solution_t<T>::solution_t( vec_t<T> &variables ) : solution_t()
{
	this->variables = variables;
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
void solution_t<T>::print()
{
	for( int i = 0; i < variables.size(); i++ )
		printf("%6.3e ",variables[i]);
	printf("\n");
}

}
