#include "gomea/src/common/solution.hpp"

namespace gomea{

template<class T>
solution_t<T>::solution_t( int number_of_variables ) : variables(vec_t<T>(number_of_variables)) {}

template<class T>
solution_t<T>::solution_t( vec_t<T> &variables ) : variables(vec_t<T>(variables)){}

template<>
solution_t<char>::solution_t(size_t numberOfVariables_, size_t alphabetSize_) : solution_t(numberOfVariables_)
{
	this->alphabetSize = alphabetSize_;
	fill(variables.begin(), variables.end(), 0);
}

template<class T>
std::ostream & operator << (std::ostream &out, const solution_t<T> &individual)
{
	for (int i = 0; i < individual.getNumberOfVariables(); ++i)
		out << +individual.variables[i];
	out << " | " << individual.getObjectiveValue();
	return out;
}


template<class T>
void solution_t<T>::randomInit(std::mt19937 *rng)
{
	char buf[256];
	sprintf(buf,"solution_t<%s> does not implement randomInit(mt19937*)\n",typeid(T).name());
	throw std::runtime_error(buf);
}

template<>
void solution_t<char>::randomInit(std::mt19937 *rng)
{
	for (int i = 0; i < getNumberOfVariables(); ++i)
	{
		variables[i] = (*rng)() % alphabetSize;
	}
}

template<class T>
void solution_t<T>::initMemory( int number_of_objectives, int number_of_fitness_buffers )
{
	if( objective_values.size() == 0 )
		initObjectiveValues( number_of_objectives );
	if( fitness_buffers.size() == 0 )
		initFitnessBuffers( number_of_fitness_buffers );
}

template<class T>
void solution_t<T>::initObjectiveValues( int number_of_objectives )
{
 	this->objective_values.resize(number_of_objectives);
	for( int i = 0; i < number_of_objectives; i++ )
		this->objective_values[i] = INFINITY;
}

template<class T>
void solution_t<T>::initFitnessBuffers( int number_of_fitness_buffers ) 
{
	this->fitness_buffers.resize(number_of_fitness_buffers);
	for( int i = 0; i < number_of_fitness_buffers; i++ )
		this->fitness_buffers[i] = INFINITY;
}

template<class T>
int solution_t<T>::getNumberOfVariables() const
{
	return variables.size();
}

template<class T>
int solution_t<T>::getNumberOfObjectives() const
{
	return objective_values.size();
}

template<class T>
size_t solution_t<T>::getAlphabetSize()
{
	char buf[256];
	sprintf(buf,"solution_t<%s> does not implement getAlphabetSize()\n",typeid(T).name());
	throw std::runtime_error(buf);
}

template<>
size_t solution_t<char>::getAlphabetSize()
{
	return alphabetSize;
}

template<class T>
double solution_t<T>::getObjectiveValue( int objective_value_index ) const
{
	return objective_values[objective_value_index];
}

template<class T>
const vec_t<double> solution_t<T>::getObjectiveValues() const
{
	return objective_values;
}

template<class T>
double solution_t<T>::getConstraintValue() const
{
	return constraint_value;
}

template<class T>
double solution_t<T>::getPartialObjectiveValue( int subfunction_index ) const
{
	return partial_objective_values[subfunction_index];
}

template<class T>
double solution_t<T>::getPartialConstraintValue( int subfunction_index ) const
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
void solution_t<T>::setObjectiveValues( const vec_t<double> &v )
{
	for( size_t i = 0; i < objective_values.size(); i++ )
		setObjectiveValue(i, v[i]);
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
double solution_t<T>::getFitnessBuffer( int buffer_index ) const
{
	return( fitness_buffers[buffer_index] );
}

template<class T>		
const vec_t<double> solution_t<T>::getFitnessBuffers() const
{
	return( fitness_buffers );
}

template<class T>
void solution_t<T>::addToFitnessBuffer( int buffer_index, double partial_fitness )
{
	fitness_buffers[buffer_index] += partial_fitness;	
}

template<class T>
void solution_t<T>::subtractFromFitnessBuffer( int buffer_index, double partial_fitness )
{
	fitness_buffers[buffer_index] -= partial_fitness;	
}

template<class T>		
void solution_t<T>::setFitnessBuffers( const vec_t<double> &buffers )
{
	assert( this->fitness_buffers.size() == buffers.size() );
	for( size_t i = 0; i < this->fitness_buffers.size(); i++ )
	{
		this->fitness_buffers[i] = buffers[i];
	}
}

template<class T>
void solution_t<T>::clearFitnessBuffers()
{
	for( size_t i = 0; i < fitness_buffers.size(); i++ )
	{
		fitness_buffers[i] = 0.0;
	}
}

template<class T>
partial_solution_t<T> solution_t<T>::getPartialCopy( const vec_t<int> &variable_indices ) const
{
	vec_t<T> backup_vars = getCopyOfVariables(variable_indices);
	partial_solution_t<T> backup_solution = partial_solution_t<T>(backup_vars, variable_indices);
	backup_solution.setObjectiveValues(getObjectiveValues());
	backup_solution.setConstraintValue(getConstraintValue());
	backup_solution.setFitnessBuffers(getFitnessBuffers());
	return( backup_solution );
}

template<class T>
const vec_t<T> solution_t<T>::getCopyOfVariables( const vec_t<int> &variable_indices ) const
{
	if( variable_indices.size() == 0 )
	{
		vec_t<T> backup = vec_t<T>(getNumberOfVariables());
		for (size_t i = 0; i < getNumberOfVariables(); i++)
		{
			backup[i] = variables[i];
		}
		return backup;
	}
	else
	{
		vec_t<T> backup = vec_t<T>(variable_indices.size());
		for (size_t i = 0; i < variable_indices.size(); i++)
		{
			int ind = variable_indices[i];
			backup[i] = variables[ind];
		}
		return backup;
	}
}

template<class T>
void solution_t<T>::insertVariables( const vec_t<T> &vars_to_insert )
{
	assert( vars_to_insert.size() == variables.size() );
	for( size_t i = 0; i < variables.size(); i++ )
	{
		variables[i] = vars_to_insert[i];
	}
}

template<class T>
void solution_t<T>::insertVariables( vec_t<T> vars_to_insert, vec_t<int> indices_to_insert )
{
	for( size_t i = 0; i < indices_to_insert.size(); i++ )
	{
		int ind = indices_to_insert[i];
		variables[ind] = vars_to_insert[i];
	}
}

template<class T>
void solution_t<T>::insertSolution( solution_t<T> *solution )
{
	insertVariables(solution->variables);
	setObjectiveValues( solution->getObjectiveValues() );
	setConstraintValue( solution->getConstraintValue() );
	setFitnessBuffers( solution->fitness_buffers );
}

template<class T>
void solution_t<T>::insertPartialSolution( partial_solution_t<T> *solution )
{
	insertVariables(solution->touched_variables,solution->touched_indices);
	setObjectiveValues( solution->getObjectiveValues() );
	setConstraintValue( solution->getConstraintValue() );
	setFitnessBuffers( solution->fitness_buffers );
}

template<class T>
void solution_t<T>::print()
{
	for( size_t i = 0; i < variables.size(); i++ )
		printf("%6.3e ",(double)variables[i]);
	printf("\n");
}

template<>
void solution_t<char>::print()
{
	for( size_t i = 0; i < variables.size(); i++ )
		printf("%c ",variables[i]);
	printf("\n");
}

template<>
void solution_t<int>::print()
{
	for( size_t i = 0; i < variables.size(); i++ )
		printf("%d ",variables[i]);
	printf("\n");
}

}
