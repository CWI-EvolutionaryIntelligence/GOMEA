#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/partial_solution.hpp"

namespace gomea{

template<class T>
class solution_t
{
	public:
		solution_t( int number_of_variables );
		solution_t( vec_t<T> &variables );
		solution_t( size_t numberOfVariables_, size_t alphabetSize_ );

		bool operator==(const solution_t<T> &solutionB)
		{
			for (int i = 0; i < getNumberOfVariables(); ++i)
			{
				if (this->variables[i] != solutionB.variables[i])
					return false;
			}
			return true;
		}

		template<class U>
		friend std::ostream &operator<<(std::ostream &out, const solution_t<U> &solution);
		void randomInit(std::mt19937 *rng);

		void initMemory(int number_of_objectives, int number_of_fitness_buffers);

		int getNumberOfVariables() const;
		int getNumberOfObjectives() const;
		double getObjectiveValue() const;
		double getObjectiveValue( int objective_value_index ) const;
		vec_t<double> getObjectiveValues() const;
		double getPartialObjectiveValue( int subfunction_index ) const;
		double getConstraintValue() const;
		double getPartialConstraintValue( int subfunction_index ) const;

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setObjectiveValues( vec_t<double> v );
		void setConstraintValue( double v );
		void setPartialObjectiveValue( int subfunction_index, double v );
		void setPartialConstraintValue( int subfunction_index, double v );

		double getFitnessBuffer( int buffer_index );
		void addToFitnessBuffer( int buffer_index, double partial_fitness );
		void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );
		void setFitnessBuffers( vec_t<double> buffers );
		void clearFitnessBuffers();

		vec_t<T> createPartialBackup(vec_t<int> variable_indices);
		void insertVariables(vec_t<T> vars_to_insert, vec_t<int> indices_to_insert);
		void insertPartialSolution( partial_solution_t<T> *solution );

		void print();

		size_t getAlphabetSize();
		
		vec_t<T> variables;
		vec_t<double> fitness_buffers;

	private:
		void initObjectiveValues(int number_of_objectives);
		void initFitnessBuffers(int number_of_fitness_buffers);

		vec_t<double> objective_values;
		double constraint_value;
		
		vec_t<double> partial_objective_values;
		vec_t<double> partial_constraint_values;

		size_t alphabetSize = 0;
};

template class solution_t<char>;
template class solution_t<int>;
template class solution_t<float>;
template class solution_t<double>;

}
