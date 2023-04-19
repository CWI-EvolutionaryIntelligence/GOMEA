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

		friend std::ostream &operator<<(std::ostream &out, const solution_t<T> &solution)
		{
			for (int i = 0; i < solution.getNumberOfVariables(); ++i)
				out << +solution.variables[i];
			out << " | " << solution.getObjectiveValue();
			return out;
		}
		
		void randomInit(std::mt19937 *rng);

		void initMemory(int number_of_objectives, int number_of_fitness_buffers);
		void initObjectiveValues(int number_of_objectives);
		void initFitnessBuffers(int number_of_fitness_buffers);

		int getNumberOfVariables() const;
		int getNumberOfObjectives() const;
		double getObjectiveValue( int objective_value_index = 0 ) const;
		const vec_t<double> getObjectiveValues() const;
		double getPartialObjectiveValue( int subfunction_index ) const;
		double getConstraintValue() const;
		double getPartialConstraintValue( int subfunction_index ) const;

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setObjectiveValues( const vec_t<double> &v );
		void setConstraintValue( double v );
		void setPartialObjectiveValue( int subfunction_index, double v );
		void setPartialConstraintValue( int subfunction_index, double v );

		double getFitnessBuffer( int buffer_index ) const;
		const vec_t<double> getFitnessBuffers() const;
		void addToFitnessBuffer( int buffer_index, double partial_fitness );
		void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );
		void setFitnessBuffers( const vec_t<double> &buffers );
		void clearFitnessBuffers();

		partial_solution_t<T> getPartialCopy( const vec_t<int> &variable_indices ) const;
		const vec_t<T> getCopyOfVariables( const vec_t<int> &variable_indices = vec_t<int>()) const;
		void insertVariables( const vec_t<T> &vars_to_insert );
		void insertVariables(vec_t<T> vars_to_insert, vec_t<int> indices_to_insert);
		void insertSolution( solution_t<T> *solution );
		void insertPartialSolution( partial_solution_t<T> *solution );

		void print();

		size_t getAlphabetSize();
		
		vec_t<T> variables;
		vec_t<double> fitness_buffers;

	private:
		vec_t<double> objective_values;
		double constraint_value;
		
		vec_t<double> partial_objective_values;
		vec_t<double> partial_constraint_values;

		size_t alphabetSize = 0;
};

}
