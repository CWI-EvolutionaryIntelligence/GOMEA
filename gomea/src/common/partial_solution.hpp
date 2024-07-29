#pragma once

#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{

template<class T>
class solution_t;

template<class T>
class partial_solution_t
{
	public:
		vec_t<int> touched_indices;
		vec_t<T> touched_variables;
		std::set<int> touched_subfunctions;
		std::unordered_map<int,T> partial_objective_values;
		vec_t<double> fitness_buffers;

		partial_solution_t( int num_touched_variables );
		partial_solution_t( const vec_t<T> &touched_variables, const vec_t<int> &touched_indices );
		partial_solution_t( vec_t<double> &touched_variables, vec_t<double> &sample_zs, vec_t<int> &touched_indices );

		void initMemory( int number_of_fitness_buffers, int number_of_objectives );
		void initObjectiveValues( int number_of_objectives );
		void initFitnessBuffers( int number_of_fitness_buffers );
		
		int getNumberOfTouchedVariables();
		
		double getObjectiveValue( int objective_value_index = 0 ) const;
		const vec_t<double> getObjectiveValues() const;
		double getConstraintValue() const;

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setObjectiveValues( vec_t<double> objective_values ); 
		void setConstraintValue( double v );

		double getFitnessBuffer(int buffer_index) const;
		const vec_t<double> getFitnessBuffers() const;
		void setFitnessBuffers( vec_t<double> buffers ); 
		void resetFitnessBuffers();
		void addToFitnessBuffer( int buffer_index, double partial_fitness );
		void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );

		void insertSolution( solution_t<T> *solution );
		
		int getTouchedIndex( int ind );

		void setSampleMean( vec_t<T> &means );

		void print();

		vec_t<T> sample_zs; // Samples z~N(0,I), later transformed to N(mu,C)
		vec_t<T> sample_means;

		bool is_accepted = false;
		bool improves_elitist = false;

	private:
		vec_t<double> objective_values;
		double constraint_value;
		
		std::unordered_map<int,int> touched_index_map;
};

}
