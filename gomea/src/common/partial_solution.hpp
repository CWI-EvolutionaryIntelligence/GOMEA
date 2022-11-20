#pragma once

#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{

template<class T>
class partial_solution_t
{
	public:
		vec_t<int> touched_indices;
		vec_t<T> touched_variables;
		std::set<int> touched_subfunctions;
		std::map<int,T> partial_objective_values;
		vec_t<double> fitness_buffers;

		partial_solution_t( int num_touched_variables );
		partial_solution_t( vec_t<T> &touched_variables, vec_t<int> &touched_indices );

		void init( int number_of_fitness_buffers, int number_of_objectives );
		int getNumberOfTouchedVariables();
		
		double getObjectiveValue();
		double getObjectiveValue( int objective_value_index );
		vec_t<double> getObjectiveValues();
		double getConstraintValue();

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setConstraintValue( double v );

		double getFitnessBuffer(int buffer_index);
		void setFitnessBuffers( vec_t<double> buffers ); 
		void resetFitnessBuffers();
		void addToFitnessBuffer( int buffer_index, double partial_fitness );
		void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );
		
		int getTouchedIndex( int ind );

		void print();

	private:
		void initObjectiveValues( int number_of_objectives );
		void initFitnessBuffers( int number_of_fitness_buffers );
		
		vec_t<double> objective_values;
		double constraint_value;
		
		std::map<int,int> touched_index_map;
};

template class partial_solution_t<char>;
template class partial_solution_t<int>;
template class partial_solution_t<float>;
template class partial_solution_t<double>;

}
