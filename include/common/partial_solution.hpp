#pragma once

#include "common/gomea_defs.hpp"

namespace gomea{

template<class T>
class partial_solution_t
{
	public:
		vec_t<int> touched_indices;
		vec_t<T> touched_variables;
		std::map<int,T> partial_objective_values;

		partial_solution_t( int num_touched_variables );
		partial_solution_t( vec_t<T> &touched_variables, vec_t<int> &touched_indices );

		int getNumberOfTouchedVariables();
		
		double getObjectiveValue();
		double getObjectiveValue( int objective_value_index );
		double getConstraintValue();

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setConstraintValue( double v );
		
		int getTouchedIndex( int ind );

		void print();

	private:
		vec_t<double> objective_values;
		double constraint_value;
		std::map<int,int> touched_index_map;
};

template class partial_solution_t<int>;
template class partial_solution_t<float>;
template class partial_solution_t<double>;

}
