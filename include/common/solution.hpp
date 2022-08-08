#pragma once

#include "common/gomea_defs.hpp"
//#include "common/genotype.hpp"
//#include "common/full_genotype.hpp"

namespace gomea{

template<class T>
class solution_t
{
	public:
		solution_t( int number_of_variables );
		solution_t( vec_t<T> &variables );
		
		vec_t<T> variables;
	
		int getNumberOfVariables();
		int getNumberOfObjectives();
		double getObjectiveValue();
		double getObjectiveValue( int objective_value_index );
		double getPartialObjectiveValue( int subfunction_index );
		double getConstraintValue();
		double getPartialConstraintValue( int subfunction_index );

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setConstraintValue( double v );
		void setPartialObjectiveValue( int subfunction_index, double v );
		void setPartialConstraintValue( int subfunction_index, double v );

		vec_t<T> createPartialBackup(vec_t<int> variable_indices);
		void insertPartialBackup(vec_t<T> backup, vec_t<int> variable_indices);

		void print();

	private:
		vec_t<double> objective_values;
		vec_t<double> partial_objective_values;
		double constraint_value;
		vec_t<double> partial_constraint_values;
};

template class solution_t<int>;
template class solution_t<float>;
template class solution_t<double>;

}
