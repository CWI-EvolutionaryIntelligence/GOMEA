#pragma once

#include <stdio.h>
#include <vector>
#include "gomea_defs.hpp"

namespace gomea{
namespace common{

template<class T>
class solution_t
{
	public:
		solution_t();
		solution_t( int number_of_variables );
		solution_t( vec_t<T> &variables );
		
		vec_t<T> variables;
	
		int getNumberOfVariables();
		int getNumberOfObjectives();
		double getObjectiveValue();
		double getObjectiveValue( int objective_value_index );
		double getConstraintValue();

		void setObjectiveValue( double v );
		void setObjectiveValue( int objective_value_index, double v );
		void setConstraintValue( double v );

		void print();

	private:
		vec_t<double> objective_values;
		vec_t<double> partial_objective_values;
		double constraint_value;
		vec_t<double> partial_constraint_values;
};

}}
