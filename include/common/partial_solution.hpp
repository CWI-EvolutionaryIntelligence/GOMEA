#pragma once

#include <stdio.h>
#include <vector>
#include <map>

namespace gomea{
namespace common{

template<class T>
class partial_solution_t
{
	public:
		std::vector<int> touched_indices;
		std::vector<T> touched_variables;

		partial_solution_t( int num_touched_variables );
		partial_solution_t( std::vector<T> &touched_variables, std::vector<int> &touched_indices );
		partial_solution_t( partial_solution_t &other);

		int getNumberOfTouchedVariables();
		int getTouchedIndex( int ind );

		void print();

	private:
		std::vector<double> objective_values;
		double constraint_value;
		std::map<int,int> touched_index_map;
};

}}
