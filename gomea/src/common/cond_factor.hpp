#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{

class cond_factor_t{
	public:
		virtual ~cond_factor_t(){};
		cond_factor_t(); 
		cond_factor_t( const vec_t<int> &variables ); 
		cond_factor_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables );
		cond_factor_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables );

		vec_t<int> variables;
		vec_t<int> variables_conditioned_on;
		vec_t<vec_t<int>> frequency_tables;

		void addGroupOfVariables( const vec_t<int> &indices, const vec_t<int> &indices_cond = vec_t<int>() );
		void addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond );
		void addGroupOfVariables( int index, const vec_t<int> &indices_cond );
		void addGroupOfVariables( int index, int index_cond );

		void initializeFrequencyTables( const vec_t<solution_t<char>*> &population );
		vec_t<char> samplePartialSolutionConditional( solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index = -1 );
		vec_t<char> samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index = -1 );
		bool isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent );
		bool isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent );
		
		void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited );
		void setOrder( const vec_t<int> &order ); 
		void print();
};

}
