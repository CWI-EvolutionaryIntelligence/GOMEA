#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/common/gomea_defs.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{

class factorization_t{
	public:
		virtual ~factorization_t(){};
		factorization_t(); 
		factorization_t( const vec_t<int> &variables ); 
		factorization_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables );
		factorization_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables );

		vec_t<int> variables;
		vec_t<int> order;
		vec_t<vec_t<int>> variable_groups;
		vec_t<vec_t<int>> variables_conditioned_on;
		vec_t<vec_t<int>> index_in_var_array;

		void addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond );
		void addGroupOfVariables( vec_t<int> indices, vec_t<int> indices_cond );
		void addGroupOfVariables( int index, const vec_t<int> &indices_cond );
		void addGroupOfVariables( int index, int index_cond );
		
		void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited );
		void setOrder( const vec_t<int> &order ); 
		void print();
};

}
