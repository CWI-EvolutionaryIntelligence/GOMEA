#include "gomea/src/common/factorization.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{

factorization_t::factorization_t() 
{
}

factorization_t::factorization_t( const vec_t<int> &variables )
{
	this->variables = variables;
}

factorization_t::factorization_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

factorization_t::factorization_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

void factorization_t::addGroupOfVariables( vec_t<int> indices, vec_t<int> indices_cond )
{
	std::sort(indices.begin(),indices.end());
	std::sort(indices_cond.begin(),indices_cond.end());
	vec_t<int> indices_map;
	for( int i : indices )
	{
		indices_map.push_back(variables.size());
		variables.push_back(i);
	}
	//std::sort(variables.begin(),variables.end());
	index_in_var_array.push_back(indices_map);
	variable_groups.push_back(indices);
	variables_conditioned_on.push_back(indices_cond);
	order.push_back(order.size());
}

void factorization_t::addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond )
{
	vec_t<int> cond;
	for( int i : indices_cond )
		cond.push_back(i);
	addGroupOfVariables(indices,cond);
}
			
void factorization_t::addGroupOfVariables( int index, const vec_t<int> &indices_cond )
{
	vec_t<int> indices;
	indices.push_back(index);
	addGroupOfVariables(indices, indices_cond);
}
			
void factorization_t::addGroupOfVariables( int index, int index_cond )
{
	vec_t<int> indices, indices_cond;
	indices.push_back(index);
	indices_cond.push_back(index_cond);
	addGroupOfVariables(indices, indices_cond);
}

void factorization_t::updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited ) 
{
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	
	variables_conditioned_on.clear();
	variables_conditioned_on.resize(variable_groups.size());
	
	for( size_t i = 0; i < variable_groups.size(); i++ )
	{
		int ind = i;
		if( order.size() > 0 )
			ind = order[i];
		vec_t<int> clique = variable_groups[ind];

		// Add FOS element of all nodes in clique, conditioned on dependent, already visited variables
		std::set<int> cond;
		for( int v : clique )
			visited[v] = IN_CLIQUE;
		for( int v : clique )
		{
			for( int x : variable_interaction_graph.at(v) )
			{
				if( visited[x] == IS_VISITED )
					cond.insert(x);
			}
		}
		for( int v : clique )
			visited[v] = IS_VISITED;
		
		vec_t<int> cond_vec;
		for( int x : cond )
			cond_vec.push_back(x);
		variables_conditioned_on[ind] = cond_vec;
	}
}

void factorization_t::setOrder( const vec_t<int> &order ) 
{
	this->order = order;
}

void factorization_t::print()
{
	for(size_t i = 0; i < variable_groups.size(); i++ )
	{
		int og = order[i];
		printf("[");
		for( int x : variable_groups[og] )
			printf(" %d",x);
		printf("]->[");
		for( int x : variables_conditioned_on[og] )
			printf(" %d",x);
		printf("],");
	}
	printf("\n");
}

}
