#include "partial_solution.hpp"

namespace gomea{
namespace realvalued{

partial_solution_t::partial_solution_t( int num_touched_variables )
{
	this->num_touched_variables = num_touched_variables;
	this->touched_variables = vec(num_touched_variables, fill::none); 
	this->sample_zs = zeros<vec>(num_touched_variables);
	objective_value = 1e308;
	constraint_value = 1e308;
}

partial_solution_t::partial_solution_t( vec touched_variables, std::vector<int> &touched_indices )
{
	this->num_touched_variables = touched_variables.n_elem;
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
	this->sample_zs = zeros<vec>(num_touched_variables);
	objective_value = 1e308;
	constraint_value = 1e308;
}

partial_solution_t::partial_solution_t( vec touched_variables, vec sample_zs, std::vector<int> &touched_indices )
{
	this->num_touched_variables = touched_variables.n_elem;
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
	this->sample_zs = sample_zs;
	objective_value = 1e308;
	constraint_value = 1e308;
}
		
partial_solution_t::partial_solution_t( partial_solution_t &other )
{
	this->num_touched_variables = other.num_touched_variables;
	for( int i = 0; i < num_touched_variables; i++ )
		this->touched_indices.push_back(other.touched_indices[i]);
	this->touched_variables = other.touched_variables;
	this->sample_zs = other.sample_zs;
	objective_value = other.objective_value;
	constraint_value = other.constraint_value;
}

int partial_solution_t::getTouchedIndex( int ind )
{
	if( touched_index_map.empty() )
	{
		for( int i = 0; i < num_touched_variables; i++ )
			touched_index_map[touched_indices[i]] = i;
	}

	auto map_ind = touched_index_map.find(ind);
	if( map_ind == touched_index_map.end() )
		return( -1 );
	else
		return map_ind->second;
}

void partial_solution_t::setSampleMean( vec means )
{
	this->sample_means = means;
}

void partial_solution_t::print()
{
	for( int i = 0; i < num_touched_variables; i++ )
		printf("[%d][%6.3e]",touched_indices[i],touched_variables[i]);
	//printf("\n");
}

}}
