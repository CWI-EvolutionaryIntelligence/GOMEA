#include "solution.hpp"

namespace gomea{
namespace realvalued{

solution_t::solution_t( int number_of_variables )
{
	this->number_of_variables = number_of_variables;
	variables = vec(number_of_variables, fill::none);
	objective_value = 1e308;
	constraint_value = 1e308;
}

solution_t::solution_t( vec &variables )
{
	this->variables = variables;
	this->number_of_variables = variables.n_elem;
	objective_value = 1e308;
	constraint_value = 1e308;
}

solution_t::~solution_t()
{
}

solution_t::solution_t( const solution_t &sol )
{
	number_of_variables = sol.number_of_variables;
	variables = sol.variables;
	objective_value = sol.objective_value;
	constraint_value = sol.constraint_value;
	NIS = sol.NIS;
}

void solution_t::print()
{
	for( int i = 0; i < number_of_variables; i++ )
		printf("%6.3e ",variables[i]);
	printf("\n");
}

}}
