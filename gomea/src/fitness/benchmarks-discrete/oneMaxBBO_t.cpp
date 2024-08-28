/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

oneMaxBBO_t::oneMaxBBO_t( int number_of_variables ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "OneMax function";
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}
		
double oneMaxBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double result = 0.0;
	for( char c : variables )
		result += c;
	return result;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
