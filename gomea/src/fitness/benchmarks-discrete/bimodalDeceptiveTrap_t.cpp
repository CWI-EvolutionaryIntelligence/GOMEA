/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

bimodalDeceptiveTrap_t::bimodalDeceptiveTrap_t( int number_of_variables, int trap_size ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	if(number_of_variables % trap_size != 0)
		throw std::runtime_error("Number of variables must be a multiple of trap size.");
	if(trap_size % 2 != 0)
		throw std::runtime_error("Bimodal deceptive trap function requires an even trap size.");
	if(trap_size < 4)
		throw std::runtime_error("Bimodal deceptive trap function requires a trap size of at least 4.");
	this->name = "Bimodal deceptive trap function";
	this->trap_size = trap_size;
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}

int bimodalDeceptiveTrap_t::getNumberOfSubfunctions() 
{
	return number_of_variables / trap_size;
}
		
double bimodalDeceptiveTrap_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	int trap_index = subfunction_index; 
	int unitation = 0;
	vec_t<int> inputs = inputsToSubfunction(subfunction_index);
	
	for( int ind : inputs )
		unitation += variables[ind];
	
	if( unitation == 0 || unitation == trap_size )
		return trap_size;
	else
		return trap_size - fabs(2*unitation - trap_size) - 2;
}

vec_t<int> bimodalDeceptiveTrap_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> inputs;
	int trap_index = subfunction_index; 
	for( int i = trap_index * trap_size; i < (trap_index+1)*trap_size; i++ )
		inputs.push_back(i);
	return( inputs );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
