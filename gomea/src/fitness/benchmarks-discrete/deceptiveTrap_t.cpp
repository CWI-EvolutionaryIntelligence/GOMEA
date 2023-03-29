/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

deceptiveTrap_t::deceptiveTrap_t( int number_of_variables, int trap_size ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "Deceptive trap function";
	this->trap_size = trap_size;
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}

int deceptiveTrap_t::getNumberOfSubfunctions() 
{
	return number_of_variables / trap_size;
}
		
double deceptiveTrap_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	int trap_index = subfunction_index; 
	int unitation = 0;
	vec_t<int> inputs = inputsToSubfunction(subfunction_index);
	
	for( int ind : inputs )
		unitation += variables[ind];
	
	if( unitation == trap_size )
		return unitation;
	else
		return trap_size - unitation - 1.0;
}

vec_t<int> deceptiveTrap_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> inputs;
	int trap_index = subfunction_index; 
	for( int i = trap_index * trap_size; i < (trap_index+1)*trap_size; i++ )
		inputs.push_back(i);
	return( inputs );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
