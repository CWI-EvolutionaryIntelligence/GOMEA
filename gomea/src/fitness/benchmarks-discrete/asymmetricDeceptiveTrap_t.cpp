/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

asymmetricDeceptiveTrap_t::asymmetricDeceptiveTrap_t( int number_of_variables, int trap_size ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	if(number_of_variables % trap_size != 0)
		throw std::runtime_error("Number of variables must be a multiple of trap size.");
	if(trap_size < 4)
		throw std::runtime_error("Asymmetric deceptive trap function requires a trap size of at least 4.");
	this->name = "Asymmetric deceptive trap function";
	this->trap_size = trap_size;
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}

int asymmetricDeceptiveTrap_t::getNumberOfSubfunctions() 
{
	return number_of_variables / trap_size;
}
		
double asymmetricDeceptiveTrap_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	int trap_index = subfunction_index; 
	int unitation = 0;
	vec_t<int> inputs = inputsToSubfunction(subfunction_index);
	
	for( int ind : inputs )
		unitation += variables[ind];
	
	if( unitation == trap_size-1 && variables[inputs[0]] == 0)
		return trap_size;
	else
		return trap_size * ((double) (trap_size - unitation))/((double) (trap_size+1));
}

vec_t<int> asymmetricDeceptiveTrap_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> inputs;
	int trap_index = subfunction_index; 
	for( int i = trap_index * trap_size; i < (trap_index+1)*trap_size; i++ )
		inputs.push_back(i);
	return( inputs );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
