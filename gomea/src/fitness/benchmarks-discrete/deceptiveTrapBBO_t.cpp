/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

deceptiveTrapBBO_t::deceptiveTrapBBO_t( int number_of_variables, int trap_size ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "Deceptive trap function";
	this->trap_size = trap_size;
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}

double deceptiveTrapBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double f = 0.0;
	for( int i = 0; i < number_of_variables; i += trap_size )
	{
		int unitation = 0;
		for( int j = 0; j < trap_size; j++ )
			unitation += variables[i+j];
		if( unitation == trap_size )
			f += unitation;
		else
			f += trap_size - unitation - 1;
	}	
	return f;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
