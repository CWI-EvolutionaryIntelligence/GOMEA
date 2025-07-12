/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

asymmetricDeceptiveTrapBBO_t::asymmetricDeceptiveTrapBBO_t( int number_of_variables, int trap_size ) : BBOFitnessFunction_t<char>(number_of_variables)
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

double asymmetricDeceptiveTrapBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double f = 0.0;
	for( int i = 0; i < number_of_variables; i += trap_size )
	{
		int unitation = 0;
		for( int j = 0; j < trap_size; j++ )
			unitation += variables[i+j];
		if( unitation == trap_size-1 && variables[i] == 0)
			return trap_size;
		else
			return trap_size * ((double) (trap_size - unitation))/((double) (trap_size+1));
	}	
	return f;
}

double asymmetricDeceptiveTrapBBO_t::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	return var_a / trap_size == var_b / trap_size ? 1.0 : 0.0;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
