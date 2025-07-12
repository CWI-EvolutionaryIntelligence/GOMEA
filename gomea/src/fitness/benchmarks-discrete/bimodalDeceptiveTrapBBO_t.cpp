/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

bimodalDeceptiveTrapBBO_t::bimodalDeceptiveTrapBBO_t( int number_of_variables, int trap_size ) : BBOFitnessFunction_t<char>(number_of_variables)
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

double bimodalDeceptiveTrapBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double f = 0.0;
	for( int i = 0; i < number_of_variables; i += trap_size )
	{
		int unitation = 0;
		for( int j = 0; j < trap_size; j++ )
			unitation += variables[i+j];
		if( unitation == 0 || unitation == trap_size )
			f += trap_size;
		else
			f += trap_size - fabs(2*unitation - trap_size) - 2;
	}	
	return f;
}

double bimodalDeceptiveTrapBBO_t::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	return var_a / trap_size == var_b / trap_size ? 1.0 : 0.0;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
