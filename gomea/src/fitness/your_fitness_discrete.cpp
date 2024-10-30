#include "gomea/src/fitness/your_fitness_discrete.hpp"

namespace gomea{
namespace fitness{

yourFitnessFunctionDiscrete::yourFitnessFunctionDiscrete( int number_of_variables, int alphabet_size, double vtr ) : GBOFitnessFunction_t(number_of_variables,vtr) 
{
	this->name = "Your own fitness function (C++)";
	this->vtr = vtr;
	this->alphabet_size = alphabet_size;
	assert( number_of_variables % trap_size == 0 );
	this->initialize();
}

int yourFitnessFunctionDiscrete::getNumberOfSubfunctions()
{
	return number_of_variables / trap_size;
}

vec_t<int> yourFitnessFunctionDiscrete::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies;
	for( int i = 0; i < trap_size; i++ )
		dependencies.push_back(subfunction_index*trap_size + i);
	return dependencies;
}

double yourFitnessFunctionDiscrete::subfunction( int subfunction_index, vec_t<char> &variables )
{
	int unitation = 0;
	for( int i = 0; i < trap_size; i++ )
		unitation += variables[trap_size * subfunction_index + i];
	if( unitation == trap_size )
		return trap_size;
	return( (double) (trap_size - 1 - unitation) );
}

}}