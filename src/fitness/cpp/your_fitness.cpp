#include "fitness/your_fitness.hpp"

namespace gomea{
namespace fitness{

yourFitnessFunction_t::yourFitnessFunction_t( int number_of_parameters, double vtr ) : customFitnessFunction_t(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (C++)";
}
		
double yourFitnessFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	return( variables[subfunction_index] * variables[subfunction_index] );
}

double yourFitnessFunction_t::getLowerRangeBound( int dimension )
{
	return( -100 );
}
		
double yourFitnessFunction_t::getUpperRangeBound( int dimension )
{
	return( 100 );
}

}}