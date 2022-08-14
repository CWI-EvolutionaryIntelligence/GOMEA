#include "fitness/your_fitness.hpp"

namespace gomea{
namespace fitness{

yourFitnessFunction_t::yourFitnessFunction_t( int number_of_parameters, double vtr ) : customFitnessFunction_t(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (C++)";
	initialize();
}

int yourFitnessFunction_t::getNumberOfSubfunctions()
{
	return number_of_parameters;
}

vec_t<int> yourFitnessFunction_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies;
	dependencies.push_back(subfunction_index);
	return dependencies;
}
		
double yourFitnessFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double x = variables[subfunction_index];
	return( x*x );
	//return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) );
}

double yourFitnessFunction_t::getLowerRangeBound( int dimension )
{
	return( -1000 );
}
		
double yourFitnessFunction_t::getUpperRangeBound( int dimension )
{
	return( 1000 );
}

}}