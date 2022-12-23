#include "gomea/src/fitness/your_fitness_realvalued.hpp"

namespace gomea{
namespace fitness{

yourFitnessFunctionRealValued::yourFitnessFunctionRealValued( int number_of_variables, double vtr ) : customFitnessFunction_t(number_of_variables,vtr)
{
	this->name = "Your own fitness function (C++)";
	initialize();
}

int yourFitnessFunctionRealValued::getNumberOfSubfunctions()
{
	return number_of_variables-1;
}

vec_t<int> yourFitnessFunctionRealValued::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies;
	dependencies.push_back(subfunction_index);
	dependencies.push_back(subfunction_index+1);
	return dependencies;
}
		
double yourFitnessFunctionRealValued::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double x = variables[subfunction_index];
	double y = variables[subfunction_index+1];
	return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) );
}

double yourFitnessFunctionRealValued::getLowerRangeBound( int dimension )
{
	return( -1000 );
}
		
double yourFitnessFunctionRealValued::getUpperRangeBound( int dimension )
{
	return( 1000 );
}

}}