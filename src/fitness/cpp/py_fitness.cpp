#include "fitness/py_fitness.hpp"

namespace gomea{
namespace fitness{

pyFitnessFunction_t::pyFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj ) : customFitnessFunction_t(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	initialize();
}

int pyFitnessFunction_t::getNumberOfSubfunctions()
{
	return number_of_parameters;
}

vec_t<int> pyFitnessFunction_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies;
	dependencies.push_back(subfunction_index);
	return dependencies;
}
		
double pyFitnessFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double subf = gomea_pyfitness_subfunction(py_class,subfunction_index,variables);
	return subf;
}

double pyFitnessFunction_t::getLowerRangeBound( int dimension )
{
	return( -1000 );
}
		
double pyFitnessFunction_t::getUpperRangeBound( int dimension )
{
	return( 1000 );
}

}}