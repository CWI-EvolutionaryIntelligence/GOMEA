#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

sphereFunction_t::sphereFunction_t( int number_of_variables, double vtr ) : GBOFitnessFunction_t(number_of_variables,vtr)
{
	this->name = "Sphere function";
	this->initialize();
}

int sphereFunction_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> sphereFunction_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec(1);
	vec.push_back(subfunction_index);
	return vec;
}
		
double sphereFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double x = variables[subfunction_index];
	return( x * x );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
