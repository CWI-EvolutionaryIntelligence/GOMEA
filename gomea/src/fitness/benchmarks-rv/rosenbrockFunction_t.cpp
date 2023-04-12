#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

using namespace gomea;

rosenbrockFunction_t::rosenbrockFunction_t( int number_of_variables, double vtr ) : GBOFitnessFunction_t(number_of_variables,vtr)
{
	this->name = "Rosenbrock function";
	this->initialize();
}

int rosenbrockFunction_t::getNumberOfSubfunctions()
{
	return number_of_variables-1;
}

vec_t<int> rosenbrockFunction_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec(2);
	vec.push_back(subfunction_index);
	vec.push_back(subfunction_index+1);
	return vec;
}

double rosenbrockFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double x = variables[subfunction_index];
	double y = variables[subfunction_index+1];
	return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
