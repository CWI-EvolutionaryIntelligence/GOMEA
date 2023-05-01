#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

rosenbrockFunctionBBO_t::rosenbrockFunctionBBO_t( int number_of_variables, double vtr ) : BBOFitnessFunction_t(number_of_variables, vtr)
{
	this->name = "Rosenbrock function (BBO)";
	this->initialize();
}

double rosenbrockFunctionBBO_t::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double result = 0.0;
	for( int i = 0; i < number_of_variables-1; i++ )
	{
		double x = variables[i];
		double y = variables[i+1];
		result += 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x);
	}
	return result;
}
		
}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
