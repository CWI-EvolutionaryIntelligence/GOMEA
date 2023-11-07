#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

sphereFunctionBBO_t::sphereFunctionBBO_t( int number_of_variables, double vtr ) : BBOFitnessFunction_t(number_of_variables, vtr)
{
	this->name = "Sphere function (BBO)";
	this->initialize();
}

double sphereFunctionBBO_t::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double result = 0.0;
	for( int i = 0; i < number_of_variables; i++ )
	{
		double x = variables[i];
		result += x * x;
	}
	return result;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
