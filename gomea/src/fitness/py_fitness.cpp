#include "gomea/src/fitness/py_fitness.hpp"

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
	int number_of_subfunctions = gomea_pyfitness_numberOfSubfunctions(py_class);
	if( number_of_subfunctions == -1 )
		throw std::runtime_error("FitnessFunction does not implement number_of_subfunctions().");
	return number_of_subfunctions;
}

vec_t<int> pyFitnessFunction_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies = gomea_pyfitness_inputsToSubfunction(py_class,subfunction_index);
	if( dependencies.size() == 0 )
		throw std::runtime_error("FitnessFunction does not implement inputsToSubfunction(int).");
	return dependencies;
}
		
double pyFitnessFunction_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double subf = gomea_pyfitness_subfunction(py_class,subfunction_index,variables);
	if( subf == 1e308 )
		throw std::runtime_error("FitnessFunction does not implement subfunction(int).");
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