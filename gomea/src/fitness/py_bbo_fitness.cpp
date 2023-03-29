#include "gomea/src/fitness/py_bbo_fitness.hpp"

namespace gomea{
namespace fitness{

template<class T>
pyBBOFitnessFunction_t<T>::pyBBOFitnessFunction_t( int number_of_parameters, PyObject *obj ) : BBOFitnessFunction_t<T>(number_of_parameters)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->initialize();
}

template<class T>
pyBBOFitnessFunction_t<T>::pyBBOFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj ) : BBOFitnessFunction_t<T>(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->initialize();
}

template<>
double pyBBOFitnessFunction_t<char>::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double result = gomea_pyfitness_objective_function_bbo_discrete(py_class,objective_index,variables);
	if( result == 1e308 )
		throw std::runtime_error("FitnessFunction does not implement objective_function(int,vector[char]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<double>::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double result = gomea_pyfitness_objective_function_bbo_realvalued(py_class,objective_index,variables);
	if( result == 1e308 )
		throw std::runtime_error("FitnessFunction does not implement objective_function(int,vector[double]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<char>::constraintFunction( int objective_index, vec_t<char> &variables )
{
	double result = gomea_pyfitness_constraint_function_bbo_discrete(py_class,variables);
	if( result == 1e308 )
		throw std::runtime_error("FitnessFunction does not implement constraint_function(int,vector[char]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<double>::constraintFunction( int objective_index, vec_t<double> &variables )
{
	double result = gomea_pyfitness_constraint_function_bbo_realvalued(py_class,variables);
	if( result == 1e308 )
		throw std::runtime_error("FitnessFunction does not implement constraint_function(int,vector[double]).");
	return result;
}

template<class T>
double pyBBOFitnessFunction_t<T>::getLowerRangeBound( int dimension )
{
	assert(0);
	return( -1 );
}

template<class T>	
double pyBBOFitnessFunction_t<T>::getUpperRangeBound( int dimension )
{
	assert(0);
	return( -1 );
}

template<>
double pyBBOFitnessFunction_t<double>::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}

template<>	
double pyBBOFitnessFunction_t<double>::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

}}