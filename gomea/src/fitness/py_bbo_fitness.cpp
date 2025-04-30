#ifndef CPP_STANDALONE

#include "gomea/src/fitness/py_bbo_fitness.hpp"

namespace gomea{
namespace fitness{

template<>
pyBBOFitnessFunction_t<char>::pyBBOFitnessFunction_t( int number_of_parameters, int alphabet_size, PyObject *obj ) : BBOFitnessFunction_t<char>(number_of_parameters)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->alphabet_size = alphabet_size;
	this->initialize();
}

template<>
pyBBOFitnessFunction_t<char>::pyBBOFitnessFunction_t( int number_of_parameters, int alphabet_size, double vtr, PyObject *obj ) : BBOFitnessFunction_t<char>(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->alphabet_size = alphabet_size;
	this->initialize();
}

template<>
pyBBOFitnessFunction_t<double>::pyBBOFitnessFunction_t( int number_of_parameters, PyObject *obj ) : BBOFitnessFunction_t<double>(number_of_parameters)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->initialize();
}

template<>
pyBBOFitnessFunction_t<double>::pyBBOFitnessFunction_t( int number_of_parameters, double vtr, PyObject *obj ) : BBOFitnessFunction_t<double>(number_of_parameters,vtr)
{
	this->name = "Your own fitness function (Python)";
	this->py_class = obj;
	this->initialize();
}

template<>
double pyBBOFitnessFunction_t<char>::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double result = gomea_pyfitness_objective_function_bbo_discrete(py_class,objective_index,variables);
	if( result == INFINITY )
		throw std::runtime_error("FitnessFunction does not implement objective_function(int,vector[char]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<double>::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double result = gomea_pyfitness_objective_function_bbo_realvalued(py_class,objective_index,variables);
	if( result == INFINITY )
		throw std::runtime_error("FitnessFunction does not implement objective_function(int,vector[double]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<char>::constraintFunction( vec_t<char> &variables )
{
	double result = gomea_pyfitness_constraint_function_bbo_discrete(py_class,variables);
	if( result == INFINITY )
		throw std::runtime_error("FitnessFunction does not implement constraint_function(vector[char]).");
	return result;
}

template<>
double pyBBOFitnessFunction_t<double>::constraintFunction( vec_t<double> &variables )
{
	double result = gomea_pyfitness_constraint_function_bbo_realvalued(py_class,variables);
	if( result == INFINITY )
		throw std::runtime_error("FitnessFunction does not implement constraint_function(vector[double]).");
	return result;
}

template<class T>
double pyBBOFitnessFunction_t<T>::getLowerRangeBound( int variable_index )
{
	return gomea_pyfitness_lower_range_bound(py_class,variable_index);
}

template<class T>	
double pyBBOFitnessFunction_t<T>::getUpperRangeBound( int variable_index )
{
	return gomea_pyfitness_upper_range_bound(py_class,variable_index);
}

template<>
double pyBBOFitnessFunction_t<char>::getLowerRangeBound( int variable_index )
{
	double result = gomea_pyfitness_lower_range_bound(py_class,variable_index);
	if( result == -INFINITY ){
		return 0;
	}
	return result;
}

template<>	
double pyBBOFitnessFunction_t<char>::getUpperRangeBound( int variable_index )
{
	double result = gomea_pyfitness_upper_range_bound(py_class,variable_index);
	if( result == INFINITY ){
		return alphabet_size-1;
	}
	return result;
}

template class pyBBOFitnessFunction_t<char>;
template class pyBBOFitnessFunction_t<double>;

}}

#endif