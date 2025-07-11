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
double pyBBOFitnessFunction_t<T>::getLowerRangeBound( int dimension )
{
	throw std::runtime_error("FitnessFunction does not implement getLowerRangeBound(int).");
}

template<class T>	
double pyBBOFitnessFunction_t<T>::getUpperRangeBound( int dimension )
{
	throw std::runtime_error("FitnessFunction does not implement getUpperRangeBound(int).");
}

template<>
double pyBBOFitnessFunction_t<char>::getLowerRangeBound( int dimension )
{
	return( 0 );
}

template<>	
double pyBBOFitnessFunction_t<char>::getUpperRangeBound( int dimension )
{
	return( alphabet_size-1 );
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

template<class T>
double pyBBOFitnessFunction_t<T>::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	double result = gomea_pyfitness_similarity_measure(py_class,var_a,var_b);
	if( result < 0.0 )
	{
		throw std::runtime_error("Fitness function does not implement similarity_measure(size_t,size_t).");
	}
	else
	{
		return result;
	}
}

template class pyBBOFitnessFunction_t<char>;
template class pyBBOFitnessFunction_t<double>;

}}

#endif