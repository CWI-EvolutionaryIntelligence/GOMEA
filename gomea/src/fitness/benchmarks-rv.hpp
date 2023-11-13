#pragma once

#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/src/fitness/bbo_fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class sphereFunction_t : public GBOFitnessFunction_t<double>
{
	public:
		sphereFunction_t( int number_of_variables, double vtr );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

class sphereFunctionBBO_t : public BBOFitnessFunction_t<double>
{
	public:
		sphereFunctionBBO_t( int number_of_variables, double vtr );
		double objectiveFunction( int objective_index, vec_t<double> &variables );
};

class rosenbrockFunction_t : public GBOFitnessFunction_t<double>
{
	public:
		rosenbrockFunction_t( int number_of_variables, double vtr );
		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

class SOREBChainStrong_t : public GBOFitnessFunction_t<double>
{
	public:
		SOREBChainStrong_t( int number_of_variables, double vtr );
		~SOREBChainStrong_t();

		int getNumberOfSubfunctions();
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<double> &variables );
};

class rosenbrockFunctionBBO_t : public BBOFitnessFunction_t<double>
{
	public:
		rosenbrockFunctionBBO_t( int number_of_variables, double vtr );
		double objectiveFunction( int objective_index, vec_t<double> &variables );
};

class SOREBChainStrongBBO_t : public BBOFitnessFunction_t<double>
{
	public:
		SOREBChainStrongBBO_t( int number_of_variables, double vtr );
		double objectiveFunction( int objective_index, vec_t<double> &variables );
};

class circlesInASquareBBO_t : public BBOFitnessFunction_t<double>
{
	public:
		circlesInASquareBBO_t( int number_of_variables, double vtr );
		double objectiveFunction( int objective_index, vec_t<double> &variables );

		double getLowerRangeBound(int dimension);
		double getUpperRangeBound(int dimension);

	private:
		vec_t<double> optima;
		void initializeOptima();

};

}}
