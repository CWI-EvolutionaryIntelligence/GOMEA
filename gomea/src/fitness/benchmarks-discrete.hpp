#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/src/fitness/bbo_fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class oneMax_t: public GBOFitnessFunction_t<char>
{
	public:
		oneMax_t( int number_of_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		double subfunction( int subfunction_index, vec_t<char> &variables );
};

class deceptiveTrap_t: public GBOFitnessFunction_t<char>
{
	public:
		deceptiveTrap_t( int number_of_variables, int trap_size );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		int trap_size;
		double subfunction( int subfunction_index, vec_t<char> &variables );
};

class deceptiveTrapBBO_t: public BBOFitnessFunction_t<char>
{
	public:
		deceptiveTrapBBO_t( int number_of_variables, int trap_size );
		double objectiveFunction( int objective_index, vec_t<char> &variables );
		
	private:
		int trap_size;
};

class maxCut_t: public GBOFitnessFunction_t<char>
{
	public:
		maxCut_t( std::string input_file, std::string vtr_file );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		
	private:
		vec_t<vec_t<int>> edges;
		vec_t<double> edge_weights;

		void readInputFile( std::string input_file );
		void readVTRFile( std::string input_file );
		double subfunction( int subfunction_index, vec_t<char> &variables );
};

class maxCutBBO_t: public BBOFitnessFunction_t<char>
{
	public:
		maxCutBBO_t( std::string input_file, std::string vtr_file );
		double objectiveFunction( int objective_index, vec_t<char> &variables );
		
	private:
		vec_t<vec_t<int>> edges;
		vec_t<double> edge_weights;

		void readInputFile( std::string input_file );
		void readVTRFile( std::string input_file );
};

}}
