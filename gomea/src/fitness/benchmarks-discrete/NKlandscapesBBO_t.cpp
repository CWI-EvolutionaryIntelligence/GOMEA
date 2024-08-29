/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <fstream>

namespace gomea{
namespace fitness{

NKlandscapesBBO_t::NKlandscapesBBO_t( int number_of_variables, int K, long long fitness_table_seed ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "NK Landscapes function";
	this->vtr = -1;
	this->use_vtr = false;
	this->K = K;
	this->fitness_table_seed = fitness_table_seed;
	this->initializeFitnessTables();
	this->initialize();
}

void NKlandscapesBBO_t::initializeFitnessTables()
{
	static std::uniform_real_distribution<double> distribution(0.0,1.0);
	std::mt19937 rng(fitness_table_seed);

	fitness_tables.resize( number_of_variables - K + 1 );
	for( int i = 0; i < fitness_tables.size(); i++ )
	{
		fitness_tables[i].resize( 1<<K );
		for( int j = 0; j < (1<<K); j++ )
			fitness_tables[i][j] = (int)(1e6*distribution(rng));
	}
}
		
double NKlandscapesBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double fitness = 0;
	for( int i = 0; i < fitness_tables.size(); i++ )
	{
		size_t fitness_key = 0; // integer representation of bits that are used as input to the subfunction
		for( int j = 0; j < K; j++ )
		{
			int bit_index = (i + j);
			fitness_key = (fitness_key<<1) + variables[bit_index];
		}
		assert((int)fitness_key < (1<<K));
		fitness += fitness_tables[i][fitness_key];
	}
	return fitness;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
