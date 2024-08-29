/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

NKlandscapes_t::NKlandscapes_t( int number_of_variables, int K, long long fitness_table_seed ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "NK Landscapes function";
	this->vtr = -1;
	this->use_vtr = false;
	this->K = K;
	this->fitness_table_seed = fitness_table_seed;
	this->initializeFitnessTables();
	this->initialize();
}

void NKlandscapes_t::initializeFitnessTables()
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

int NKlandscapes_t::getNumberOfSubfunctions() 
{
	return number_of_variables - K + 1;
}
		
double NKlandscapes_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	size_t fitness_key = 0; // integer representation of bits that are used as input to the subfunction
	for( int i = 0; i < K; i++ )
	{
		int bit_index = (subfunction_index + i);
		fitness_key = (fitness_key<<1) + variables[bit_index];
	}
	assert((int)fitness_key < (1<<K));
	return fitness_tables[subfunction_index][fitness_key];
}

vec_t<int> NKlandscapes_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> inputs;
	for( int i = 0; i < K; i++ )
		inputs.push_back( (subfunction_index + i) );
	return( inputs );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
