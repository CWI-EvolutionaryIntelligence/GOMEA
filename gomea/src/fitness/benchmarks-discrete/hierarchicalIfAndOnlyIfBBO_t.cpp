/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

hierarchicalIfAndOnlyIfBBO_t::hierarchicalIfAndOnlyIfBBO_t( int number_of_variables ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	int symbols = number_of_variables;
	while( symbols > 1 )
	{
		if( symbols % 2 != 0 )
			throw std::runtime_error("Number of variables at each level must be a multiple of 2 for HIFF function.");
		symbols /= 2;
	}
	this->name = "Hierarchical if and only if (HIFF)";
	this->vtr = computeVTR();
	this->use_vtr = true;
	this->initialize();
}

double hierarchicalIfAndOnlyIfBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double result = 0.0;
	int block_size = 1;
	while( block_size <= number_of_variables )
	{
		for( int i = 0; i < number_of_variables; i += block_size )
		{
		bool same = true;
		for( int j = 0; j < block_size; j++ )
		{
			if( variables[i+j] != variables[i] )
			{
				same = false;
				break;
			}
		}
		if( same )
			result += block_size;
		}
		block_size *= 2;
	}

	return result;
}

double hierarchicalIfAndOnlyIfBBO_t::computeVTR()
{
	int logl = 0;
	int k = number_of_variables;
	while( k > 1 ){
		logl++;
		k /= 2;
	}
	return (logl+1) * number_of_variables;
}

double hierarchicalIfAndOnlyIfBBO_t::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	if( var_a == var_b )
		return 1.0;
	else{
		double result = 1.0;
		int k = 2;
		while( k < number_of_variables )
		{
			if(var_a / k == var_b / k)
				result *= 2;
			k *= 2;
		}
		return result;
	}
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
