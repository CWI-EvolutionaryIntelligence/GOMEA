#pragma once

#include <chrono>
#include <exception>
#include "gomea/src/utils/gomea_defs.hpp"

namespace gomea{
namespace utils{

	vec_t<int> randomPermutation( int size );
	vecE random1DNormalUnitVector( int length );
	double randomRealUniform01();
	int randomInt( int max );

	void initializeRandomNumberGenerator();
	void initializeRandomNumberGenerator( long long seed );

	extern std::mt19937 rng;
	extern long long random_seed;

	}
}

