#pragma once

#include <chrono>
#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{
	namespace utils{

	double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
	double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
	double vectorDotProduct( double *vector0, double *vector1, int n0 );

	vec_t<int> getSortedOrder( vec_t<int> &data );

	bool isPowerOfK(int n, int k);

	vec_t<int> randomPermutation( int size );
	double randomRealUniform01();
	int randomInt( int max );

	void initializeRandomNumberGenerator();
	void initializeRandomNumberGenerator( long long seed );

	static std::mt19937 rng;
	extern long long random_seed;

	}
}

