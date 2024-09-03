#pragma once

#include <chrono>
#include <exception>
#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{
	namespace utils{

		class terminationException : public std::runtime_error
		{
			public:
				terminationException(std::string message_) : runtime_error(message_) {}
		};

	double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
	double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
	double vectorDotProduct( double *vector0, double *vector1, int n0 );

	vec_t<int> getSortedOrder( vec_t<int> &data );

	bool isPowerOfK(int n, int k);

	double distanceEuclidean( double *solution_a, double *solution_b, int n );
	double distanceEuclidean( vec_t<double> &x, vec_t<double> &y );
	double distanceEuclidean2D( double x1, double y1, double x2, double y2 );

	vec_t<int> randomPermutation( int size );
	double randomRealUniform01();
	int randomInt( int max );

	void initializeRandomNumberGenerator();
	void initializeRandomNumberGenerator( long long seed );

	extern std::mt19937 rng;
	extern long long random_seed;

	}
}

