#pragma once

#include <chrono>
#include <exception>
#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{
namespace utils{

	class customException : public std::exception
	{
		private:
			std::string message;

		public:
			customException(std::string message_) : message(message_) {}
			const char *what() const throw()
			{
				return message.c_str();
			}
	};

	double **matrixNew( int n, int m );
	double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
	double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
	double vectorDotProduct( double *vector0, double *vector1, int n0 );

	void *Malloc( long size );
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

	int *mergeSort( double *array, int array_size );
	void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
	void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q );
	void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
	int *mergeSortInt( int *array, int array_size );
	void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q );

	double normalize( double v, double min, double max );

	int *greedyScatteredSubsetSelection( double **points, int number_of_points, int number_of_dimensions, int number_to_select );

	extern std::mt19937 rng;
	extern long long random_seed;

	}
}

