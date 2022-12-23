#include "gomea/src/utils/tools.hpp"

namespace gomea{
	namespace utils{

long long random_seed = 0;
std::mt19937 rng;
    
vec_t<int> getSortedOrder( vec_t<int> &data ){
    vec_t<int> order(data.size());
    iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&data](size_t i1, size_t i2) {return data[i1] < data[i2];});
    return order;
}

/**
 * Computes the matrix multiplication of two matrices A and B
 * of dimensions A: n0 x n1 and B: n1 x n2.
 */
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 )
{
    double **result = new double*[n0];
    for( int i = 0; i < n0; i++ )
        result[i] = new double[n2];

    for( int i = 0; i < n0; i++ )
    {
        for( int j = 0; j < n2; j++ )
        {
            result[i][j] = 0;
            for( int k = 0; k < n1; k++ )
                result[i][j] += matrix0[i][k]*matrix1[k][j];
        }
    }
    return( result );
}

double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 )
{
    double *result = new double[n0];
    for( int i = 0; i < n0; i++ )
        result[i] = vectorDotProduct( matrix[i], vector, n1 );
    return( result );
}

double vectorDotProduct( double *vector0, double *vector1, int n0 )
{
    double result = 0.0;
    for( int i = 0; i < n0; i++ )
        result += vector0[i]*vector1[i];
    return( result );
}

bool isPowerOfK(int n, int k)
{
    double logNBaseK = log(n) / log(k);
    return (ceil(logNBaseK) == floor(logNBaseK));
}

vec_t<int> randomPermutation( int size )
{
    vec_t<int> perm(size);
    iota(perm.begin(), perm.end(), 0);
    std::shuffle( perm.begin(), perm.end(), rng );
    return( perm );
}

double randomRealUniform01()
{
    static std::uniform_real_distribution<double> distribution(0.0,1.0);
	return distribution(rng);
}

int randomInt( int max )
{
    std::uniform_int_distribution<int> distribution(0,max);
	return distribution(rng);
}

void initializeRandomNumberGenerator()
{
	utils::random_seed = static_cast<long long>(std::chrono::system_clock::now().time_since_epoch().count());
	rng.seed(utils::random_seed);
}

void initializeRandomNumberGenerator( long long seed )
{
    utils::random_seed = seed;
	rng.seed(utils::random_seed);
}

}}
