#include "utils/tools.hpp"

namespace gomea{
	namespace utils{
    
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

}}
