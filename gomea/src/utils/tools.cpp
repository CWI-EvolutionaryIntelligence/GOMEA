#include "gomea/src/utils/tools.hpp"

namespace gomea{
	namespace utils{

long long random_seed = 0;
std::mt19937 rng;

/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *Malloc( long size )
{
    void *result;

    result = (void *) malloc( size );
    if( !result )
    {
        printf("\n");
        printf("Error while allocating memory in Malloc( %ld ), aborting program.", size);
        printf("\n");

        exit( 0 );
    }

    return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Matrix -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Creates a new matrix with dimensions n x m.
 */
double **matrixNew( int n, int m )
{
    int      i;
    double **result;

    result = (double **) malloc( n*( sizeof( double * ) ) );
    for( i = 0; i < n; i++ )
        result[i] = (double *) malloc( m*( sizeof( double ) ) );

    return( result );
}

/**
 * Computes the dot product of two vectors of the same dimensionality n0.
 */
double vectorDotProduct( double *vector0, double *vector1, int n0 )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i < n0; i++ )
        result += vector0[i]*vector1[i];
    
    return( result );
}

/**
 * Computes the Euclidean norm of a given vector.
 */
double vectorNorm( double *vector0, int n0 )
{
    return( sqrt(vectorDotProduct( vector0, vector0, n0 )) );
}

/**
 * Computes the multiplication Av of a matrix A and a vector v
 * where matrix A has dimensions n0 x n1 and vector v has
 * dimensionality n1.
 */
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 )
{
    int     i;
    double *result;

    result = (double *) malloc( n0*sizeof( double ) );
    for( i = 0; i < n0; i++ )
        result[i] = vectorDotProduct( matrix[i], vector, n1 );

    return( result );
}

double *matrixVectorPartialMultiplication( double **matrix, double *vector, int n0, int n1, int number_of_elements, int *element_indices )
{
    int i,j,index;
    double *result;

    result = (double *) malloc( n0*sizeof( double ) );
    for( i = 0; i < n0; i++ )
        result[i] = 0;

    for( j = 0; j < number_of_elements; j++)
    {
        index = element_indices[j];
        for( i = 0; i < n0; i++ )
            result[i] += ( vector[index] * matrix[i][index] );
    }

    return result;
}
/**
 * Computes the matrix multiplication of two matrices A and B
 * of dimensions A: n0 x n1 and B: n1 x n2.
 */
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 )
{
    int     i, j, k;
    double **result;

    result = (double **) malloc( n0*sizeof( double * ) );
    for( i = 0; i < n0; i++ )
        result[i] = (double *) malloc( n2*sizeof( double ) );

    for( i = 0; i < n0; i++ )
    {
        for( j = 0; j < n2; j++ )
        {
            result[i][j] = 0;
            for( k = 0; k < n1; k++ )
                result[i][j] += matrix0[i][k]*matrix1[k][j];
        }
    }

    return( result );
}

vec_t<int> getSortedOrder( vec_t<int> &data ){
    vec_t<int> order(data.size());
    iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&data](size_t i1, size_t i2) {return data[i1] < data[i2];});
    return order;
}

/**
 * Computes the distance between two solutions a and b as
 * the Euclidean distance in parameter space.
 */
double distanceEuclidean( double *x, double *y, int number_of_dimensions )
{
    int    i;
    double value, result;

    result = 0.0;
    for( i = 0; i < number_of_dimensions; i++ )
    {
        value   = y[i] - x[i];
        result += value*value;
    }
    result = sqrt( result );

    return( result );
}

double distanceEuclidean( vec_t<double> &x, vec_t<double> &y ) 
{
    assert( x.size() == y.size() );
    return( distanceEuclidean(x.data(), y.data(), x.size()) );
}

/**
 * Computes the Euclidean distance between two points.
 */
double distanceEuclidean2D( double x1, double y1, double x2, double y2 )
{
    double result;

    result = (y1 - y2)*(y1-y2) + (x1-x2)*(x1-x2);
    result = sqrt( result );

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

/*vec random1DNormalUnitVector( int length )
{
    std::uniform_int_distribution<int> distribution(0,max);
	return randn<vec>(length);
}*/

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

/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of doubles and returns the sort-order (small to large).
 */
int *mergeSort( double *array, int array_size )
{
    int i, *sorted, *tosort;

    sorted = (int *) Malloc( array_size * sizeof( int ) );
    tosort = (int *) Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ )
        tosort[i] = i;

    if( array_size == 1 )
        sorted[0] = 0;
    else
        mergeSortWithinBounds( array, sorted, tosort, 0, array_size-1 );

    free( tosort );

    return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
    int r;

    if( p < q )
    {
        r = (p + q) / 2;
        mergeSortWithinBounds( array, sorted, tosort, p, r );
        mergeSortWithinBounds( array, sorted, tosort, r+1, q );
        mergeSortMerge( array, sorted, tosort, p, r+1, q );
    }
}
void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q )
{
    int r;

    if( p < q )
    {
        r = (p + q) / 2;
        mergeSortWithinBoundsInt( array, sorted, tosort, p, r );
        mergeSortWithinBoundsInt( array, sorted, tosort, r+1, q );
        mergeSortMergeInt( array, sorted, tosort, p, r+1, q );
    }
}
/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
    int i, j, k, first;

    i = p;
    j = r;
    for( k = p; k <= q; k++ )
    {
        first = 0;
        if( j <= q )
        {
            if( i < r )
            {
                if( array[tosort[i]] < array[tosort[j]] )
                    first = 1;
            }
        }
        else
            first = 1;

        if( first )
        {
            sorted[k] = tosort[i];
            i++;
        }
        else
        {
            sorted[k] = tosort[j];
            j++;
        }
    }

    for( k = p; k <= q; k++ )
        tosort[k] = sorted[k];
}

int *mergeSortInt( int *array, int array_size )
{
    int i, *sorted, *tosort;

    sorted = (int *) Malloc( array_size * sizeof( int ) );
    tosort = (int *) Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ )
        tosort[i] = i;

    if( array_size == 1 )
        sorted[0] = 0;
    else
        mergeSortWithinBoundsInt( array, sorted, tosort, 0, array_size-1 );

    free( tosort );

    return( sorted );
}


void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q )
{
    int i, j, k, first;

    i = p;
    j = r;
    for( k = p; k <= q; k++ )
    {
        first = 0;
        if( j <= q )
        {
            if( i < r )
            {
                if( array[tosort[i]] < array[tosort[j]] )
                    first = 1;
            }
        }
        else
            first = 1;

        if( first )
        {
            sorted[k] = tosort[i];
            i++;
        }
        else
        {
            sorted[k] = tosort[j];
            j++;
        }
    }

    for( k = p; k <= q; k++ )
        tosort[k] = sorted[k];
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


double normalize( double v, double min, double max )
{
	return( (v-min)/(max-min) );
}

int *greedyScatteredSubsetSelection( double **points, int number_of_points, int number_of_dimensions, int number_to_select )
{
  int     i, index_of_farthest, random_dimension_index, number_selected_so_far,
         *indices_left, *result;
  double *nn_distances, distance_of_farthest, value;

  if( number_to_select > number_of_points )
  {
    printf("\n");
    printf("Error: greedyScatteredSubsetSelection asked to select %d solutions from set of size %d.", number_to_select, number_of_points);
    printf("\n\n");

    exit( 0 );
  }

  result = (int *) Malloc( number_to_select*sizeof( int ) );

  indices_left = (int *) Malloc( number_of_points*sizeof( int ) );
  for( i = 0; i < number_of_points; i++ )
    indices_left[i] = i;

  /* Find the first point: maximum value in a randomly chosen dimension */
  random_dimension_index = utils::randomInt( number_of_dimensions );

  index_of_farthest    = 0;
  distance_of_farthest = points[indices_left[index_of_farthest]][random_dimension_index];
  for( i = 1; i < number_of_points; i++ )
  {
    if( points[indices_left[i]][random_dimension_index] > distance_of_farthest )
    {
      index_of_farthest    = i;
      distance_of_farthest = points[indices_left[i]][random_dimension_index];
    }
  }

  number_selected_so_far          = 0;
  result[number_selected_so_far]  = indices_left[index_of_farthest];
  indices_left[index_of_farthest] = indices_left[number_of_points-number_selected_so_far-1];
  number_selected_so_far++;

  /* Then select the rest of the solutions: maximum minimum
   * (i.e. nearest-neighbour) distance to so-far selected points */
  nn_distances = (double *) Malloc( number_of_points*sizeof( double ) );
  for( i = 0; i < number_of_points-number_selected_so_far; i++ )
    nn_distances[i] = utils::distanceEuclidean( points[indices_left[i]], points[result[number_selected_so_far-1]], number_of_dimensions );

  while( number_selected_so_far < number_to_select )
  {
    index_of_farthest    = 0;
    distance_of_farthest = nn_distances[0];
    for( i = 1; i < number_of_points-number_selected_so_far; i++ )
    {
      if( nn_distances[i] > distance_of_farthest )
      {
        index_of_farthest    = i;
        distance_of_farthest = nn_distances[i];
      }
    }

    result[number_selected_so_far]  = indices_left[index_of_farthest];
    indices_left[index_of_farthest] = indices_left[number_of_points-number_selected_so_far-1];
    nn_distances[index_of_farthest] = nn_distances[number_of_points-number_selected_so_far-1];
    number_selected_so_far++;

    for( i = 0; i < number_of_points-number_selected_so_far; i++ )
    {
      value = utils::distanceEuclidean( points[indices_left[i]], points[result[number_selected_so_far-1]], number_of_dimensions );
      if( value < nn_distances[i] )
        nn_distances[i] = value;
    }
  }

  free( nn_distances );
  free( indices_left );

  return( result );
}


}}
