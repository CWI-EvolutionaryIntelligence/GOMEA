#include "gomea/src/utils/tools.hpp"
#include <queue>

namespace gomea{
namespace utils{


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

template<class T>
void reorder( vec_t<T> &vec, const vec_t<int> &_order )
{
    assert(vec.size() == _order.size());
    vec_t<int> order = _order;
    for (size_t i = 0; i < order.size(); i++)
        if (i < order[i])
        {
            int alt = order[i];
            std::swap( vec[i], vec[alt] );
            std::swap( order[i], order[alt] );
        }
}
template void reorder<int>( vec_t<int> &vec, const vec_t<int> &order );
template void reorder<char>( vec_t<char> &vec, const vec_t<int> &order );
template void reorder<double>( vec_t<double> &vec, const vec_t<int> &order );
template void reorder<float>( vec_t<float> &vec, const vec_t<int> &order );
template void reorder<vec_t<int>>( vec_t<vec_t<int>> &vec, const vec_t<int> &order );
template void reorder<vec_t<char>>( vec_t<vec_t<char>> &vec, const vec_t<int> &order );
template void reorder<vec_t<double>>( vec_t<vec_t<double>> &vec, const vec_t<int> &order );
template void reorder<vec_t<float>>( vec_t<vec_t<float>> &vec, const vec_t<int> &order );

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

std::vector<int> getGraphOrderBreadthFirst( const graph_t &graph )
{
	const int UNVISITED = 0;
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	const int IN_QUEUE = 3;

	int num_nodes = graph.size();
	std::vector<int> visited(num_nodes,0);
	vec_t<int> var_order = gomea::utils::randomPermutation( num_nodes );

	std::vector<int> graph_order;
	for( int i = 0; i < num_nodes; i++ )
	{
		int ind = var_order[i];
		if( visited[ind] == IS_VISITED )
			continue;
		visited[ind] = IN_CLIQUE;
	
		std::queue<int> q;
		q.push(ind);

		while( !q.empty() )
		{
			ind = q.front();
			q.pop();

			if( visited[ind] == IS_VISITED )
				continue;
			visited[ind] = IS_VISITED;

			graph_order.push_back(ind);

			for( int x : graph.at(ind) ) 
			{
				if( visited[x] == UNVISITED )
				{
					q.push(x);
					visited[x] = IN_QUEUE;
					//printf("Q[ %d ]\n",x);
				}
			}
		}
	}
	return( graph_order );
}

int *hungarianAlgorithm( int **similarity_matrix, int dim )
{
	int x,y,ty;

	int *lx = (int*) Malloc(dim*sizeof(int));
	int *ly = (int*) Malloc(dim*sizeof(int));
	int *xy = (int*) Malloc(dim*sizeof(int));
	int *yx = (int*) Malloc(dim*sizeof(int));
	int *slack = (int*) Malloc(dim*sizeof(int));
	int *slackx = (int*) Malloc(dim*sizeof(int));
	int *prev = (int*) Malloc(dim*sizeof(int));
	bool *S = (bool*) Malloc(dim*sizeof(bool));
	bool *T = (bool*) Malloc(dim*sizeof(bool));

	int root = -1;
	int max_match = 0;
	for(int i = 0; i < dim; i++ )
	{
		lx[i] = 0;
		ly[i] = 0;
		xy[i] = -1;
		yx[i] = -1;
	}
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			if(similarity_matrix[i][j] > lx[i])
				lx[i] = similarity_matrix[i][j];

	bool terminated = false;
	while(!terminated)
	{
		if (max_match == dim) break;

		int wr = 0;
		int rd = 0;
		int *q = (int*) Malloc(dim*sizeof(int));
		for(int i = 0; i < dim; i++ )
		{
			S[i] = false;
			T[i] = false;
			prev[i] = -1;
		}

		for (x = 0; x < dim; x++)
		{
			if (xy[x] == -1)
			{
				q[wr++] = root = x;
				prev[x] = -2;
				S[x] = true;
				break;
			}
		}

		for (y = 0; y < dim; y++)
		{
			slack[y] = lx[root] + ly[y] - similarity_matrix[root][y];
			slackx[y] = root;
		}

		while ( 1 )
		{
			while (rd < wr)
			{
				x = q[rd++];
				for (y = 0; y < dim; y++)
				{
					if (similarity_matrix[x][y] == lx[x] + ly[y] && !T[y])
					{
						if (yx[y] == -1) break;
						T[y] = true;
						q[wr++] = yx[y];
						hungarianAlgorithmAddToTree(yx[y], x, S, prev, slack, slackx, lx, ly, similarity_matrix, dim);
					}
				}
				if (y < dim) break;
			}
			if (y < dim) break;

			int delta = 100000000;
			for(y = 0; y < dim; y++)
				if(!T[y] && slack[y] < delta)
					delta = slack[y];
			for(x = 0; x < dim; x++)
				if(S[x])
					lx[x] -= delta;
			for(y = 0; y < dim; y++)
				if(T[y])
					ly[y] += delta;
			for(y = 0; y < dim; y++)
				if(!T[y])
					slack[y] -= delta;

			wr = 0;
			rd = 0;
			for (y = 0; y < dim; y++)
			{
				if (!T[y] && slack[y] == 0)
				{
					if (yx[y] == -1)
					{
						x = slackx[y];
						break;
					}
					else
					{
						T[y] = true;
						if (!S[yx[y]])
						{
							q[wr++] = yx[y];
							hungarianAlgorithmAddToTree(yx[y], slackx[y], S, prev, slack, slackx, lx, ly, similarity_matrix, dim);
						}
					}
				}
			}
			if (y < dim) break;
		}

		if (y < dim)
		{
			max_match++;
			for (int cx = x, cy = y; cx != -2; cx = prev[cx], cy = ty)
			{
				ty = xy[cx];
				yx[cy] = cx;
				xy[cx] = cy;
			}
		}
		else terminated = true;

		free( q );
	}

	free( lx );
	free( ly );
	free( yx );
	free( slack );
	free( slackx );
	free( prev );
	free( S );
	free( T );

	return xy;
}

void hungarianAlgorithmAddToTree(int x, int prevx, bool *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim) 
{
	S[x] = true;
	prev[x] = prevx;
	for (int y = 0; y < dim; y++)
	{
		if (lx[x] + ly[y] - similarity_matrix[x][y] < slack[y])
		{
			slack[y] = lx[x] + ly[y] - similarity_matrix[x][y];
			slackx[y] = x;
		}
	}
}


}}
