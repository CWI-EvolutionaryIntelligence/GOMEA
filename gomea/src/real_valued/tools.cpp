/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/real_valued/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

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

/**
 * BLAS subroutine.
 */
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy )
{
    double dtmp;

    if (n > 0)
    {
        incx *= sizeof( double );
        incy *= sizeof( double );

        dtmp  = (*dx);
        *dx   = (*dy);
        *dy   = dtmp;

        while( (--n) > 0 )
        {
            dx = (double *) ((char *) dx + incx);
            dy = (double *) ((char *) dy + incy);
            dtmp = (*dx); *dx = (*dy); *dy = dtmp;
        }
    }

    return( 0 );
}

/**
 * BLAS subroutine.
 */
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy)
{
    double dtmp0, dtmp, *dx0, *dy0;

    if( n > 0 && da != 0. )
    {
        incx *= sizeof(double);
        incy *= sizeof(double);
        *dy  += da * (*dx);

        if( (n & 1) == 0 )
        {
            dx   = (double *) ((char *) dx + incx);
            dy   = (double *) ((char *) dy + incy);
            *dy += da * (*dx);
            --n;
        }
        n = n >> 1;
        while( n > 0 )
        {
            dy0   = (double *) ((char *) dy + incy);
            dy    = (double *) ((char *) dy0 + incy);
            dtmp0 = (*dy0);
            dtmp  = (*dy);
            dx0   = (double *) ((char *) dx + incx);
            dx    = (double *) ((char *) dx0 + incx);
            *dy0  = dtmp0 + da * (*dx0);
            *dy   = dtmp + da * (*dx);
            --n;
        }
    }

    return( 0 );
}

/**
 * BLAS subroutine.
 */
void blasDSCAL( int n, double sa, double x[], int incx )
{
    int i, ix, m;

    if( n <= 0 )
    {
    }
    else if( incx == 1 )
    {
        m = n % 5;

        for( i = 0; i < m; i++ )
        {
            x[i] = sa * x[i];
        }

        for( i = m; i < n; i = i + 5 )
        {
            x[i]   = sa * x[i];
            x[i+1] = sa * x[i+1];
            x[i+2] = sa * x[i+2];
            x[i+3] = sa * x[i+3];
            x[i+4] = sa * x[i+4];
        }
    }
    else
    {
        if( 0 <= incx )
        {
            ix = 0;
        }
        else
        {
            ix = ( - n + 1 ) * incx;
        }

        for( i = 0; i < n; i++ )
        {
            x[ix] = sa * x[ix];
            ix = ix + incx;
        }
    }
}

/**
 * LINPACK subroutine.
 */
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] )
{
    int    info, j, jp, k, l, maxl, pl, pu;
    double maxdia, temp;

    pl   = 1;
    pu   = 0;
    info = p;
    for( k = 1; k <= p; k++ )
    {
        maxdia = a[k-1+(k-1)*lda];
        maxl   = k;
        if( pl <= k && k < pu )
        {
            for( l = k+1; l <= pu; l++ )
            {
                if( maxdia < a[l-1+(l-1)*lda] )
                {
                    maxdia = a[l-1+(l-1)*lda];
                    maxl   = l;
                }
            }
        }

        if( maxdia <= 0.0 )
        {
            info = k - 1;

            return( info );
        }

        if( k != maxl )
        {
            blasDSWAP( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );

            a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
            a[k-1+(k-1)*lda]       = maxdia;
            jp                     = ipvt[maxl-1];
            ipvt[maxl-1]           = ipvt[k-1];
            ipvt[k-1]              = jp;
        }
        work[k-1]        = sqrt( a[k-1+(k-1)*lda] );
        a[k-1+(k-1)*lda] = work[k-1];

        for( j = k+1; j <= p; j++ )
        {
            if( k != maxl )
            {
                if( j < maxl )
                {
                    temp                = a[k-1+(j-1)*lda];
                    a[k-1+(j-1)*lda]    = a[j-1+(maxl-1)*lda];
                    a[j-1+(maxl-1)*lda] = temp;
                }
                else if ( maxl < j )
                {
                    temp                = a[k-1+(j-1)*lda];
                    a[k-1+(j-1)*lda]    = a[maxl-1+(j-1)*lda];
                    a[maxl-1+(j-1)*lda] = temp;
                }
            }
            a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
            work[j-1]        = a[k-1+(j-1)*lda];
            temp             = -a[k-1+(j-1)*lda];

            blasDAXPY( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
        }
    }

    return( info );
}

/**
 * Computes the lower-triangle Cholesky Decomposition
 * of a square, symmetric and positive-definite matrix.
 * Subroutines from LINPACK and BLAS are used.
 */
double **choleskyDecomposition( double **matrix, int n )
{
    int     i, j, k, info, *ipvt;
    double *a, *work, **result;

    a    = (double *) Malloc( n*n*sizeof( double ) );
    work = (double *) Malloc( n*sizeof( double ) );
    ipvt = (int *) Malloc( n*sizeof( int ) );

    k = 0;
    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            a[k] = matrix[i][j];
            k++;
        }
        ipvt[i] = 0;
    }

    info = linpackDCHDC( a, n, n, work, ipvt );

    result = matrixNew( n, n );
    if( info != n ) /* Matrix is not positive definite */
    {
        k = 0;
        for( i = 0; i < n; i++ )
        {
            for( j = 0; j < n; j++ )
            {
                result[i][j] = i != j ? 0.0 : sqrt( matrix[i][j] );
                k++;
            }
        }
    }
    else
    {
        k = 0;
        for( i = 0; i < n; i++ )
        {
            for( j = 0; j < n; j++ )
            {
                result[i][j] = i < j ? 0.0 : a[k];
                k++;
            }
        }
    }

    free( ipvt );
    free( work );
    free( a );

    return( result );
}

/*mat choleskyDecomposition( mat matrix )
{
	mat result;
	if( !chol(result, matrix, "lower") )
	{
		printf("Warning: cholesky decomposition failed. Diag = [ ");
		//printf("Warning: cholesky decomposition failed.\n");
		//if( matrix(0,0) < 0 )
			//printf("AAAA\n");
		result = mat( matrix.n_cols, matrix.n_cols, fill::zeros );
		for( int i = 0; i < matrix.n_cols; i++ )
		{
			printf("%6.3e ",matrix(i,i));
			result(i,i) = sqrt(matrix(i,i));
		}
		printf("]\n");
	}
	return( result );
}*/

/**
 * LINPACK subroutine.
 */
int linpackDTRDI( double t[], int ldt, int n )
{
    int    j, k, info;
    double temp;

    info = 0;
    for( k = n; 1 <= k; k-- )
    {
        if ( t[k-1+(k-1)*ldt] == 0.0 )
        {
            info = k;
            break;
        }

        t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
        temp = -t[k-1+(k-1)*ldt];

        if ( k != n )
        {
            blasDSCAL( n-k, temp, t+k+(k-1)*ldt, 1 );
        }

        for( j = 1; j <= k-1; j++ )
        {
            temp = t[k-1+(j-1)*ldt];
            t[k-1+(j-1)*ldt] = 0.0;
            blasDAXPY( n-k+1, temp, t+k-1+(k-1)*ldt, 1, t+k-1+(j-1)*ldt, 1 );
        }
    }

    return( info );
}

/**
 * Computes the inverse of a matrix that is of
 * lower triangular form.
 */
double **matrixLowerTriangularInverse( double **matrix, int n )
{
    int     i, j, k, info;
    double *t, **result;

    t = (double *) Malloc( n*n*sizeof( double ) );

    k = 0;
    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            t[k] = matrix[j][i];
            k++;
        }
    }

    info = linpackDTRDI( t, n, n );

    result = matrixNew( n, n );
    k = 0;
    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            result[j][i] = i > j ? 0.0 : t[k];
            k++;
        }
    }

    free( t );

    return( result );
}

void eigenDecomposition( double **matrix, int n, double **D, double **Q )
{
    int     i, j;
    double *rgtmp, *diag;

    rgtmp = (double *) Malloc( n*sizeof( double ) );
    diag  = (double *) Malloc( n*sizeof( double ) );

    for( i = 0; i < n; i++ )
    {
        for( j = 0; j <= i; j++ )
        {
            Q[j][i] = matrix[j][i];
            Q[i][j] = Q[j][i];
        }
    }

    eigenDecompositionHouseholder2( n, Q, diag, rgtmp );
    eigenDecompositionQLalgo2( n, Q, diag, rgtmp );

    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            D[i][j] = 0.0;
        }
        D[i][i] = diag[i];
    }

    free( diag );
    free( rgtmp );
}


void eigenDecompositionQLalgo2( int n, double **V, double *d, double *e )
{
    int i, k, l, m;
    double f = 0.0;
    double tst1 = 0.0;
    double eps = 2.22e-16; /* Math.pow(2.0,-52.0);  == 2.22e-16 */

    /* shift input e */
    for (i = 1; i < n; i++) {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0; /* never changed again */

    for (l = 0; l < n; l++) {

        /* Find small subdiagonal element */

        if (tst1 < fabs(d[l]) + fabs(e[l]))
            tst1 = fabs(d[l]) + fabs(e[l]);
        m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) {
                /* if (fabs(e[m]) + fabs(d[m]+d[m+1]) == fabs(d[m]+d[m+1])) { */
                break;
            }
            m++;
        }

        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */

        if (m > l) {
            int iter = 0;
            do { /* while (fabs(e[l]) > eps*tst1); */
                double dl1, h;
                double g = d[l];
                double p = (d[l+1] - g) / (2.0 * e[l]);
                double r = myhypot(p, 1.);

                iter = iter + 1;  /* Could check iteration count here */

                /* Compute implicit shift */

                if (p < 0) {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                dl1 = d[l+1];
                h = g - d[l];
                for (i = l+2; i < n; i++) {
                    d[i] -= h;
                }
                f = f + h;

                /* Implicit QL transformation. */

                p = d[m];
                {
                    double c = 1.0;
                    double c2 = c;
                    double c3 = c;
                    double el1 = e[l+1];
                    double s = 0.0;
                    double s2 = 0.0;
                    for (i = m-1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = myhypot(p, e[i]);
                        e[i+1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i+1] = h + s * (c * g + s * d[i]);

                        /* Accumulate transformation. */

                        for (k = 0; k < n; k++) {
                            h = V[k][i+1];
                            V[k][i+1] = s * V[k][i] + c * h;
                            V[k][i] = c * V[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;
                }

                /* Check for convergence. */

            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    /* Sort eigenvalues and corresponding vectors. */
#if 1
    /* TODO: really needed here? So far not, but practical and only O(n^2) */
    {
        int j;
        double p;
        for (i = 0; i < n-1; i++) {
            k = i;
            p = d[i];
            for (j = i+1; j < n; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (j = 0; j < n; j++) {
                    p = V[j][i];
                    V[j][i] = V[j][k];
                    V[j][k] = p;
                }
            }
        }
    }
#endif
} /* QLalgo2 */

double myhypot( double a, double b )
{
    double r = 0;
    if( fabs(a) > fabs(b) )
    {
        r = b/a;
        r = fabs(a)*sqrt(1+r*r);
    }
    else if (b != 0)
    {
        r = a/b;
        r = fabs(b)*sqrt(1+r*r);
    }

    return r;
}

void eigenDecompositionHouseholder2( int n, double **V, double *d, double *e )
{
    int i,j,k;

    for (j = 0; j < n; j++) {
        d[j] = V[n-1][j];
    }

    /* Householder reduction to tridiagonal form */

    for (i = n-1; i > 0; i--) {

        /* Scale to avoid under/overflow */

        double scale = 0.0;
        double h = 0.0;
        for (k = 0; k < i; k++) {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 0; j < i; j++) {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        } else {

            /* Generate Householder vector */

            double f, g, hh;

            for (k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            f = d[i-1];
            g = sqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            /* Apply similarity transformation to remaining columns */

            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            hh = f / (h + h);
            for (j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for (j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (k = j; k <= i-1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    /* Accumulate transformations */

    for (i = 0; i < n-1; i++) {
        double h;
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++) {
                d[k] = V[k][i+1] / h;
            }
            for (j = 0; j <= i; j++) {
                double g = 0.0;
                for (k = 0; k <= i; k++) {
                    g += V[k][i+1] * V[k][j];
                }
                for (k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for (k = 0; k <= i; k++) {
            V[k][i+1] = 0.0;
        }
    }
    for (j = 0; j < n; j++) {
        d[j] = V[n-1][j];
        V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;

}

/**
 * Writes the contents of a matrix of dimensions n0 x n1 to a file.
 */
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 )
{
    int  i, j;
    char line_for_output[10000];

    sprintf( line_for_output, "[" );
    fputs( line_for_output, file );
    for( i = 0; i < n0; i++ )
    {
        sprintf( line_for_output, "[" );
        fputs( line_for_output, file );
        for( j = 0; j < n1; j++ )
        {
            sprintf( line_for_output, "%lf", matrix[i][j] );
            fputs( line_for_output, file );
            if( j < n1-1 )
            {
                sprintf( line_for_output, ", " );
                fputs( line_for_output, file );
            }
        }
        if( i == n0-1 )
            sprintf( line_for_output, "]" );
        else
            sprintf( line_for_output, "];" );
        fputs( line_for_output, file );
    }
    sprintf( line_for_output, "]\n" );
    fputs( line_for_output, file );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

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

int *getRanks( double *array, int array_size )
{
    int i, *sorted, *ranks;

    sorted = mergeSort( array, array_size );
    ranks = (int *) Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ ) ranks[sorted[i]] = i;

    free( sorted );
    return( ranks );
}

int *getRanksFromSorted( int *sorted, int array_size )
{
    int i, *ranks;

    ranks = (int *) Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ ) ranks[sorted[i]] = i;

    return( ranks );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

vec random1DNormalUnitVector( int length )
{
    vec result = vec(length);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for( int i = 0; i < length; i++ )
        result(i) = distribution(utils::rng);
    return result;
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

double min( double x, double y )
{
    if( x <= y )
        return x;
    return y;
}

double max( double x, double y )
{
    if( x >= y )
        return x;
    return y;
}

// method for calculating the pseudo-Inverse as recommended by Eigen developers
// Source: https://gist.github.com/pshriwise/67c2ae78e5db3831da38390a8b2a209f
mat pinv(const mat &a, double epsilon)
{
	Eigen::JacobiSVD<mat> svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
        // For a non-square matrix
        // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}


/*
 * Bhattacharyya distance
 */
/*double normalDistributionDistance( vec mu1, vec mu2, mat cov1, mat cov2 )
{
	vec mudiff = mu1 - mu2;
	mat covavg = (cov1 + cov2)/2.0;
	//double det1 = det(cov1);
	//double det2 = det(cov2);
	double det1 = cov1.determinant();
	double det2 = cov2.determinant();
	double detavg = covavg.determinant();

	mat a = mudiff.transpose() * pinv(covavg) * mudiff; 
	double dist = a(0)/8.0 + 0.5*log(detavg/sqrt(det1*det2));
	return( dist );
}*/

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
