#include "gomea/src/utils/linalg.hpp"

namespace gomea{
namespace utils{

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
matE choleskyDecomposition( const matE &matrix )
{
    int     i, j, k, info, *ipvt;
    double *a, *work;

	int n = matrix.rows();
    a    = (double *) Malloc( n*n*sizeof( double ) );
    work = (double *) Malloc( n*sizeof( double ) );
    ipvt = (int *) Malloc( n*sizeof( int ) );

    k = 0;
    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            a[k] = matrix(i,j);
            k++;
        }
        ipvt[i] = 0;
    }

    info = linpackDCHDC( a, n, n, work, ipvt );

    matE result = matE(n,n);
    if( info != n ) /* Matrix is not positive definite */
    {
        k = 0;
        for( i = 0; i < n; i++ )
        {
            for( j = 0; j < n; j++ )
            {
                result(i,j) = i != j ? 0.0 : sqrt( matrix(i,j) );
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
                result(i,j) = i < j ? 0.0 : a[k];
                k++;
            }
        }
    }

    free( ipvt );
    free( work );
    free( a );

    return( result );
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

        if (tst1 < std::abs(d[l]) + std::abs(e[l]))
            tst1 = std::abs(d[l]) + std::abs(e[l]);
        m = l;
        while (m < n) {
            if (std::abs(e[m]) <= eps*tst1) {
                /* if (std::abs(e[m]) + std::abs(d[m]+d[m+1]) == std::abs(d[m]+d[m+1])) { */
                break;
            }
            m++;
        }

        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */

        if (m > l) {
            int iter = 0;
            do { /* while (std::abs(e[l]) > eps*tst1); */
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

            } while (std::abs(e[l]) > eps*tst1);
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
    if( std::abs(a) > std::abs(b) )
    {
        r = b/a;
        r = std::abs(a)*sqrt(1+r*r);
    }
    else if (b != 0)
    {
        r = a/b;
        r = std::abs(b)*sqrt(1+r*r);
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
            scale = scale + std::abs(d[k]);
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

// method for calculating the pseudo-Inverse as recommended by Eigen developers
// Source: https://gist.github.com/pshriwise/67c2ae78e5db3831da38390a8b2a209f
matE pinv(const matE &a, double epsilon)
{
	Eigen::JacobiSVD<matE> svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
        // For a non-square matrix
        // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
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

}}