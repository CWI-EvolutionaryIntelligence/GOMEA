/**
 *
 * RV-GOMEA
 *
 * If you use this software for any purpose, please cite the most recent publication:
 * A. Bouter, C. Witteveen, T. Alderliesten, P.A.N. Bosman. 2017.
 * Exploiting Linkage Information in Real-Valued Optimization with the Real-Valued
 * Gene-pool Optimal Mixing Evolutionary Algorithm. In Proceedings of the Genetic
 * and Evolutionary Computation Conference (GECCO 2017).
 * DOI: 10.1145/3071178.3071272
 *
 * Copyright (c) 1998-2017 Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jörn Grahl
 * - Anton Bouter
 *
 */

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <map>
#include <memory>
#include <cassert>
#include <armadillo>

#include "gomea/src/utils/tools.hpp"

namespace gomea{
namespace realvalued{

using namespace arma;

template<class T>
using vec_t = std::vector<T>;
template<class T>
using vec_pt = std::shared_ptr<vec_t<T>>;

void *Malloc( long size );
double **matrixNew( int n, int m );
double vectorDotProduct( double *vector0, double *vector1, int n0 );
double vectorNorm( double *vector0, int n0 );
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
void blasDSCAL( int n, double sa, double x[], int incx );
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
double **choleskyDecomposition( double **matrix, int n );
mat choleskyDecomposition( mat matrix );
int linpackDTRDI( double t[], int ldt, int n );
double **matrixLowerTriangularInverse( double **matrix, int n );
void eigenDecomposition( double **matrix, int n, double **D, double **Q );
void eigenDecompositionQLalgo2( int n, double **V, double *d, double *e );
double myhypot( double a, double b );
void eigenDecompositionHouseholder2( int n, double **V, double *d, double *e );
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 );

int *mergeSort( double *array, int array_size );
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q );

void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortInt( int *array, int array_size );
void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q );

int *getRanks(double *array, int array_size );
int *getRanksFromSorted(int *sorted, int array_size );

vec random1DNormalUnitVector( int length );

int *greedyScatteredSubsetSelection( double **points, int number_of_points, int number_of_dimensions, int number_to_select );

double max( double x, double y );
double min( double x, double y );
double distanceEuclidean( double *solution_a, double *solution_b, int n );
double distanceEuclidean( vec_t<double> &x, vec_t<double> &y );
double distanceEuclidean2D( double x1, double y1, double x2, double y2 );
double normalize( double v, double min, double max );

double *matrixVectorPartialMultiplication( double **matrix, double *vector, int n0, int n1, int number_of_elements, int *element_indices );

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Constants -=-=-=-=-=-=-=-=-=-=-=-=-=-*/
#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798
#endif
#define FALSE 0
#define TRUE 1
#define OMP_NUM_THREADS 1
#define OPENBLAS_NUM_THREADS 1

}}
