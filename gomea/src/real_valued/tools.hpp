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
#include <limits>
#include <Eigen>

#include "gomea/src/utils/tools.hpp"

namespace gomea{
namespace realvalued{

using mat = Eigen::MatrixXd;
using vec = Eigen::VectorXd;

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
mat pinv(const mat &a, double epsilon = std::numeric_limits<double>::epsilon());
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
