#pragma once

#include <Eigen>

#include "gomea/src/utils/tools.hpp"

using matE = Eigen::MatrixXd;
using vecE = Eigen::VectorXd;

namespace gomea{
	namespace utils{
		int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
		int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
		void blasDSCAL( int n, double sa, double x[], int incx );
		int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );

		matE choleskyDecomposition( const matE &matrix );
		double **matrixLowerTriangularInverse( double **matrix, int n );
		void eigenDecomposition( double **matrix, int n, double **D, double **Q );
		void eigenDecompositionQLalgo2( int n, double **V, double *d, double *e );
		double myhypot( double a, double b );
		void eigenDecompositionHouseholder2( int n, double **V, double *d, double *e );
		matE pinv(const matE &a, double epsilon = std::numeric_limits<double>::epsilon());

		void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 );
	}
}