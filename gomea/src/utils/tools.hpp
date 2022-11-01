#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{
	namespace utils{

	double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
	double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
	double vectorDotProduct( double *vector0, double *vector1, int n0 );

	vec_t<int> getSortedOrder( vec_t<int> &data );
	//template vec_t<int> getSortedOrder( vec_t<float> &data );
	//template vec_t<int> getSortedOrder( vec_t<double> &data );
	
	}
}

