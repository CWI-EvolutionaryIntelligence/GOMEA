#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

using namespace gomea;

SOREBChainStrongBBO_t::SOREBChainStrongBBO_t( int number_of_variables, double vtr ) : BBOFitnessFunction_t(number_of_variables, vtr)
{
	this->name = "Sum of rotated ellipsoid blocks (SOREB) with strong dependencies -- 1 overlap in chain structure";
	this->initialize();
	this->rotation_angle = -45;
	this->rotation_block_size = 2;
	this->rotation_matrix = initializeObjectiveRotationMatrix( this->rotation_angle, this->rotation_block_size );
}

double SOREBChainStrongBBO_t::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double result = 0.0;
	for( int i = 0; i < number_of_variables-1; i++ )
	{
		double vars[] = {variables[i],variables[i+1]};
		double *rotated_variables = rotateVariables(&vars[0],rotation_block_size,this->rotation_matrix);
		for( int j = 0; j < rotation_block_size; j++)
		{
			result += pow( 10.0, 6.0*(((double) (j))/((double) (rotation_block_size-1))) )*rotated_variables[j]*rotated_variables[j];
		}
		delete[] rotated_variables;
	}
	return result;
}
		
}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
