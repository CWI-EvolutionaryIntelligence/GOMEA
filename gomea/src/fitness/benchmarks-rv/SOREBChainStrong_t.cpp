#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

SOREBChainStrong_t::SOREBChainStrong_t( int number_of_variables, double vtr ) : GBOFitnessFunction_t(number_of_variables,vtr)
{
	this->name = "Sum of rotated ellipsoid blocks (SOREB) with strong dependencies -- 1 overlap in chain structure";
	this->initialize();
	this->rotation_angle = -45;
	this->rotation_block_size = 2;
	this->rotation_matrix = initializeObjectiveRotationMatrix( this->rotation_angle, this->rotation_block_size );
}

SOREBChainStrong_t::~SOREBChainStrong_t()
{
	for( int i = 0; i < rotation_block_size; i++ )
		delete[] rotation_matrix[i];
	delete[] rotation_matrix;
}

int SOREBChainStrong_t::getNumberOfSubfunctions()
{
	return number_of_variables-1;
}

vec_t<int> SOREBChainStrong_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec(2);
	vec.push_back(subfunction_index);
	vec.push_back(subfunction_index+1);
	return vec;
}

double SOREBChainStrong_t::subfunction( int subfunction_index, vec_t<double> &variables )
{
	double vars[] = {variables[subfunction_index],variables[subfunction_index+1]};
	double *rotated_variables = rotateVariables(&vars[0],rotation_block_size,this->rotation_matrix);
	double result = 0.0;
	for( int i = 0; i < rotation_block_size; i++)
	{
		result += pow( 10.0, 6.0*(((double) (i))/((double) (rotation_block_size-1))) )*rotated_variables[i]*rotated_variables[i];
	}
	delete[] rotated_variables;
	return result;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
