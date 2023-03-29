/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

oneMax_t::oneMax_t( int number_of_variables ) : GBOFitnessFunction_t(number_of_variables)
{
	this->name = "OneMax function";
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->initialize();
}
		
int oneMax_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> oneMax_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}
		
double oneMax_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	return( variables[subfunction_index] );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
