#include "common/full_genotype.hpp"

namespace gomea{

template<class T>
full_genotype_t<T>::full_genotype_t()
{
	this->variables = vec_t<T>(1);
}

template<class T>
full_genotype_t<T>::full_genotype_t( int number_of_variables )
{
	this->variables = vec_t<T>(number_of_variables);
}

template<class T>
full_genotype_t<T>::full_genotype_t( vec_t<T> &variables )
{
	this->variables = variables;
}
		
template<class T>
T& full_genotype_t<T>::operator[](std::size_t idx) {
	return variables[idx];
}

template<class T>
const T& full_genotype_t<T>::operator[](std::size_t idx) const {
	return variables[idx];
}

template<class T>
std::size_t full_genotype_t<T>::size()
{
	return variables.size();
}

template<class T>
vec_t<T> &full_genotype_t<T>::asVector()
{
	return variables;
}

}
