#include "common/partial_genotype.hpp"

namespace gomea{

template<class T>
partial_genotype_t<T>::partial_genotype_t( full_genotype_t<T> *parent_genotype, vec_t<int> touched_indices, vec_t<T> touched_variables )
{
	this->parent_genotype = parent_genotype;
	this->touched_indices = touched_indices;
	this->touched_variables = touched_variables;
}

template<class T>
T& partial_genotype_t<T>::operator[](std::size_t idx) {
	return touched_variables[idx];
}

template<class T>
const T& partial_genotype_t<T>::operator[](std::size_t idx) const {
	return touched_variables[idx];
}


}
