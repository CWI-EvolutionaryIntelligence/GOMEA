#pragma once

#include "common/gomea_defs.hpp"
#include "common/genotype.hpp"
#include "common/full_genotype.hpp"

namespace gomea{

template<class T>
class partial_genotype_t : public genotype_t<T>
{
	public:
		partial_genotype_t( full_genotype_t<T> *parent_genotype, vec_t<int> touched_indices, vec_t<T> touched_variables );
		
		T& operator[](std::size_t idx);
  		const T& operator[](std::size_t idx) const;

	private:
		//bool checkParentIdentifierConsistency();

		full_genotype_t<T> *parent_genotype;
		vec_t<int> touched_indices;
		vec_t<T> touched_variables;
		//int parent_identifier = -1;
};

template class partial_genotype_t<int>;
template class partial_genotype_t<float>;
template class partial_genotype_t<double>;

}
