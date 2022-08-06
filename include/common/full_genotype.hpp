#pragma once

#include "common/genotype.hpp"
#include "common/gomea_defs.hpp"

namespace gomea{

template<class T>
class full_genotype_t : public genotype_t<T>
{
	public:
		full_genotype_t();
		full_genotype_t(int number_of_variables);
		full_genotype_t( vec_t<T> &variables );
		
		T& operator[](std::size_t idx);
  		const T& operator[](std::size_t idx) const;

		std::size_t size();
		vec_t<T> &asVector();

	private:
		vec_t<T> variables;
		int identifier = 0; // TODO : use identifier to check if parent has not changed
};

template class full_genotype_t<int>;
template class full_genotype_t<float>;
template class full_genotype_t<double>;

}
