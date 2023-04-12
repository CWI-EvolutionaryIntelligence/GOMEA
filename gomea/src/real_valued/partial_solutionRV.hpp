#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "src/common/gomea_defs.hpp"
#include "src/common/solution.hpp"
#include "src/common/partial_solution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

template<class T>
class partial_solution_t : public gomea::partial_solution_t<T>
{
	using gomea::partial_solution_t<T>::partial_solution_t;
	
	public:
		partial_solution_t( vec_t<T> &touched_variables, vec_t<T> &sample_zs, vec_t<int> &touched_indices );

		vec_t<T> sample_zs; // Samples z~N(0,I), later transformed to N(mu,C)
		vec_t<T> sample_means;

		bool is_accepted = false;
		bool improves_elitist = false;

		void setSampleMean( vec_t<T> &means );
};

}}
