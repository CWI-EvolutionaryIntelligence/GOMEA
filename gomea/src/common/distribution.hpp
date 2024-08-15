#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/utils/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/cond_factor.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/utils/rng.hpp"
#include "gomea/src/utils/linalg.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{

class distribution_t{
	public:
		virtual ~distribution_t(){};

	protected:
		int num_variables;
};

class distribution_Dt : public distribution_t{
	public:
		void initializeFrequencyTables( const vec_t<solution_t<char>*> &population );

	private:
		vec_t<vec_t<int>> frequency_tables;
};

}
