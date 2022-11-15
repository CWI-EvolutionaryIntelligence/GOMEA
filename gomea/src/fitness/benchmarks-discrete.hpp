#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class oneMax_t: public fitness_t<char>
{
	public:
		oneMax_t( int number_of_variables );
		
	private:
		void evaluationFunction( solution_t<char> *solution );
		void partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution );
		double subfunction( char x );
};

}}
