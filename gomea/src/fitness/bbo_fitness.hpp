#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/embed.hpp"
#include "gomea/src/utils/tools.hpp"

namespace gomea{
namespace fitness{

template<class T>
class BBOFitnessFunction_t : public fitness_t<T> 
{
	public:
		BBOFitnessFunction_t( int number_of_parameters );
		BBOFitnessFunction_t( int number_of_parameters, double vtr );

		void initialize();

		double objectiveFunction( int objective_index, solution_t<T> *solution );
		virtual double objectiveFunction( int objective_index, vec_t<T> &variables ) = 0; 
		
		double constraintFunction( solution_t<T> *solution );
		virtual double constraintFunction( vec_t<T> &variables );

	private:
		void evaluationFunction( solution_t<T> *solution );
		void partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution );
};

}}
