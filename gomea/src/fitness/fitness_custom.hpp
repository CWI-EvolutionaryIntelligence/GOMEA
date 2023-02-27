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
class customFitnessFunction_t : public fitness_t<T> 
{
	public:
		customFitnessFunction_t( int number_of_parameters );
		customFitnessFunction_t( int number_of_parameters, double vtr );
		
		void initializeVariableInteractionGraph(); 
		virtual double getSimilarityMeasure( size_t var_a, size_t var_b );

	private:
		void evaluationFunction( solution_t<T> *solution );
		void partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution );
		virtual double subfunction( int subfunction_index, vec_t<T> &variables ) = 0;

		double objectiveFunction( int objective_index, solution_t<T> *solution );
		double objectiveFunction( int objective_index, partial_solution_t<T> *solution );
		virtual double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers ); 
		
		double constraintFunction( solution_t<T> *solution );
		double constraintFunction( partial_solution_t<T> *solution );
		virtual double constraintFunction( vec_t<double> &fitness_buffers ); 
};

template class customFitnessFunction_t<char>;
template class customFitnessFunction_t<double>;

}}
