#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/common/distribution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class distribution_Rt : public distribution_t {
	public:
		distribution_Rt( vec_t<int> variables );
		distribution_Rt( const vec_t<int> &variables, const vec_t<int> &conditioned_variables ) : variables(variables), variables_conditioned_on(conditioned_variables){};
		distribution_Rt( const vec_t<int> &variables, const std::set<int> &conditioned_variables );
		virtual ~distribution_Rt(){};

		// Parameter settings and default values
		int samples_drawn = 0;
		int out_of_bounds_draws = 0;
			
		vec_t<int> variables;
		vec_t<int> variables_conditioned_on;

		vec_t<double> mean_vector;
		vec_t<double> mean_vector_conditioned_on;
		matE covariance_matrix;
		matE cholesky_decomposition; 
		matE rho_matrix;

		double getStandardDeviationRatio( partial_solution_t<double> **partial_solutions, int num_solutions );

		static matE estimateFullCovarianceMatrixML( solution_t<double> **selection, int selection_size, double distribution_multiplier );
		static double estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size );
		static double estimateMean( int var, solution_t<double> **selection, int selection_size );
		
		vec_t<double> estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		matE estimateRegularCovarianceMatrixML( vec_t<int> &variables, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size, double distribution_multiplier );
		matE estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size, double distribution_multiplier );
		matE estimateUnivariateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size, double distribution_multiplier );

		bool regularizeCovarianceMatrix( matE &cov_mat, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size );

		vecE sample();
		vecE sample(const vec_t<double> &sample_means);
		vec_t<double> getConditionalSampleMeans(vec_t<double> solution_conditioned_on);
		vec_t<double> getConditionalSampleMeans(solution_t<double> *solution_conditioned_on);

		void estimateDistribution( solution_t<double> **selection, int selection_size, double distribution_multiplier );

		void print();
};

}}
