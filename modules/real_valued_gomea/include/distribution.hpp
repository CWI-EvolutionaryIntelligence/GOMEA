/**
 *
 * RV-GOMEA
 *
 * If you use this software for any purpose, please cite the most recent publication:
 * A. Bouter, C. Witteveen, T. Alderliesten, P.A.N. Bosman. 2017.
 * Exploiting Linkage Information in Real-Valued Optimization with the Real-Valued
 * Gene-pool Optimal Mixing Evolutionary Algorithm. In Proceedings of the Genetic
 * and Evolutionary Computation Conference (GECCO 2017).
 * DOI: 10.1145/3071178.3071272
 *
 * Copyright (c) 1998-2017 Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jörn Grahl
 * - Anton Bouter
 *
 */

#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "tools.hpp"
#include "partial_solution.hpp"
#include "solution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class distribution_t {
	public:
		virtual ~distribution_t();

		// Parameter settings and default values
		double st_dev_ratio_threshold = 1.0;
		double distribution_multiplier_decrease = 0.9;
		double distribution_multiplier_increase = 1.0/0.9; // 1.0/c_dec

		// Variables
		double distribution_multiplier = 1.0;
		int samples_drawn = 0;
		int out_of_bounds_draws = 0;
			
		std::vector<int> variables;

		void adaptDistributionMultiplier( partial_solution_t** partial_solutions, int num_solutions );
		void adaptDistributionMultiplierMaximumStretch( partial_solution_t** partial_solutions, int num_solutions );
		virtual short generationalImprovementForOnePopulationForFOSElement( partial_solution_t** partial_solutions, int num_solutions, double *st_dev_ratio ) = 0;

		double estimateMean( int var, solution_t **selection, int selection_size );
		double estimateCovariance( int vara, int varb, solution_t **selection, int selection_size );
		vec estimateMeanVectorML( std::vector<int> variables, solution_t **selection, int selection_size );
		mat estimateRegularCovarianceMatrixML( std::vector<int> variables, vec mean_vector, solution_t **selection, int selection_size );
		mat estimateCovarianceMatrixML( std::vector<int> variables, solution_t **selection, int selection_size );
		mat estimateUnivariateCovarianceMatrixML( std::vector<int> variables, solution_t **selection, int selection_size );
		bool regularizeCovarianceMatrix( mat &cov_mat, vec &mean_vector, solution_t **selection, int selection_size );
		mat pseudoInverse( const mat &matrix );
		mat choleskyDecomposition( const mat &matrix );
		int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
		void blasDSCAL( int n, double sa, double x[], int incx );
		int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
		int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
			
		virtual void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, int visited[] );
		virtual void setOrder( const std::vector<int> &order ); 
		virtual void estimateDistribution( solution_t **selection, int selection_size ) = 0;	
		virtual partial_solution_t *generatePartialSolution( solution_t *parent ) = 0;
		virtual void print();
};

class normal_distribution_t : public distribution_t {
		public:
			normal_distribution_t( std::vector<int> variables );

			vec mean_vector;
			mat covariance_matrix;
			mat cholesky_decomposition; 
			
			void estimateDistribution( solution_t **selection, int selection_size );
			partial_solution_t *generatePartialSolution( solution_t *parent = NULL );
			
			short generationalImprovementForOnePopulationForFOSElement( partial_solution_t** partial_solutions, int num_solutions, double *st_dev_ratio );
};

class conditional_distribution_t : public distribution_t {
		public:
			conditional_distribution_t(); 
			conditional_distribution_t( const std::vector<int> &variables, const std::vector<int> &conditioned_variables );
			conditional_distribution_t( const std::vector<int> &variables, const std::set<int> &conditioned_variables );

			std::vector<int> order;
			std::vector<std::vector<int>> variable_groups;
			std::vector<std::vector<int>> neighboring_groups;
			std::vector<std::vector<int>> variables_conditioned_on;
			std::vector<std::vector<int>> index_in_var_array; // variables[index_in_var_array[a][b]] == variable_groups[a][b]

			std::vector<vec> mean_vectors;
			std::vector<vec> mean_vectors_conditioned_on;
			std::vector<mat> covariance_matrices;
			std::vector<mat> rho_matrices;
			std::vector<mat> cholesky_decompositions;

			void addGroupOfVariables( const std::vector<int> &indices, const std::set<int> &indices_cond );
			void addGroupOfVariables( std::vector<int> indices, std::vector<int> indices_cond );
			void addGroupOfVariables( int index, const std::vector<int> &indices_cond );
			void addGroupOfVariables( int index, int index_cond );
			void estimateDistribution( solution_t **selection, int selection_size );

			void setOrder( const std::vector<int> &order );
			void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, int visited[] );

			partial_solution_t *generatePartialSolution( solution_t *solution_conditioned_on = NULL ); 
			short generationalImprovementForOnePopulationForFOSElement( partial_solution_t** partial_solutions, int num_solutions, double *st_dev_ratio );
		private:
			void initializeMemory();
			void estimateConditionalGaussianML( int variable_group_index, solution_t **selection, int selection_size );
			void print();
};

}}
