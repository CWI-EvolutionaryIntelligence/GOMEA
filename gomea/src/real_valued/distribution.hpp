#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/real_valued/partial_solutionRV.hpp"
#include "gomea/src/real_valued/solutionRV.hpp"
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
			
		vec_t<int> variables;

		void adaptDistributionMultiplier( partial_solution_t<double>** partial_solutions, int num_solutions );
		void adaptDistributionMultiplierMaximumStretch( partial_solution_t<double>** partial_solutions, int num_solutions );
		virtual bool generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio ) = 0;

		static mat estimateFullCovarianceMatrixML( solution_t<double> **selection, int selection_size );
		static double estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size );
		static double estimateMean( int var, solution_t<double> **selection, int selection_size );
		
		vec_t<double> estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		mat estimateRegularCovarianceMatrixML( vec_t<int> &variables, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size );
		mat estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		mat estimateUnivariateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		bool regularizeCovarianceMatrix( mat &cov_mat, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size );
		mat pseudoInverse( const mat &matrix );
		mat choleskyDecomposition( const mat &matrix );
		int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
		void blasDSCAL( int n, double sa, double x[], int incx );
		int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
		int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
			
		virtual void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited );
		virtual void setOrder( const vec_t<int> &order ); 
		virtual void estimateDistribution( solution_t<double> **selection, int selection_size ) = 0;	
		virtual partial_solution_t<double> *generatePartialSolution( solution_t<double> *parent, fitness::fitness_generic_t *fitness_function = NULL ) = 0;
		virtual void print();
};

class normal_distribution_t : public distribution_t {
		public:
			normal_distribution_t( vec_t<int> variables );

			vec_t<double> mean_vector;
			mat covariance_matrix;
			mat cholesky_decomposition; 
			
			void estimateDistribution( solution_t<double> **selection, int selection_size );
			partial_solution_t<double> *generatePartialSolution( solution_t<double> *parent = NULL, fitness::fitness_generic_t *fitness_function = NULL );
			
			bool generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio );
};

class conditional_distribution_t : public distribution_t {
		public:
			conditional_distribution_t(); 
			conditional_distribution_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables );
			conditional_distribution_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables );

			vec_t<int> order;
			vec_t<vec_t<int>> variable_groups;
			vec_t<vec_t<int>> neighboring_groups;
			vec_t<vec_t<int>> variables_conditioned_on;
			vec_t<vec_t<int>> index_in_var_array; // variables[index_in_var_array[a][b]] == variable_groups[a][b]

			vec_t<vec_t<double>> mean_vectors;
			vec_t<vec_t<double>> mean_vectors_conditioned_on;
			vec_t<mat> covariance_matrices;
			vec_t<mat> rho_matrices;
			vec_t<mat> cholesky_decompositions;

			void addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond );
			void addGroupOfVariables( vec_t<int> indices, vec_t<int> indices_cond );
			void addGroupOfVariables( int index, const vec_t<int> &indices_cond );
			void addGroupOfVariables( int index, int index_cond );
			void estimateDistribution( solution_t<double> **selection, int selection_size );

			void setOrder( const vec_t<int> &order );
			void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited );

			partial_solution_t<double> *generatePartialSolution( solution_t<double> *solution_conditioned_on = NULL, fitness::fitness_generic_t *fitness_function = NULL ); 
			bool generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio );
		private:
			void initializeMemory();
			void estimateConditionalGaussianML( int variable_group_index, solution_t<double> **selection, int selection_size );
			void print();
};

}}
