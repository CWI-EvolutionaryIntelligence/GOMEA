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
#include "distribution.hpp"
#include "partial_solution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

class fos_t {

		public:
			fos_t( int number_of_variables );
			fos_t( int number_of_variables, int FOS_element_size );
			fos_t( int problem_index, int number_of_variables, int FOS_type );
			fos_t( int problem_index, double **covariance_matrix, int n );
			fos_t( FILE *file );
			//fos_t( const std::map<int,std::set<int>> &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element, int VIG_order );
			fos_t( int number_of_variables, const std::map<int,std::set<int>> &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			fos_t( const fos_t &f );
			~fos_t();
			
			int getLength();
			std::vector<int> getSet( int element_index );
			int getSetLength( int element_index );
			double getAcceptanceRate(); 
			int getDistributionMultiplier( int element_index );
			
			void addGroup( int var_index );
			void addGroup( const std::set<int> &group );
			void addGroup( std::vector<int> group );
			void addGroup( distribution_t *dist );
			void addConditionedGroup( std::vector<int> variables );
			void addConditionedGroup( std::vector<int> variables, std::set<int> conditioned_variables );
			void randomizeOrder();
			void randomizeOrder( const std::map<int,std::set<int>> &variable_interaction_graph ); 
			std::vector<int> getVIGOrderBreadthFirst( const std::map<int,std::set<int>> &variable_interaction_graph );

			double getSimilarity( int a, int b, int *mpm_num_ind );
			double **computeMIMatrix( double **covariance_matrix, int n );
			void inheritDistributionMultipliers( fos_t *other, double *multipliers );
			int *matchFOSElements( fos_t *other );
			int *hungarianAlgorithm( int** similarity_matrix, int dim );
			void hungarianAlgorithmAddToTree(int x, int prevx, short *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim );
			int determineNearestNeighbour(int index, double **S_matrix, int *mpm_num_ind, int mpm_length );
			void ezilaitini();

			void initializeNormalDistribution(int FOS_index);
			void initializeConditionalDistribution( int FOS_index );

			partial_solution_t *generatePartialSolution( int FOS_index, solution_t *solution_conditioned_on );
			void estimateDistributions( solution_t **selection, int selection_size );
			void estimateDistribution( int FOS_index, solution_t **selection, int selection_size );
			void adaptDistributionMultiplier( int FOS_index, partial_solution_t **solutions, int num_solutions );

			int number_of_variables = -1;
			std::vector<distribution_t*> distributions;
			int no_improvement_stretch = 0;
			int maximum_no_improvement_stretch = 100;
			
			double p_accept = 0.00;
			std::vector<std::vector<int>> sets;
			std::vector<uvec> variables_conditioned_on; 

			bool is_conditional = false;
			int max_clique_size;
			bool include_cliques_as_fos_elements;
			bool include_full_fos_element;

			void print();

			uvec order;
			int *next_variable_to_sample = NULL;

			double **S_matrix;
			double *S_vector;                             /* Avoids quadratic memory requirements when a linkage tree is learned based on a random distance measure. */
};

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
extern int FOS_element_ub,                       /* Cut-off value for bounded fixed linkage tree (BFLT). */
          use_univariate_FOS,                   /* Whether a univariate FOS is used. */
          learn_linkage_tree,                   /* Whether the FOS is learned at the start of each generation. */
          static_linkage_tree,                  /* Whether the FOS is fixed throughout optimization. */
          random_linkage_tree;                  /* Whether the fixed linkage tree is learned based on a random distance measure. */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
