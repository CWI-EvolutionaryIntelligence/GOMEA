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
#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/real_valued/distribution.hpp"
#include "gomea/src/real_valued/partial_solutionRV.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class linkage_model_rv_t : public linkage_model_t {
		public:
			linkage_model_rv_t( int problem_index, double **covariance_matrix, int n );
			linkage_model_rv_t( int number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			linkage_model_rv_t( const linkage_model_rv_t &f );
			~linkage_model_rv_t();

			static std::shared_ptr<linkage_model_rv_t> createFOSInstance( const linkage_config_t &config, size_t numberOfVariables, const graph_t &VIG = graph_t() );

			static std::shared_ptr<linkage_model_rv_t> univariate(size_t numberOfvariables_);
			static std::shared_ptr<linkage_model_rv_t> full(size_t numberOfvariables_);
			static std::shared_ptr<linkage_model_rv_t> linkage_tree(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_ );
			static std::shared_ptr<linkage_model_rv_t> marginal_product_model(size_t numberOfvariables_, size_t block_size);
			static std::shared_ptr<linkage_model_rv_t> conditional( size_t number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			static std::shared_ptr<linkage_model_rv_t> custom_fos(size_t numberOfvariables_, const vec_t<vec_t<int>> &FOS);
			static std::shared_ptr<linkage_model_rv_t> from_file(std::string filename);

			void initializeDistributions();

			double getAcceptanceRate(); 
			int getDistributionMultiplier( int element_index );
			
			void addGroup( vec_t<int> group );
			void addGroup( distribution_t *dist );
			void addConditionedGroup( vec_t<int> variables );
			void addConditionedGroup( vec_t<int> variables, std::set<int> conditioned_variables );
			void randomizeOrder( const graph_t &variable_interaction_graph ); 
			vec_t<int> getVIGOrderBreadthFirst( const graph_t &variable_interaction_graph );

			void learnLinkageTreeFOS(mat covariance_matrix);
			vec_t<vec_t<double>> computeMIMatrix( mat covariance_matrix, int n );
			void inheritDistributionMultipliers( linkage_model_rv_t *other, double *multipliers );
			int *matchFOSElements( linkage_model_rv_t *other );
			int *hungarianAlgorithm( int** similarity_matrix, int dim );
			void hungarianAlgorithmAddToTree(int x, int prevx, short *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim );
			int determineNearestNeighbour(int index, double **S_matrix, int *mpm_num_ind, int mpm_length );
			void ezilaitini();

			void initializeNormalDistribution(int FOS_index);
			void initializeConditionalDistribution( int FOS_index );

			partial_solution_t<double> *generatePartialSolution( int FOS_index, solution_t<double> *solution_conditioned_on );
			void estimateDistributions( solution_t<double> **selection, int selection_size );
			void estimateDistribution( int FOS_index, solution_t<double> **selection, int selection_size );
			void adaptDistributionMultiplier( int FOS_index, partial_solution_t<double> **solutions, int num_solutions );

			int number_of_variables = -1;
			vec_t<distribution_t*> distributions;
			int no_improvement_stretch = 0;
			int maximum_no_improvement_stretch = 100;
			
			double p_accept = 0.05;
			vec_t<uvec> variables_conditioned_on; 

			bool is_conditional = false;
			int max_clique_size;
			bool include_cliques_as_fos_elements;
			bool include_full_fos_element;

			void print();

			int *next_variable_to_sample = NULL;

			double **S_matrix;
			double *S_vector;                             /* Avoids quadratic memory requirements when a linkage tree is learned based on a random distance measure. */
		
		private:
			linkage_model_rv_t(size_t numberOfVariables_) : linkage_model_t(numberOfVariables_) {initializeDistributions();};
 	  	   	linkage_model_rv_t(size_t numberOfVariables_, size_t block_size ) : linkage_model_t(numberOfVariables_,block_size){initializeDistributions();};
		    linkage_model_rv_t(size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS ) : linkage_model_t(numberOfVariables_, FOS){initializeDistributions();};
			linkage_model_rv_t(size_t numberOfVariables_, int similarityMeasure, bool filtered, int maximumSetSize) : linkage_model_t(numberOfVariables_,similarityMeasure,filtered,maximumSetSize){initializeDistributions();};
			linkage_model_rv_t(size_t number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			linkage_model_rv_t(std::string filename) : linkage_model_t(filename){initializeDistributions();};
	};

typedef std::shared_ptr<linkage_model_rv_t> linkage_model_rv_pt;

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

}}
