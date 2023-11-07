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
			static std::shared_ptr<linkage_model_rv_t> linkage_tree(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static_ ); 
			static std::shared_ptr<linkage_model_rv_t> marginal_product_model(size_t numberOfvariables_, size_t block_size);
			static std::shared_ptr<linkage_model_rv_t> conditional( size_t number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			static std::shared_ptr<linkage_model_rv_t> custom_fos(size_t numberOfvariables_, const vec_t<vec_t<int>> &FOS);
			static std::shared_ptr<linkage_model_rv_t> from_file(std::string filename);

			void initializeDistributions();

			double getAcceptanceRate(); 
			double getDistributionMultiplier( int element_index );
			
			void addGroup( vec_t<int> group );
			void addGroup( distribution_t *dist );
			void addConditionedGroup( vec_t<int> variables, std::set<int> conditioned_variables );
			void randomizeOrder( const graph_t &variable_interaction_graph ); 

			void learnLinkageTreeFOS(matE covariance_matrix);
			vec_t<vec_t<double>> computeMIMatrix( matE covariance_matrix, int n );
			void inheritDistributionMultipliers( linkage_model_rv_t *other, double *multipliers );
			int *matchFOSElements( linkage_model_rv_t *other );
			int *hungarianAlgorithm( int** similarity_matrix, int dim );
			void hungarianAlgorithmAddToTree(int x, int prevx, bool *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim );
			int determineNearestNeighbour(int index, double **S_matrix, int *mpm_num_ind, int mpm_length );
			void ezilaitini();

			void initializeNormalDistribution(int FOS_index);
			void initializeConditionalDistribution( int FOS_index );

			partial_solution_t<double> *generatePartialSolution( int FOS_index, solution_t<double> *solution_conditioned_on, fitness::fitness_generic_t *fitness_function = NULL );
			void estimateDistributions( solution_t<double> **selection, int selection_size );
			void estimateDistribution( int FOS_index, solution_t<double> **selection, int selection_size );
			void adaptDistributionMultiplier( int FOS_index, partial_solution_t<double> **solutions, int num_solutions );

			vec_t<distribution_t*> distributions;
			int no_improvement_stretch = 0;
			int maximum_no_improvement_stretch = 100;
			
			double p_accept = 0.05;
			//vec_t<uvec> variables_conditioned_on; 

			void print();

			int *next_variable_to_sample = NULL;

			double **S_matrix;
			double *S_vector;                             /* Avoids quadratic memory requirements when a linkage tree is learned based on a random distance measure. */
		
		private:
			linkage_model_rv_t(size_t numberOfVariables_) : linkage_model_t(numberOfVariables_) {};
 	  	   	linkage_model_rv_t(size_t numberOfVariables_, size_t block_size ) : linkage_model_t(numberOfVariables_,block_size){};
		    linkage_model_rv_t(size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS ) : linkage_model_t(numberOfVariables_, FOS){};
			linkage_model_rv_t(size_t numberOfVariables_, int similarityMeasure, bool filtered, int maximumSetSize, bool is_static) : linkage_model_t(numberOfVariables_,similarityMeasure,filtered,maximumSetSize,is_static){};
			linkage_model_rv_t(size_t number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element );
			linkage_model_rv_t(std::string filename) : linkage_model_t(filename){};
	};

typedef std::shared_ptr<linkage_model_rv_t> linkage_model_rv_pt;

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

}}
