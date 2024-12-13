#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/real_valued/Config.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/real_valued/linkage_model.hpp"
#include "gomea/src/real_valued/solutionRV.hpp"
#include "gomea/src/common/output_statistics.hpp"

#include <deque>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class population_t {
	public:
		/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
	 	population_t( fitness_t<double> *fitness, int population_size, double lower_init, double upper_init, int population_index = -1 );
	 	~population_t();
		
		void runGeneration();
		void makeSelection();
		solution_t<double> *getElitist();
		void updateElitist();
		void computeRanks();
		void makeSelectionUsingDiversityOnRank0();
		void estimateDistribution( int FOS_index );
		void estimateDistributions();
		double estimateMean( int var );
		void updateAMSMeans();
		void copyBestSolutionsToPopulation();
		void getBestInPopulation( int *individual_index );
		void evaluateCompletePopulation();
		void generateAndEvaluateNewSolutions();
		bool checkForImprovement( solution_t<double> *solution, partial_solution_t<double> *part );
		partial_solution_t<double> *generatePartialSolution( int FOS_index );
		void applyPartialAMS( partial_solution_t<double> *solution, double cmul );
		void generateAndEvaluateNewSolutionsToFillPopulation();
		partial_solution_t<double> *generateNewSolutionFromFOSElement( int FOS_index, int individual_index, bool apply_AMS );
		bool applyAMS( int individual_index );
		void applyForcedImprovements( int individual_index, int donor_index );
		double getFitnessMean();
		double getFitnessVariance();
		double getConstraintValueMean();
		double getConstraintValueVariance();
		solution_t<double> *getBestSolution();
		solution_t<double> *getWorstSolution();
		
		void initialize();
		void initializeDefaultParameters();
		void initializeNewPopulationMemory();
		void initializeFOS();
		void initializeFOSFromIndex( int FOSIndex );
		void initializeFOS( linkage_config_t *linkage_config );
		void initializeParameterRangeBounds( double lower_user_range, double upper_user_range );
		void initializeProblem( int problem_index );
		void initializePopulationAndFitnessValues();
		void writeGenerationalStatistics( bool is_final = false );
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
		
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
		double eta_cov,
			   tau,
			   st_dev_ratio_threshold,
			   distribution_multiplier_increase,
			   distribution_multiplier_decrease,
			  *lower_init_ranges,
			  *upper_init_ranges;
		int    maximum_no_improvement_stretch,
			   num_elitists_to_copy = 1,
			   FOS_element_size,
			   problem_index,
			   population_index;
		double delta_AMS;
		bool update_elitist_during_gom, selection_during_gom;
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
		partial_solution_t<double>*** sampled_solutions; 
		int    number_of_generations,
			   population_size,                                /* The size of the first population in the multi-start scheme. */
			   selection_size,                                     /* The size of the selection for each population. */
			   *individual_NIS;                                      /* The number of generations a solution has not improved. */
		solution_t<double> **individuals,
			   		**selection;                                          /* Selected solutions, one for each population. */
		fitness_t<double> *fitness;
		double *ranks,                                               /* Ranks of population members. */
			   objective_value_elitist,                         /* Objective values of selected solutions. */
			   constraint_value_elitist,                        /* Sum of all constraint violations of selected solutions. */
			   *mean_shift_vector,                                   /* The mean vectors of the previous generation, one for each population. */
			   *prev_mean_vector;                                   /* The mean vectors of the previous generation, one for each population. */
		bool  population_terminated;
		linkage_config_t *linkage_config;
		linkage_model_rv_pt linkage_model;
		output_statistics_t *output;
		Config *config;
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
};

}}
