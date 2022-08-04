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
#include "real_valued_gomea/tools.hpp"
#include "fitness/fitness_basic.hpp"
#include "real_valued_gomea/fos.hpp"
#include "real_valued_gomea/solutionRV.hpp"

#include <deque>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class population_t {
	public:
		/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
	 	population_t( gomea::fitness::fitness_t *fitness, int population_size, double lower_init, double upper_init );
	 	~population_t();
		
		void runGeneration();
		void makeSelection();
		void updateElitist();
		void computeRanks();
		void makeSelectionUsingDiversityOnRank0();
		void estimateDistribution( int FOS_index );
		void estimateDistribution();
		double estimateMean( int var );
		void updateAMSMeans();
		void copyBestSolutionsToPopulation();
		void getBestInPopulation( int *individual_index );
		void evaluateCompletePopulation();
		void generateAndEvaluateNewSolutions();
		void insertImprovement( solution_t<double> *solution, partial_solution_t<double> *part );
		short checkForImprovement( solution_t<double> *solution, partial_solution_t<double> *part );
		partial_solution_t<double> *generatePartialSolution( int FOS_index );
		void applyPartialAMS( partial_solution_t<double> *solution, double cmul );
		void generateAndEvaluateNewSolutionsToFillPopulation();
		partial_solution_t<double> *generateNewSolutionFromFOSElement( int FOS_index, int individual_index, short apply_AMS );
		short applyAMS( int individual_index );
		void applyForcedImprovements( int individual_index, int donor_index );
		double getFitnessMean();
		double getFitnessVariance();
		
		void initialize();
		void initializeDefaultParameters();
		void initializeNewPopulationMemory();
		void initializeFOS();
		void initializeParameterRangeBounds( double lower_user_range, double upper_user_range );
		void initializeProblem( int problem_index );
		void initializePopulationAndFitnessValues();
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
			   FOS_element_size;
		int    problem_index;
		double delta_AMS;
		short update_elitist_during_gom, selection_during_gom;
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
		partial_solution_t<double>*** sampled_solutions; 
		int    number_of_generations,
			   population_size,                                /* The size of the first population in the multi-start scheme. */
			   selection_size,                                     /* The size of the selection for each population. */
			   *individual_NIS;                                      /* The number of generations a solution has not improved. */
		solution_t<double> **individuals,
			   		**selection;                                          /* Selected solutions, one for each population. */
		gomea::fitness::fitness_t *fitness;
		double *ranks,                                               /* Ranks of population members. */
			   objective_value_elitist,                         /* Objective values of selected solutions. */
			   constraint_value_elitist,                        /* Sum of all constraint violations of selected solutions. */
			   *mean_shift_vector,                                   /* The mean vectors of the previous generation, one for each population. */
			   *prev_mean_vector;                                   /* The mean vectors of the previous generation, one for each population. */
		short  population_terminated;
		fos_t *linkage_model;
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
};

}}
