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
#include "gomea/src/real_valued/population.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/real_valued/linkage_model.hpp"
#include "gomea/src/real_valued/Config.hpp"
#include "gomea/src/utils/embed.hpp"
#include "gomea/real_valued.h"
#include "gomea/fitness.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

class rvg_t {
	public:
		/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
		rvg_t(); 
		rvg_t(Config *config); 
		rvg_t( int argc, char **argv );
		~rvg_t();

		void run( void );
		void runAllPopulations();
		void generationalStepAllPopulations();
		void generationalStepAllPopulationsRecursiveFold( int population_index_smallest, int population_index_biggest );
		void parseCommandLine( int argc, char **argv );
		void parseOptions( int argc, char **argv, int *index );
		void parseFOSIndex( int *index, int argc, char** argv );
		void optionError( char **argv, int index );
		void parseParameters( int argc, char **argv, int *index );
		void printUsage( void );
		void checkOptions( void );
		void printVerboseOverview( void );
		void initialize( void );
		void initializeParameterRangeBounds( void );
		void initializeMemory( void );
		void initializeNewPopulation( void );
		void initializeProblem(); 
		void restartLargestPopulation();
		fitness_t *getFitnessClass( int problem_index );
		void writeGenerationalStatisticsForOnePopulation( int population_index );
		void writeGenerationalSolutions( bool final );
		void writeGenerationalSolutionsBest( bool final );
		bool checkTerminationCondition( void );
		bool checkSubgenerationTerminationConditions();
		bool checkPopulationTerminationConditions( int population_index );
		bool checkTimeLimitTerminationCondition( void );
		bool checkNumberOfEvaluationsTerminationCondition( void );
		bool checkVTRTerminationCondition( void );
		void checkAverageFitnessTerminationConditions( void );
		void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best );
		bool checkFitnessVarianceTermination( int population_index );
		bool checkNumberOfGenerationsTerminationCondition( int population_index );
		bool checkDistributionMultiplierTerminationCondition( int population_index );
		std::vector<double> getOverallBestFitness();
		void ezilaitini();
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

		/*-=-=-=-=-=-=-=-=-=-=-=- Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
		std::vector<population_t*> populations;
		fitness_t *fitness;
		int total_number_of_writes = 0;                              /* Total number of times a statistics file has been written. */
		output_statistics_t output;
		/*-=-=-=-=-=-=-=-=-=-=-=- Options -=-=-=-=-=-=-=-=-=-=-=-=-*/
		Config *config;
		double eta_ams = 1.0,
			   eta_cov = 1.0;
		bool	use_guidelines = false;	 /* Whether to override parameters with guidelines (for those that exist). */
		double rotation_angle = 0.0; /* The angle of rotation to be applied to the problem. */
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
};

}}
