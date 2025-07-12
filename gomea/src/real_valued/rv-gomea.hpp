#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/real_valued/population.hpp"
#include "gomea/src/real_valued/tools.hpp"
#include "gomea/src/real_valued/linkage_model.hpp"
#include "gomea/src/real_valued/Config.hpp"
#include "gomea/src/fitness/fitness.hpp"
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
		void writeGenerationalStatisticsForOnePopulation( int population_index, bool is_final = false );
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
		solution_t<double> *getElitist();
		void ezilaitini();
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

		/*-=-=-=-=-=-=-=-=-=-=-=- Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
		std::vector<population_t*> populations;
		fitness_t<double> *fitness;
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
