#pragma once
#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-rv.hpp"

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace realvalued{

using gomea::fitness::fitness_t;
using gomea::fitness::fitness_pt;

class Config {
	public:
		/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
		Config(){}; 

		void parseCommandLine( int argc, char **argv );
		void parseOptions( int argc, char **argv, int *index );
		void parseFOSIndex( int *index, int argc, char** argv );
		void checkOptions( void );
		void printUsage( void );
		void printVerboseOverview( void );
		void printInstalledProblems( void );
		void optionError( char **argv, int index );
		void parseParameters( int argc, char **argv, int *index );
		void initializeFOSFromIndex( int FOSIndex );
		/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

		/*-=-=-=-=-=-=-=-=-=-=-=- Options -=-=-=-=-=-=-=-=-=-=-=-=-*/
		int problem_index;
		bool   use_guidelines,
			   use_vtr,
			   print_verbose_overview,                              /* Whether to print a overview of settings (0 = no). */
			   fix_seed,                                            /* Whether a fixed seed is used. */
			   black_box_evaluations,                         /* Whether full (black-box) evaluations must always be performed. */
			   write_generational_statistics,                 /* Whether to compute and write statistics every generation (0 = no). */
			   write_generational_solutions,                  /* Whether to write the population every generation (0 = no). */
			   selection_during_gom,
			   update_elitist_during_gom,
			   verbose;
		int    number_of_variables,									/* The number of variables in the optimization problem. */
			   FOSIndex = -1,											/* The index of the FOS to be used. */					
			   base_population_size,                                /* The size of the first population in the multi-start scheme. */
			   maximum_number_of_populations = 25,                  /* The maximum number of populations in the multi-start scheme. */
			   maximum_number_of_generations = -1,                  /* The maximum number of generations before any population in the multi-start scheme will be terminated. */
			   number_of_subgenerations_per_population_factor = 8,  /* The subgeneration factor in the multi-start scheme. */
			   maximum_no_improvement_stretch = 100;                /* The maximum number of subsequent generations without an improvement while the distribution multiplier is <= 1.0. */
		double maximum_number_of_evaluations = -1,                  /* The maximum number of evaluations. */
			   maximum_number_of_seconds = -1,                      /* The maximum number of seconds. */
			   tau,                                                 /* The selection truncation percentile (in [1/population_size,1]). */
			   distribution_multiplier_increase,                    /* The multiplicative distribution multiplier increase. */
			   distribution_multiplier_decrease,                    /* The multiplicative distribution multiplier decrease. */
			   st_dev_ratio_threshold,                              /* The maximum ratio of the distance of the average improvement to the mean compared to the distance of one standard deviation before triggering AVS (SDR mechanism). */
			   fitness_variance_tolerance;                          /* The minimum fitness variance level that is allowed. */
		double vtr,                                           /* The value-to-reach (function value of best solution that is feasible). */
			   lower_user_range,                              /* The initial lower range-bound indicated by the user (same for all dimensions). */
			   upper_user_range;                              /* The initial upper range-bound indicated by the user (same for all dimensions). */
		fitness_pt<double> fitness;
        linkage_config_t *linkage_config;
		long long random_seed;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
};

}}
