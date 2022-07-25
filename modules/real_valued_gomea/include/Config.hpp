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

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

class Config {
	public:
		/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
		Config(); 

		/*-=-=-=-=-=-=-=-=-=-=-=- Options -=-=-=-=-=-=-=-=-=-=-=-=-*/
		int problem_index, number_of_parameters;
		short  print_verbose_overview,                              /* Whether to print a overview of settings (0 = no). */
			   fix_seed,                                            /* Whether a fixed seed is used. */
			   use_vtr,
			   black_box_evaluations,                         /* Whether full (black-box) evaluations must always be performed. */
			   write_generational_statistics,                 /* Whether to compute and write statistics every generation (0 = no). */
			   write_generational_solutions,                  /* Whether to write the population every generation (0 = no). */
			   selection_during_gom,
			   update_elitist_during_gom,
			   verbose;
		int    base_population_size,                                /* The size of the first population in the multi-start scheme. */
			   maximum_number_of_populations,                       /* The maximum number of populations in the multi-start scheme. */
			   number_of_subgenerations_per_population_factor,      /* The subgeneration factor in the multi-start scheme. */
			   maximum_no_improvement_stretch,                      /* The maximum number of subsequent generations without an improvement while the distribution multiplier is <= 1.0. */
			   use_conditional_sampling,
			   FOS_element_size;
		double maximum_number_of_evaluations,                       /* The maximum number of evaluations. */
			   maximum_number_of_seconds,                           /* The maximum number of seconds. */
			   tau,                                                 /* The selection truncation percentile (in [1/population_size,1]). */
			   distribution_multiplier_increase,                    /* The multiplicative distribution multiplier increase. */
			   distribution_multiplier_decrease,                    /* The multiplicative distribution multiplier decrease. */
			   st_dev_ratio_threshold,                              /* The maximum ratio of the distance of the average improvement to the mean compared to the distance of one standard deviation before triggering AVS (SDR mechanism). */
			   fitness_variance_tolerance;                          /* The minimum fitness variance level that is allowed. */
		double vtr,                                           /* The value-to-reach (function value of best solution that is feasible). */
			   lower_user_range,                              /* The initial lower range-bound indicated by the user (same for all dimensions). */
			   upper_user_range;                              /* The initial upper range-bound indicated by the user (same for all dimensions). */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
};

