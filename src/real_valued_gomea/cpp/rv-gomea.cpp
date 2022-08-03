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

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "rv-gomea.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

std::vector<population_t*> populations;

/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
rvg_t::rvg_t()
{}

rvg_t::rvg_t( Config *config )
{
	this->config = config;
    this->fitness = config->fitness;
    use_univariate_FOS = 0;
    learn_linkage_tree = 0;
    static_linkage_tree = 0;
    random_linkage_tree = 0;
    haveNextNextGaussian = 0;

    printf("TEST : vtr=%10.3e\n",config->fitness->vtr);
    if( use_guidelines )
    {
        config->tau                              = 0.35;
        if( config->maximum_number_of_populations == 1 )
            config->base_population_size         = (int) (36.1 + 7.58*log2((double) fitness->number_of_parameters));
        else
            config->base_population_size         = 10;
        //config->base_population_size           = (int) (10.0*pow((double) number_of_parameters,0.5));
        //config->base_population_size           = (int) (17.0 + 3.0*pow((double) number_of_parameters,1.5));
        //config->base_population_size           = (int) (4.0*pow((double) number_of_parameters,0.5));
        config->distribution_multiplier_decrease = 0.9;
        config->st_dev_ratio_threshold           = 1.0;
        config->maximum_no_improvement_stretch   = 25 + fitness->number_of_parameters;
    }
    FOS_element_ub = fitness->number_of_parameters;
    if( config->FOS_element_size == -1 ) config->FOS_element_size = fitness->number_of_parameters;
    if( config->FOS_element_size == -2 ) learn_linkage_tree = 1;
    if( config->FOS_element_size == -3 ) static_linkage_tree = 1;
    if( config->FOS_element_size == -4 ) {static_linkage_tree = 1; FOS_element_ub = 100;}
    if( config->FOS_element_size == -5 ) {random_linkage_tree = 1; static_linkage_tree = 1; FOS_element_ub = 100;}
    if( config->FOS_element_size == 1 ) use_univariate_FOS = 1;

    checkOptions();
}

rvg_t::rvg_t( int argc, char **argv )
{
	this->config = new Config();

    use_univariate_FOS = 0;
    learn_linkage_tree = 0;
    static_linkage_tree = 0;
    random_linkage_tree = 0;
    config->FOS_element_size = -1;
    haveNextNextGaussian = 0;

    parseCommandLine( argc, argv );

    if( use_guidelines )
    {
        config->tau                              = 0.35;
        if( config->maximum_number_of_populations == 1 )
            config->base_population_size         = (int) (36.1 + 7.58*log2((double) fitness->number_of_parameters));
        else
            config->base_population_size         = 10;
        //config->base_population_size           = (int) (10.0*pow((double) number_of_parameters,0.5));
        //config->base_population_size           = (int) (17.0 + 3.0*pow((double) number_of_parameters,1.5));
        //config->base_population_size           = (int) (4.0*pow((double) number_of_parameters,0.5));
        config->distribution_multiplier_decrease = 0.9;
        config->st_dev_ratio_threshold           = 1.0;
        config->maximum_no_improvement_stretch   = 25 + fitness->number_of_parameters;
    }
    FOS_element_ub = fitness->number_of_parameters;
    if( config->FOS_element_size == -1 ) config->FOS_element_size = fitness->number_of_parameters;
    if( config->FOS_element_size == -2 ) learn_linkage_tree = 1;
    if( config->FOS_element_size == -3 ) static_linkage_tree = 1;
    if( config->FOS_element_size == -4 ) {static_linkage_tree = 1; FOS_element_ub = 100;}
    if( config->FOS_element_size == -5 ) {random_linkage_tree = 1; static_linkage_tree = 1; FOS_element_ub = 100;}
    if( config->FOS_element_size == 1 ) use_univariate_FOS = 1;

    checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void rvg_t::parseCommandLine( int argc, char **argv )
{
    int index;

    index = 1;

    parseOptions( argc, argv, &index );

    parseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
void rvg_t::parseOptions( int argc, char **argv, int *index )
{
	config = new Config();
    double dummy;

    config->write_generational_statistics = 0;
    config->write_generational_solutions  = 0;
    config->print_verbose_overview        = 0;
    config->use_vtr                       = 0;
	use_guidelines                = 0;
	config->black_box_evaluations         = 0;
	config->fix_seed					  = 0;

    for( ; (*index) < argc; (*index)++ )
    {
        if( argv[*index][0] == '-' )
        {
            /* If it is a negative number, the option part is over */
            if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
                break;

            if( argv[*index][1] == '\0' )
                optionError( argv, *index );
            else if( argv[*index][2] != '\0' )
                optionError( argv, *index );
            else
            {
                switch( argv[*index][1] )
                {
                case 'h': printUsage(); break;
				case 'P': fitness_t::printAllInstalledProblems(); break;
                case 's': config->write_generational_statistics = 1; break;
                case 'w': config->write_generational_solutions  = 1; break;
                case 'v': config->print_verbose_overview        = 1; break;
                case 'r': config->use_vtr                       = 1; break;
                case 'g': use_guidelines                = 1; break;
                case 'b': config->black_box_evaluations         = 1; break;
                case 'f': parseFOSElementSize( index, argc, argv ); break;
                case 'S': config->fix_seed                      = 1; break;
                default : optionError( argv, *index );
                }
            }
        }
        else /* Argument is not an option, so option part is over */
            break;
    }
}

void rvg_t::parseFOSElementSize( int *index, int argc, char** argv )
{
    short noError = 1;

    (*index)++;
    noError = noError && sscanf( argv[*index], "%d", &config->FOS_element_size );

    if( !noError )
    {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Checks whether the selected options are feasible.
 */
void rvg_t::checkOptions( void )
{
    if( fitness->number_of_parameters < 1 )
    {
        printf("\n");
        printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", fitness->number_of_parameters);
        printf("\n\n");

        exit( 0 );
    }

    if( ((int) (config->tau*config->base_population_size)) <= 0 || config->tau >= 1 )
    {
        printf("\n");
        printf("Error: tau not in range (read: %e). Require tau in [1/pop,1] (read: [%e,%e]).", config->tau, 1.0/((double) config->base_population_size), 1.0);
        printf("\n\n");

        exit( 0 );
    }

    if( config->base_population_size < 1 )
    {
        printf("\n");
        printf("Error: population size < 1 (read: %d). Require population size >= 1.", config->base_population_size);
        printf("\n\n");

        exit( 0 );
    }

    if( config->maximum_number_of_populations < 1 )
    {
        printf("\n");
        printf("Error: number of populations < 1 (read: %d). Require number of populations >= 1.", config->maximum_number_of_populations );
        printf("\n\n");

        exit( 0 );
    }

    if( fitness_t::installedProblemName( config->problem_index ) == NULL )
    {
        printf("\n");
        printf("Error: unknown index for problem (read index %d).", config->problem_index );
        printf("\n\n");

        exit( 0 );
    }

    /*if( rotation_angle > 0 && ( !learn_linkage_tree && config->FOS_element_size > 1 && config->FOS_element_size != block_size && config->FOS_element_size != fitness->number_of_parameters) )
    {
        printf("\n");
        printf("Error: invalid FOS element size (read %d). Must be %d, %d or %d.", config->FOS_element_size, 1, block_size, fitness->number_of_parameters );
        printf("\n\n");

        exit( 0 );
    }*/
}


/**
 * Informs the user of an illegal option and exits the program.
 */
void rvg_t::optionError( char **argv, int index )
{
    printf("Illegal option: %s\n\n", argv[index]);

    printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void rvg_t::parseParameters( int argc, char **argv, int *index )
{
    if( (argc - *index) < 15 )
    {
        printf("Number of parameters is incorrect, require 15 parameters (you provided %d).\n\n", (argc - *index));

        printUsage();
    }

	config->selection_during_gom = 1;
	config->update_elitist_during_gom = 1;

	int a,b;

    int noError = 1;
    noError = noError && sscanf( argv[*index+0], "%d", &config->problem_index );
    noError = noError && sscanf( argv[*index+1], "%d", &fitness->number_of_parameters );
    noError = noError && sscanf( argv[*index+2], "%lf", &config->lower_user_range );
    noError = noError && sscanf( argv[*index+3], "%lf", &config->upper_user_range );
    noError = noError && sscanf( argv[*index+4], "%lf", &rotation_angle );
    noError = noError && sscanf( argv[*index+5], "%lf", &config->tau );
    noError = noError && sscanf( argv[*index+6], "%d", &config->base_population_size );
    noError = noError && sscanf( argv[*index+7], "%d", &config->maximum_number_of_populations );
    noError = noError && sscanf( argv[*index+8], "%lf", &config->distribution_multiplier_decrease );
    noError = noError && sscanf( argv[*index+9], "%lf", &config->st_dev_ratio_threshold );
    noError = noError && sscanf( argv[*index+10], "%lf", &config->maximum_number_of_evaluations );
    noError = noError && sscanf( argv[*index+11], "%lf", &config->vtr );
    noError = noError && sscanf( argv[*index+12], "%d", &config->maximum_no_improvement_stretch );
    noError = noError && sscanf( argv[*index+13], "%lf", &config->fitness_variance_tolerance );
    noError = noError && sscanf( argv[*index+14], "%lf", &config->maximum_number_of_seconds );
	if( argc-*index > 15 )
	{
		noError = noError && sscanf( argv[*index+15], "%d", &a );
    	noError = noError && sscanf( argv[*index+16], "%d", &b );
		config->selection_during_gom = (short) a;
		config->update_elitist_during_gom = (short) b;
	}

    if( !noError )
    {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Prints the settings as read from the command line.
 */
void rvg_t::printVerboseOverview( void )
{
    int i;

    printf("### Settings ######################################\n");
    printf("#\n");
    printf("# Statistics writing every generation: %s\n", config->write_generational_statistics ? "enabled" : "disabled");
    printf("# Population file writing            : %s\n", config->write_generational_solutions ? "enabled" : "disabled");
    printf("# Use of value-to-reach (vtr)        : %s\n", config->use_vtr ? "enabled" : "disabled");
    printf("#\n");
    printf("###################################################\n");
    printf("#\n");
    printf("# Problem                 = %s\n", fitness_t::installedProblemName( config->problem_index ));
    printf("# Number of parameters    = %d\n", fitness->number_of_parameters);
    printf("# Initialization ranges   = [%e;%e]\n", config->lower_user_range, config->upper_user_range );
    printf("# Boundary ranges         = ");
    for( i = 0; i < fitness->number_of_parameters; i++ )
    {
        printf("x_%d: [%e;%e]", i, fitness->getLowerRangeBound(i), fitness->getUpperRangeBound(i) );
        if( i < fitness->number_of_parameters-1 )
            printf("\n#                           ");
    }
    printf("\n");
    printf("# Rotation angle          = %e\n", rotation_angle);
    printf("# Tau                     = %e\n", config->tau);
    printf("# Population size/normal  = %d\n", config->base_population_size);
    printf("# FOS element size        = %d\n", config->FOS_element_size);
    printf("# Max num of populations  = %d\n", config->maximum_number_of_populations);
    printf("# Dis. mult. decreaser    = %e\n", config->distribution_multiplier_decrease);
    printf("# St. dev. rat. threshold = %e\n", config->st_dev_ratio_threshold);
    printf("# Maximum numb. of eval.  = %lf\n", config->maximum_number_of_evaluations);
    printf("# Value to reach (vtr)    = %e\n", config->vtr);
    printf("# Max. no improv. stretch = %d\n", config->maximum_no_improvement_stretch);
    printf("# Fitness var. tolerance  = %e\n", config->fitness_variance_tolerance);
    printf("# Random seed             = %ld\n", random_seed);
    printf("#\n");
    printf("###################################################\n");
}

/**
 * Prints usage information and exits the program.
 */
void rvg_t::printUsage( void )
{
    printf("Usage: RV-GOMEA [-?] pro dim low upp rot tau pop nop dmd srt eva vtr imp tol\n");
    printf(" -h: Prints out this usage information.\n");
    printf(" -P: Prints out a list of all installed optimization problems.\n");
    printf(" -s: Enables computing and writing of statistics every generation.\n");
    printf(" -w: Enables writing of solutions and their fitnesses every generation.\n");
    printf(" -v: Enables verbose mode. Prints the settings before starting the run.\n");
    printf(" -r: Enables use of vtr in termination condition (value-to-reach).\n");
    printf(" -b: Enables counting every partial evaluation as a full evaluation.\n");
    printf(" -f %%d: Sets linkage model that is used. Positive: Use a FOS with elements of %%d consecutive variables.\n");
    printf("     Use -1 for full linkage model, -2 for dynamic linkage tree learned from the population, -3 for fixed linkage tree learned from distance measure,\n");
    printf("     -4 for bounded fixed linkage tree learned from distance measure, -5 for fixed bounded linkage tree learned from random distance measure.\n");
    printf(" -g: Uses guidelines to override parameter settings for those parameters\n");
    printf("     for which a guideline is known in literature. These parameters are:\n");
    printf("     tau pop dmd srt imp\n");
    ;printf(" -S: A fixed random seed is used.\n");

    printf("\n");
    printf("  pro: Index of optimization problem to be solved (minimization).\n");
    printf("  dim: Number of parameters.\n");
    printf("  low: Overall initialization lower bound.\n");
    printf("  upp: Overall initialization upper bound.\n");
    printf("  rot: The angle by which to rotate the problem.\n");
    printf("  tau: Selection percentile (tau in [1/pop,1], truncation selection).\n");
    printf("  pop: Population size per normal.\n");
    printf("  nop: The number of populations (parallel runs that initially partition the search space).\n");
    printf("  dmd: The distribution multiplier decreaser (in (0,1), increaser is always 1/dmd).\n");
    printf("  srt: The standard-devation ratio threshold for triggering variance-scaling.\n");
    printf("  eva: Maximum number of evaluations allowed.\n");
    printf("  vtr: The value to reach. If the objective value of the best feasible solution reaches\n");
    printf("       this value, termination is enforced (if -r is specified).\n");
    printf("  imp: Maximum number of subsequent generations without an improvement while the\n");
    printf("       the distribution multiplier is <= 1.0.\n");
    printf("  tol: The tolerance level for fitness variance (i.e. minimum fitness variance)\n");
    printf("  sec: The time limit in seconds.\n");
    exit( 0 );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void rvg_t::initialize( void )
{
	if( is_initialized )
		return;
	is_initialized = true;

    total_number_of_writes = 0;
    config->number_of_subgenerations_per_population_factor = 8;

    if( config->fix_seed ) random_seed = 14627;
    initializeRandomNumberGenerator();

	initializeProblem();
}

void rvg_t::restartLargestPopulation()
{
	int pop_size = populations[populations.size()-1]->population_size;
	population_t *new_pop = new population_t( fitness, pop_size, config->lower_user_range, config->upper_user_range );
	new_pop->maximum_no_improvement_stretch = config->maximum_no_improvement_stretch;
	new_pop->st_dev_ratio_threshold = config->st_dev_ratio_threshold;
	new_pop->distribution_multiplier_decrease = config->distribution_multiplier_decrease;
	new_pop->maximum_no_improvement_stretch = config->maximum_no_improvement_stretch;
	new_pop->tau = config->tau;
	new_pop->initialize();
	delete( populations[populations.size()-1] );
	populations[populations.size()-1] = new_pop;
}


void rvg_t::initializeNewPopulation()
{
	int new_pop_size = config->base_population_size;
	if( populations.size() > 0 ) new_pop_size = 2 * populations[populations.size()-1]->population_size;
	population_t *new_pop = new population_t( fitness, new_pop_size, config->lower_user_range, config->upper_user_range );
	new_pop->maximum_no_improvement_stretch = config->maximum_no_improvement_stretch;
	new_pop->st_dev_ratio_threshold = config->st_dev_ratio_threshold;
	new_pop->distribution_multiplier_decrease = config->distribution_multiplier_decrease;
	new_pop->maximum_no_improvement_stretch = config->maximum_no_improvement_stretch;
	new_pop->tau = config->tau;
	new_pop->selection_during_gom = config->selection_during_gom;
	new_pop->update_elitist_during_gom = config->update_elitist_during_gom;
	new_pop->FOS_element_size = config->FOS_element_size;
	new_pop->initializeProblem(config->problem_index);
	new_pop->initialize();
	populations.push_back(new_pop);
}

void rvg_t::initializeProblem()
{
	fitness = fitness_t::getFitnessClass( config->problem_index, fitness->number_of_parameters, config->vtr );
	if( fitness == NULL )
	{
		printf("Unknown problem index.\n");
		exit(0);
	}
	fitness->use_vtr = config->use_vtr;
	fitness->black_box_optimization = config->black_box_evaluations;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void rvg_t::writeGenerationalStatisticsForOnePopulation( int population_index )
{
    char    string[1000];
    FILE   *file;

    /* Average, best and worst */
    double population_objective_avg    = 0.0;
    double population_constraint_avg   = 0.0;
    double population_objective_best   = populations[population_index]->individuals[0]->getObjectiveValue();
    double population_constraint_best  = populations[population_index]->individuals[0]->getConstraintValue();
    double population_objective_worst  = populations[population_index]->individuals[0]->getObjectiveValue();
    double population_constraint_worst = populations[population_index]->individuals[0]->getConstraintValue();
    for(int j = 0; j < populations[population_index]->population_size; j++ )
    {
        population_objective_avg  += populations[population_index]->individuals[j]->getObjectiveValue();
        population_constraint_avg += populations[population_index]->individuals[j]->getConstraintValue();
        if( fitness_t::betterFitness( population_objective_worst, population_constraint_worst, populations[population_index]->individuals[j]->getObjectiveValue(), populations[population_index]->individuals[j]->getConstraintValue() ) )
        {
            population_objective_worst = populations[population_index]->individuals[j]->getObjectiveValue();
            population_constraint_worst = populations[population_index]->individuals[j]->getConstraintValue();
        }
        if( fitness_t::betterFitness( populations[population_index]->individuals[j]->getObjectiveValue(), populations[population_index]->individuals[j]->getConstraintValue(), population_objective_best, population_constraint_best ) )
        {
            population_objective_best = populations[population_index]->individuals[j]->getObjectiveValue();
            population_constraint_best = populations[population_index]->individuals[j]->getConstraintValue();
        }
    }
    population_objective_avg  = population_objective_avg / ((double) populations[population_index]->population_size);
    population_constraint_avg = population_constraint_avg / ((double) populations[population_index]->population_size);

    /* Variance */
    double population_objective_var    = 0.0;
    double population_constraint_var   = 0.0;
    for(int j = 0; j < populations[population_index]->population_size; j++ )
    {
        population_objective_var  += (populations[population_index]->individuals[j]->getObjectiveValue() - population_objective_avg)*(populations[population_index]->individuals[j]->getObjectiveValue() - population_objective_avg);
        population_constraint_var += (populations[population_index]->individuals[j]->getConstraintValue() - population_constraint_avg)*(populations[population_index]->individuals[j]->getConstraintValue() - population_constraint_avg);
    }
    population_objective_var  = population_objective_var / ((double) populations[population_index]->population_size);
    population_constraint_var = population_constraint_var / ((double) populations[population_index]->population_size);

    if( population_objective_var <= 0.0 )
        population_objective_var = 0.0;
    if( population_constraint_var <= 0.0 )
        population_constraint_var = 0.0;

    /* Then write them */
    file = NULL;
    if( total_number_of_writes == 0 )
    {
        file = fopen( "statistics.dat", "w" );

        sprintf( string, "# Generation  Evaluations  Time(s)  Best-obj. Best-cons. [Pop.index  Subgen.  Pop.size  Dis.mult.[0]  Pop.best.obj. Pop.avg.obj.  Pop.var.obj. Pop.worst.obj.  Pop.best.con. Pop.avg.con.  Pop.var.con. Pop.worst.con.]\n" );
        fputs( string, file );
    }
    else
        file = fopen( "statistics.dat", "a" );

    sprintf( string, "%10d %11lf %11.3lf %20.15e %13e  ", populations.size(), fitness->number_of_evaluations, getTimer(), fitness->elitist_objective_value, fitness->elitist_constraint_value );
    fputs( string, file );

    //sprintf( string, "[ %4d %6d %10d %13e %13e %13e %13e %13e %13e %13e %13e %13e ]", population_index, number_of_generations[population_index], population_sizes[population_index], distribution_multipliers[population_index][0], population_objective_best, population_objective_avg, population_objective_var, population_objective_worst, population_constraint_best, population_constraint_avg, population_constraint_var, population_constraint_worst );
    //fputs( string, file );

    sprintf( string, "\n");
    fputs( string, file );

    fclose( file );

    total_number_of_writes++;
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".
 *
 * all_populations_generation_xxxxx.dat : all populations combined
 * population_xxxxx_generation_xxxxx.dat: the individual populations
 * selection_xxxxx_generation_xxxxx.dat : the individual selections
 */
void rvg_t::writeGenerationalSolutions( short final )
{
    int   i, j, k;
    char  string[1000];
    FILE *file_all, *file_population, *file_selection;

    file_selection = NULL;
    if( final )
        sprintf( string, "all_populations_generation_final.dat" );
    else
        sprintf( string, "all_populations_generation_%05d.dat", populations.size() );
    file_all = fopen( string, "w" );

    for( i = 0; i < populations.size(); i++ )
    {
        if( final )
            sprintf( string, "population_%05d_generation_final.dat", i );
        else
            sprintf( string, "population_%05d_generation_%05d.dat", i, populations[i]->number_of_generations );
        file_population = fopen( string, "w" );

        //if( populations[i]->number_of_generations > 0 && !final )
        {
			populations[i]->makeSelection();
            sprintf( string, "selection_%05d_generation_%05d.dat", i, populations[i]->number_of_generations );
            file_selection = fopen( string, "w" );
        }

        /* Populations */
        for( j = 0; j < populations[i]->population_size; j++ )
        {
            for( k = 0; k < fitness->number_of_parameters; k++ )
            {
                sprintf( string, "%13e", populations[i]->individuals[j]->variables[k] );
                fputs( string, file_all );
                fputs( string, file_population );
                if( k < fitness->number_of_parameters-1 )
                {
                    sprintf( string, " " );
                    fputs( string, file_all );
                    fputs( string, file_population );
                }
            }
            sprintf( string, "     " );
            fputs( string, file_all );
            fputs( string, file_population );
            sprintf( string, "%13e %13e", populations[i]->individuals[j]->getObjectiveValue(), populations[i]->individuals[j]->getConstraintValue() );
            fputs( string, file_all );
            fputs( string, file_population );
            sprintf( string, "\n" );
            fputs( string, file_all );
            fputs( string, file_population );
        }

        fclose( file_population );

        /* Selections */
        /*if( populations[i]->number_of_generations > 0 && !final )
        {
            for( j = 0; j < populations[i]->selection_size; j++ )
            {
                for( k = 0; k < fitness->number_of_parameters; k++ )
                {
                    sprintf( string, "%13e", populations[i]->selection[j][k] );
                    fputs( string, file_selection );
                    if( k < fitness->number_of_parameters-1 )
                    {
                        sprintf( string, " " );
                        fputs( string, file_selection );
                    }
                    sprintf( string, "     " );
                    fputs( string, file_selection );
                }
                sprintf( string, "%13e %13e", populations[i]->objective_value_selection[j], populations[i]->constraint_values_selection[j] );
                fputs( string, file_selection );
                sprintf( string, "\n" );
                fputs( string, file_selection );
            }
            fclose( file_selection );
        }*/
    }

    fclose( file_all );

    writeGenerationalSolutionsBest( final );
}


/**
 * Writes the best solution (measured in the single
 * available objective) to a file named
 * best_generation_xxxxx.dat where xxxxx is the
 * generation number. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".The output
 * file contains the solution values with the
 * dimensions separated by a single white space,
 * followed by five white spaces and then the
 * single objective value for that solution
 * and its sum of constraint violations.
 */
void rvg_t::writeGenerationalSolutionsBest( short final )
{
    int   i, population_index_best, individual_index_best;
    char  string[1000];
    FILE *file;
    static int *c = NULL;
    if( c == NULL ){c = (int *) Malloc( sizeof( int ) ); c[0]=0;}

    /* First find the best of all */
    determineBestSolutionInCurrentPopulations( &population_index_best, &individual_index_best );

    /* Then output it */
    if( final )
        sprintf( string, "best_generation_final.dat" );
    else
        sprintf( string, "best_generation_%05d.dat", c[0] );
    file = fopen( string, "w" );

    for( i = 0; i < fitness->number_of_parameters; i++ )
    {
        sprintf( string, "%13e", populations[population_index_best]->individuals[individual_index_best]->variables[i] );
        fputs( string, file );
        if( i < fitness->number_of_parameters-1 )
        {
            sprintf( string, " " );
            fputs( string, file );
        }
    }
    sprintf( string, "     " );
    fputs( string, file );
    sprintf( string, "%13e %13e", populations[population_index_best]->individuals[individual_index_best]->getObjectiveValue(), populations[population_index_best]->individuals[individual_index_best]->getConstraintValue() );
    fputs( string, file );
    sprintf( string, "\n" );
    fputs( string, file );
    c[0]++;

    fclose( file );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced, 0 otherwise.
 */
short rvg_t::checkTerminationCondition( void )
{
	short allTrue;
	int   i;

	if( checkNumberOfEvaluationsTerminationCondition() )
	{
		return( 1 );
	}

	if( checkVTRTerminationCondition() )
	{
		return( 1 );
	}

	if( checkTimeLimitTerminationCondition() )
	{
		return( 1 );
	}

	checkAverageFitnessTerminationConditions();

	if( populations.size() < config->maximum_number_of_populations )
		return( 0 );

	allTrue = 1;
	for( i = 0; i < populations.size(); i++ )
	{
		if( !populations[i]->population_terminated )
		{
			allTrue = 0;
			break;
		}
	}
	if( allTrue )
	{
		restartLargestPopulation();
		allTrue = 0;
	}

	return( allTrue );
}

short rvg_t::checkPopulationTerminationConditions( int population_index )
{
	if( checkFitnessVarianceTermination(population_index) )
		return( 1 );
    
	if( checkDistributionMultiplierTerminationCondition(population_index) )
		return( 1 );

    return( 0 );
}

short rvg_t::checkSubgenerationTerminationConditions()
{
	if( checkNumberOfEvaluationsTerminationCondition() )
        return( 1 );

    if( checkVTRTerminationCondition() )
        return( 1 );

    if( checkTimeLimitTerminationCondition() )
        return( 1 );
    
	return( 0 );
}

short rvg_t::checkTimeLimitTerminationCondition( void )
{
    return( config->maximum_number_of_seconds > 0 && getTimer() > config->maximum_number_of_seconds );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
short rvg_t::checkNumberOfEvaluationsTerminationCondition( void )
{
    if( fitness->number_of_evaluations >= config->maximum_number_of_evaluations && config->maximum_number_of_evaluations > 0 )
        return( 1 );

    return( 0 );
}

/**
 * Returns 1 if the value-to-reach has been reached (in any population).
 */
short rvg_t::checkVTRTerminationCondition( void )
{
    return( config->use_vtr && fitness->vtr_hit_status );
}

void rvg_t::checkAverageFitnessTerminationConditions( void )
{
    double *average_objective_values = (double*) Malloc( populations.size() * sizeof(double) );
    double *average_constraint_values = (double*) Malloc( populations.size() * sizeof(double) );
    for(int i = populations.size()-1; i >= 0; i-- )
    {
        average_objective_values[i] = 0;
        average_constraint_values[i] = 0;
        for(int j = 0; j < populations[i]->population_size; j++ )
        {
            average_objective_values[i] += populations[i]->individuals[j]->getObjectiveValue();
            average_constraint_values[i] += populations[i]->individuals[j]->getConstraintValue();
        }
        average_objective_values[i] /= populations[i]->population_size;
        average_constraint_values[i] /= populations[i]->population_size;
        if( i < populations.size()-1 && fitness_t::betterFitness(average_objective_values[i+1], average_constraint_values[i+1], average_objective_values[i], average_constraint_values[i]) )
        {
            for(int j = i; j >= 0; j-- )
                populations[j]->population_terminated = 1;
            break;
        }
    }
    free( average_objective_values );
    free( average_constraint_values );
}

/**
 * Determines which solution is the best of all solutions
 * in all current populations.
 */
void rvg_t::determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best )
{
    (*population_of_best) = 0;
    (*index_of_best)      = 0;
    for(int i = 0; i < populations.size(); i++ )
    {
        for(int j = 0; j < populations[i]->population_size; j++ )
        {
            if( fitness_t::betterFitness( populations[i]->individuals[j]->getObjectiveValue(), populations[i]->individuals[j]->getConstraintValue(),
                               populations[(*population_of_best)]->individuals[(*index_of_best)]->getObjectiveValue(), populations[(*population_of_best)]->individuals[(*index_of_best)]->getConstraintValue() ) )
            {
                (*population_of_best) = i;
                (*index_of_best)      = j;
            }
        }
    }
}

/**
 * Checks whether the fitness variance in any population
 * has become too small (user-defined tolerance).
 */
short rvg_t::checkFitnessVarianceTermination( int population_index )
{
	if( populations[population_index]->getFitnessVariance() < config->fitness_variance_tolerance * populations[population_index]->getFitnessMean() )
		return( 1 );
	return( 0 );
}


/**
 * Checks whether the distribution multiplier in any population
 * has become too small (1e-10).
 */
short rvg_t::checkDistributionMultiplierTerminationCondition( int population_index )
{
	int i = population_index;
	if( !populations[i]->population_terminated )
	{
		short converged = 1;
		for(int j = 0; j < populations[i]->linkage_model->getLength(); j++ )
		{
			if( populations[i]->linkage_model->getDistributionMultiplier(j) > 1e-10 )
			{
				converged = 0;
				break;
			}
		}

		if( converged )
			return( 1 );
	}
	return( 0 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
rvg_t::~rvg_t()
{
	if( is_initialized )
	{
		for( int i = 0; i < populations.size(); i++ )
			delete( populations[i] );
		delete( fitness );
	}
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void rvg_t::generationalStepAllPopulations()
{
    int population_index_smallest, population_index_biggest;

    population_index_biggest  = populations.size()-1;
    population_index_smallest = 0;
    while( population_index_smallest <= population_index_biggest )
    {
        if( !populations[population_index_smallest]->population_terminated )
            break;

        population_index_smallest++;
    }

    generationalStepAllPopulationsRecursiveFold( population_index_smallest, population_index_biggest );
}

void rvg_t::generationalStepAllPopulationsRecursiveFold( int population_index_smallest, int population_index_biggest )
{
    for(int i = 0; i < config->number_of_subgenerations_per_population_factor-1; i++ )
    {
        for(int population_index = population_index_smallest; population_index <= population_index_biggest; population_index++ )
        {
            if( !populations[population_index]->population_terminated )
            {
				populations[population_index]->runGeneration();

				if( populations.size() == 1 && config->write_generational_statistics )
					writeGenerationalStatisticsForOnePopulation( 0 );

				if( populations.size() == 1 && config->write_generational_solutions )
					writeGenerationalSolutions( 0 );

                if( checkSubgenerationTerminationConditions() )
                {
                    for(int j = 0; j < populations.size(); j++ )
                        populations[j]->population_terminated = 1;
                    return;
                }

				if( checkPopulationTerminationConditions(population_index) )
					populations[population_index]->population_terminated = 1;
            }
        }

        for(int population_index = population_index_smallest; population_index < population_index_biggest; population_index++ )
            generationalStepAllPopulationsRecursiveFold( population_index_smallest, population_index );
    }
}

void rvg_t::runAllPopulations()
{
    while( !checkTerminationCondition() )
    {
        if( populations.size() < config->maximum_number_of_populations )
        {
            initializeNewPopulation();
            
			if( populations.size() == 1 && config->write_generational_statistics )
				writeGenerationalStatisticsForOnePopulation( 0 );

            if( populations.size() == 1 && config->write_generational_solutions )
                writeGenerationalSolutions( 0 );
        }

        generationalStepAllPopulations();

        if( populations.size() > 1 && config->write_generational_statistics )
            writeGenerationalStatisticsForOnePopulation( populations.size()-1 );

        if( populations.size() > 1 && config->write_generational_solutions )
            writeGenerationalSolutions( 0 );
    }
}

/**
 * Runs the IDEA.
 */
void rvg_t::run( void )
{
    startTimer();
	initialize();

	if( config->print_verbose_overview )
		printVerboseOverview();

	runAllPopulations();

	printf("evals %f ", fitness->number_of_evaluations);

	printf("obj_val %6.2e ", fitness->elitist_objective_value);

	printf("time %lf ", getTimer());
	printf("generations ");
	if( populations.size() > 1 )
	{
		printf("[ ");
		for(int i = 0; i < populations.size(); i++ )
			printf("%d ",populations[i]->number_of_generations);
		printf("]");
	}
	else    
		printf("%d", populations[0]->number_of_generations );
	printf(" random_seed %ld",random_seed);
	printf("\n");
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main( int argc, char **argv )
{
	rvg_t rvgomea = rvg_t(argc, argv);

	rvgomea.run();
	
	return( 0 );
}

}}
