#include "gomea/src/real_valued/Config.hpp"

namespace gomea{
namespace realvalued{

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void Config::parseCommandLine( int argc, char **argv )
{
    int index;

    index = 1;

    parseOptions( argc, argv, &index );

    parseParameters( argc, argv, &index );

    checkOptions();
}

/**
 * Parses only the options from the command line.
 */
void Config::parseOptions( int argc, char **argv, int *index )
{
    double dummy;

    generational_statistics       = false;
    generational_solution         = false;
    print_verbose_overview        = false;
    use_vtr                       = false;
	use_guidelines                = false;
	black_box_evaluations         = false;
	fix_seed					  = false;

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
                case 'P': printInstalledProblems(); break;
                case 's': generational_statistics       = true; break;
                case 'w': generational_solution         = true; break;
                case 'v': print_verbose_overview        = true; break;
                case 'r': use_vtr                       = true; break;
                case 'g': use_guidelines                = true; break;
                case 'b': black_box_evaluations         = true; break;
                case 'f': parseFOSIndex( index, argc, argv ); break;
                case 'S': fix_seed                      = true; break;
                default : optionError( argv, *index );
                }
            }
        }
        else /* Argument is not an option, so option part is over */
            break;
    }
}

void Config::parseFOSIndex( int *index, int argc, char** argv )
{
    bool noError = true;

    (*index)++;
    noError = noError && sscanf( argv[*index], "%d", &FOSIndex );

    if( !noError )
    {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Checks whether the selected options are feasible.
 */
void Config::checkOptions( void )
{
    if( use_guidelines )
    {
        tau                              = 0.35;
        if( maximum_number_of_populations == 1 )
            base_population_size         = (int) (36.1 + 7.58*log2((double) fitness->number_of_variables));
        else
            base_population_size         = 10;
        //base_population_size           = (int) (10.0*pow((double) number_of_variables,0.5));
        //base_population_size           = (int) (17.0 + 3.0*pow((double) number_of_variables,1.5));
        //base_population_size           = (int) (4.0*pow((double) number_of_variables,0.5));
        distribution_multiplier_decrease = 0.9;
        st_dev_ratio_threshold           = 1.0;
        maximum_no_improvement_stretch   = 25 + fitness->number_of_variables;
    }

    if( fitness->number_of_variables < 1 )
    {
        printf("\n");
        printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", fitness->number_of_variables);
        printf("\n\n");

        exit( 0 );
    }

    if( ((int) (tau*base_population_size)) <= 0 || tau >= 1 )
    {
        printf("\n");
        printf("Error: tau not in range (read: %e). Require tau in [1/pop,1] (read: [%e,%e]).", tau, 1.0/((double) base_population_size), 1.0);
        printf("\n\n");

        exit( 0 );
    }

    if( base_population_size < 1 )
    {
        printf("\n");
        printf("Error: population size < 1 (read: %d). Require population size >= 1.", base_population_size);
        printf("\n\n");

        exit( 0 );
    }

    if( maximum_number_of_populations < 1 )
    {
        printf("\n");
        printf("Error: number of populations < 1 (read: %d). Require number of populations >= 1.", maximum_number_of_populations );
        printf("\n\n");

        exit( 0 );
    }
}


/**
 * Informs the user of an illegal option and exits the program.
 */
void Config::optionError( char **argv, int index )
{
    printf("Illegal option: %s\n\n", argv[index]);

    printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void Config::parseParameters( int argc, char **argv, int *index )
{
    if( (argc - *index) < 14 )
    {
        printf("Number of parameters is incorrect, require 14 parameters (you provided %d).\n\n", (argc - *index));

        printUsage();
    }

	selection_during_gom = true;
	update_elitist_during_gom = true;

    bool noError = true;
    noError = noError && sscanf( argv[*index+0], "%d", &problem_index );
    noError = noError && sscanf( argv[*index+1], "%d", &number_of_variables );
    noError = noError && sscanf( argv[*index+2], "%lf", &lower_user_range );
    noError = noError && sscanf( argv[*index+3], "%lf", &upper_user_range );
    noError = noError && sscanf( argv[*index+4], "%lf", &tau );
    noError = noError && sscanf( argv[*index+5], "%d", &base_population_size );
    noError = noError && sscanf( argv[*index+6], "%d", &maximum_number_of_populations );
    noError = noError && sscanf( argv[*index+7], "%lf", &distribution_multiplier_decrease );
    noError = noError && sscanf( argv[*index+8], "%lf", &st_dev_ratio_threshold );
    noError = noError && sscanf( argv[*index+9], "%lf", &maximum_number_of_evaluations );
    noError = noError && sscanf( argv[*index+10], "%lf", &vtr );
    noError = noError && sscanf( argv[*index+11], "%d", &maximum_no_improvement_stretch );
    noError = noError && sscanf( argv[*index+12], "%lf", &fitness_variance_tolerance );
    noError = noError && sscanf( argv[*index+13], "%lf", &maximum_number_of_seconds );

    // Initialize fitness function
    if(black_box_evaluations){
        switch(problem_index){
            case 0:
                fitness = new fitness::sphereFunctionBBO_t(number_of_variables,vtr);
                break;
            case 1:
                fitness = new fitness::rosenbrockFunctionBBO_t(number_of_variables,vtr);
                break;
            case 2:
                fitness = new fitness::SOREBChainStrongBBO_t(number_of_variables,vtr);
                break;
            case 10:
                fitness = new fitness::circlesInASquareBBO_t(number_of_variables,vtr);
                break;
            default:
                throw std::invalid_argument("Invalid problem index.");
        }
    }
    else{
        switch(problem_index){
            case 0:
                fitness = new fitness::sphereFunction_t(number_of_variables,vtr);
                break;
            case 1:
                fitness = new fitness::rosenbrockFunction_t(number_of_variables,vtr);
                break;
            case 2:
                fitness = new fitness::SOREBChainStrong_t(number_of_variables,vtr);
                break;
            default:
                throw std::invalid_argument("Invalid problem index.");
        }
    }

    // Initialize linkage model
    initializeFOSFromIndex(FOSIndex);

    if( !noError )
    {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

void Config::initializeFOSFromIndex( int FOSIndex )
{
	int max_clique_size;
	bool include_cliques_as_fos_elements, include_full_fos_element, filtered_lt;
    int lt_similarity_measure = 0, lt_max_set_size = -1;
	if( FOSIndex == 1 ) // UNIVARIATE
    {
        linkage_config = new linkage_config_t();
    }
    else if( FOSIndex > 1 ) // MPM
	{
		linkage_config = new linkage_config_t(true, FOSIndex);
	}
	else if( FOSIndex == -1 ) // FULL
	{
		linkage_config = new linkage_config_t(true, number_of_variables);
	}
	else if( FOSIndex == -2 ) // Dynamic LT
	{
		linkage_config = new linkage_config_t(lt_similarity_measure, filtered_lt, lt_max_set_size, false);
	}
    else if( FOSIndex == -3 ) // Static LT
	{
		linkage_config = new linkage_config_t(lt_similarity_measure, filtered_lt, lt_max_set_size, true);
	}
	else if( FOSIndex <= -10 ) // CONDITIONAL
	{
		int id = -1 * FOSIndex;
		id /= 10;
		include_full_fos_element = (id%10) == 1;
		id /= 10;
		include_cliques_as_fos_elements = (id%10) == 1;
		id /= 10;
		max_clique_size = id;
        linkage_config = new linkage_config_t(max_clique_size, include_cliques_as_fos_elements, include_full_fos_element);
	}
	else
	{
		throw std::invalid_argument("Invalid FOS index.");
	}
}

/**
 * Prints the settings as read from the command line.
 */
void Config::printVerboseOverview( void )
{
    int i;

    printf("### Settings ######################################\n");
    printf("#\n");
    printf("# Statistics writing every generation: %s\n", generational_statistics ? "enabled" : "disabled");
    printf("# Population file writing            : %s\n", generational_solution ? "enabled" : "disabled");
    printf("# Use of value-to-reach (vtr)        : %s\n", use_vtr ? "enabled" : "disabled");
    printf("#\n");
    printf("###################################################\n");
    printf("#\n");
    printf("# Problem                 = %s\n", fitness->name.c_str());
    printf("# Number of parameters    = %d\n", fitness->number_of_variables);
    printf("# Initialization ranges   = [%e;%e]\n", lower_user_range, upper_user_range );
    printf("# Boundary ranges         = ");
    for( i = 0; i < fitness->number_of_variables; i++ )
    {
        printf("x_%d: [%e;%e]", i, fitness->getLowerRangeBound(i), fitness->getUpperRangeBound(i) );
        if( i < fitness->number_of_variables-1 )
            printf("\n#                           ");
    }
    printf("\n");
    printf("# Tau                     = %e\n", tau);
    printf("# Population size/normal  = %d\n", base_population_size);
    printf("# Max num of populations  = %d\n", maximum_number_of_populations);
    printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
    printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
    printf("# Maximum numb. of eval.  = %lf\n", maximum_number_of_evaluations);
    printf("# Value to reach (vtr)    = %e\n", vtr);
    printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
    printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
    printf("# Random seed             = %ld\n", utils::random_seed);
    printf("#\n");
    printf("###################################################\n");
}

/**
 * Prints usage information and exits the program.
 */
void Config::printUsage( void )
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

void Config::printInstalledProblems( void )
{
    printf("Installed optimization problems:\n");
    printf("  0: Sphere function [GBO,BBO]\n");
    printf("  1: Rosenbrock's function [GBO,BBO]\n");
    printf("  2: SOREB chain function [GBO,BBO]\n");
    printf(" 10: Circles in a square function [BBO]\n");
    exit( 0 );
}

}}
