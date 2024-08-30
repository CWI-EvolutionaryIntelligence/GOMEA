/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/real_valued/rv-gomea.hpp"
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
    this->fitness->maximum_number_of_evaluations = config->maximum_number_of_evaluations;
    this->fitness->maximum_number_of_seconds = config->maximum_number_of_seconds;

    config->checkOptions();
}

rvg_t::rvg_t( int argc, char **argv )
{
	this->config = new Config();
    config->parseCommandLine( argc, argv );

    this->fitness = config->fitness;
    this->fitness->maximum_number_of_evaluations = config->maximum_number_of_evaluations;
    this->fitness->maximum_number_of_seconds = config->maximum_number_of_seconds;

    run();
}

/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void rvg_t::initialize( void )
{
    total_number_of_writes = 0;
    config->number_of_subgenerations_per_population_factor = 8;
    fitness->number_of_evaluations = 0;
    output = output_statistics_t();

    if( config->fix_seed )
    {
        utils::initializeRandomNumberGenerator(config->random_seed);
    }
    else
    {
        utils::initializeRandomNumberGenerator();
    }
    fitness->initializeRun();
	
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
	new_pop->selection_during_gom = config->selection_during_gom;
	new_pop->update_elitist_during_gom = config->update_elitist_during_gom;
    new_pop->linkage_config = config->linkage_config;
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
    new_pop->linkage_config = config->linkage_config;
	new_pop->initialize();
	populations.push_back(new_pop);
}

void rvg_t::initializeProblem()
{
	fitness = config->fitness;
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
std::vector<double> rvg_t::getOverallBestFitness()
{
	double best_obj_val = populations[0]->objective_value_elitist;
	double best_cons_val = populations[0]->constraint_value_elitist;
	for( int i = 1; i < populations.size(); i++ )
	{
		if( fitness->betterFitness( populations[i]->objective_value_elitist, populations[i]->constraint_value_elitist, best_obj_val, best_cons_val ) )
        {
            best_obj_val = populations[i]->objective_value_elitist;
            best_cons_val = populations[i]->constraint_value_elitist;
        }
	}

    std::vector<double> result(2);
    result[0] = best_obj_val;
    result[1] = best_cons_val;
    return( result );
}

void rvg_t::writeGenerationalStatisticsForOnePopulation( int population_index )
{
    /* Average, best and worst */
    /*double population_objective_avg  = populations[population_index]->getFitnessMean();
    double population_constraint_avg = populations[population_index]->getConstraintValueMean();
    double population_objective_var  = populations[population_index]->getFitnessVariance();
    double population_constraint_var = populations[population_index]->getConstraintValueVariance();
    solution_t<double> *best_solution = populations[population_index]->getBestSolution();
    solution_t<double> *worst_solution = populations[population_index]->getWorstSolution();*/

    int key = total_number_of_writes;
    output.addMetricValue("generation",key,populations[population_index]->number_of_generations);
    output.addMetricValue("evaluations",key,fitness->number_of_evaluations);
    output.addMetricValue("time",key,utils::getElapsedTimeSinceStartSeconds());
    output.addMetricValue("eval_time",key,utils::getTimer("eval_time"));
    output.addMetricValue("population_index",key,population_index);
    output.addMetricValue("population_size",key,populations[population_index]->population_size);
    output.addMetricValue("best_genotype",key,getElitist()->variablesToString());
    output.addMetricValue("best_obj_val",key,fitness->elitist_objective_value);
    output.addMetricValue("best_cons_val",key,fitness->elitist_constraint_value);
    //output.addMetricValue("subpop_best_obj_val",key,best_solution->getObjectiveValue());
    //output.addMetricValue("subpop_best_cons_val",key,best_solution->getConstraintValue());
    //output.addMetricValue("subpop_obj_val_avg",key,population_objective_avg);
    //output.addMetricValue("subpop_obj_val_var",key,population_objective_var);
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
void rvg_t::writeGenerationalSolutions( bool final )
{
    char  string[1000];
    FILE *file_all, *file_population, *file_selection;

    file_selection = NULL;
    if( final )
        sprintf( string, "all_populations_generation_final.dat" );
    else
        sprintf( string, "all_populations_generation_%05d.dat", (int) populations.size() );
    file_all = fopen( string, "w" );

    for( size_t i = 0; i < populations.size(); i++ )
    {
        if( final )
            sprintf( string, "population_%05d_generation_final.dat", (int) i );
        else
            sprintf( string, "population_%05d_generation_%05d.dat", (int) i, populations[i]->number_of_generations );
        file_population = fopen( string, "w" );

        //if( populations[i]->number_of_generations > 0 && !final )
        {
			populations[i]->makeSelection();
            sprintf( string, "selection_%05d_generation_%05d.dat", (int) i, populations[i]->number_of_generations );
            file_selection = fopen( string, "w" );
        }

        /* Populations */
        for( int j = 0; j < populations[i]->population_size; j++ )
        {
            for( int k = 0; k < fitness->number_of_variables; k++ )
            {
                sprintf( string, "%13e", populations[i]->individuals[j]->variables[k] );
                fputs( string, file_all );
                fputs( string, file_population );
                if( k < fitness->number_of_variables-1 )
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
            for( int j = 0; j < populations[i]->selection_size; j++ )
            {
                for( int k = 0; k < fitness->number_of_variables; k++ )
                {
                    sprintf( string, "%13e", populations[i]->selection[j][k] );
                    fputs( string, file_selection );
                    if( k < fitness->number_of_variables-1 )
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
void rvg_t::writeGenerationalSolutionsBest( bool final )
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

    for( i = 0; i < fitness->number_of_variables; i++ )
    {
        sprintf( string, "%13e", populations[population_index_best]->individuals[individual_index_best]->variables[i] );
        fputs( string, file );
        if( i < fitness->number_of_variables-1 )
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
bool rvg_t::checkTerminationCondition( void )
{
	if( checkNumberOfEvaluationsTerminationCondition() )
	{
		return( true );
	}

	if( checkVTRTerminationCondition() )
	{
		return( true );
	}

	if( checkTimeLimitTerminationCondition() )
	{
		return( true );
	}

	checkAverageFitnessTerminationConditions();

	if((int) populations.size() < config->maximum_number_of_populations )
    {
		return( false );
    }

	bool allTrue = true;
	for( size_t i = 0; i < populations.size(); i++ )
	{
		if( !populations[i]->population_terminated )
		{
			allTrue = false;
			break;
		}
	}
	if( allTrue )
	{
		restartLargestPopulation();
		allTrue = false;
	}

	return( allTrue );
}

bool rvg_t::checkPopulationTerminationConditions( int population_index )
{
	if( checkNumberOfGenerationsTerminationCondition(population_index) )
	{
        return( true );
    }
	
    if( checkFitnessVarianceTermination(population_index) )
	{
        return( true );
    }
    
	if( checkDistributionMultiplierTerminationCondition(population_index) )
	{
        return( true );
    }

    return( false );
}

bool rvg_t::checkSubgenerationTerminationConditions()
{
	if( checkNumberOfEvaluationsTerminationCondition() )
    {
        return( true );
    }

    if( checkVTRTerminationCondition() )
    {
        return( true );
    }

    if( checkTimeLimitTerminationCondition() )
    {
        return( true );
    }

	return( false );
}

bool rvg_t::checkTimeLimitTerminationCondition( void )
{
    return( config->maximum_number_of_seconds > 0 && utils::getElapsedTimeSinceStartSeconds() > config->maximum_number_of_seconds );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
bool rvg_t::checkNumberOfEvaluationsTerminationCondition( void )
{
    if( fitness->number_of_evaluations >= config->maximum_number_of_evaluations && config->maximum_number_of_evaluations > 0 )
        return( true );

    return( false );
}

/**
 * Returns 1 if the value-to-reach has been reached (in any population).
 */
bool rvg_t::checkVTRTerminationCondition( void )
{
    return( config->use_vtr && fitness->vtr_hit_status );
}

void rvg_t::checkAverageFitnessTerminationConditions( void )
{
    double *average_objective_values = (double*) Malloc( populations.size() * sizeof(double) );
    double *average_constraint_values = (double*) Malloc( populations.size() * sizeof(double) );
    for(int i = ((int)populations.size())-1; i >= 0; i-- )
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
        if( i < ((int)populations.size())-1 && fitness->betterFitness(average_objective_values[i+1], average_constraint_values[i+1], average_objective_values[i], average_constraint_values[i]) )
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
    for(size_t i = 0; i < populations.size(); i++ )
    {
        for(int j = 0; j < populations[i]->population_size; j++ )
        {
            if( fitness->betterFitness( populations[i]->individuals[j]->getObjectiveValue(), populations[i]->individuals[j]->getConstraintValue(),
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
bool rvg_t::checkFitnessVarianceTermination( int population_index )
{
	if( populations[population_index]->getFitnessVariance() < config->fitness_variance_tolerance * populations[population_index]->getFitnessMean() )
		return( true );
	return( false );
}

bool rvg_t::checkNumberOfGenerationsTerminationCondition( int population_index )
{
    if( config->maximum_number_of_generations > 0 && populations[population_index]->number_of_generations >= config->maximum_number_of_generations )
        return( true );
    return( false );
}

/**
 * Checks whether the distribution multiplier in any population
 * has become too small (1e-10).
 */
bool rvg_t::checkDistributionMultiplierTerminationCondition( int population_index )
{
	int i = population_index;
	if( !populations[i]->population_terminated )
	{
		bool converged = true;
		for(int j = 0; j < populations[i]->linkage_model->size(); j++ )
		{
			if( populations[i]->linkage_model->getDistributionMultiplier(j) > 1e-10 )
			{
				converged = false;
				break;
			}
		}

		if( converged )
			return( true );
	}
	return( false );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
rvg_t::~rvg_t()
{
    for (size_t i = 0; i < populations.size(); i++)
        delete (populations[i]);
}

void rvg_t::ezilaitini()
{
    for (size_t i = 0; i < populations.size(); i++)
        delete (populations[i]);
    populations.clear();
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

				//if( config->write_generational_statistics )
				//if( populations.size() == 1 && config->write_generational_statistics )
				if( populations[population_index]->number_of_generations % 10 == 1 && config->write_generational_statistics )
					writeGenerationalStatisticsForOnePopulation( population_index );

				if( populations.size() == 1 && config->write_generational_solutions )
					writeGenerationalSolutions( 0 );

                if( checkSubgenerationTerminationConditions() )
                {
                    for(size_t j = 0; j < populations.size(); j++ )
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
        if( (int) populations.size() < config->maximum_number_of_populations )
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

solution_t<double> *rvg_t::getElitist()
{
    solution_t<double> *best_so_far = populations[0]->getElitist();
    for( int i = 1; i < populations.size(); i++ )
    {
        solution_t<double> *pop_elitist = populations[i]->getElitist();
        if( fitness->betterFitness( pop_elitist, best_so_far ) )
            best_so_far = pop_elitist;
    }
    return best_so_far;
}

/**
 * Runs the IDEA.
 */
void rvg_t::run( void )
{
    //printf("Running RV-GOMEA\n");

    //int out = gomea::utils::initializePythonEmbedding("gomea", PyInit_real_valued);
    //assert(out == 0);

    utils::initStartTime();
	utils::clearTimers();

	initialize();

	if( config->print_verbose_overview )
		config->printVerboseOverview();

	try
    {
        runAllPopulations();
    }
	catch( utils::customException const& ){
        for( auto &p : populations )
            p->updateElitist();
    }
	writeGenerationalStatisticsForOnePopulation( populations.size()-1 );
    ezilaitini();

	/*printf("evals %f ", fitness->number_of_evaluations);
	printf("obj_val %6.2e ", fitness->elitist_objective_value);
	printf("time %lf ", utils::getElapsedTimeSinceStartSeconds());
	printf("generations ");
	if( populations.size() > 1 )
	{
		printf("[ ");
		for(size_t i = 0; i < populations.size(); i++ )
			printf("%d ",populations[i]->number_of_generations);
		printf("]");
	}
	else    
		printf("%d", populations[0]->number_of_generations );
	printf(" random_seed %ld",utils::random_seed);
	printf("\n");*/
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
