#include <iostream>
#include <fstream>
using namespace std;

#include "gomea/src/discrete/gomeaIMS.hpp"
#include "gomea/src/discrete/utils.hpp"

namespace gomea{
namespace discrete{

gomeaIMS::gomeaIMS()
{
	return;
}

gomeaIMS::gomeaIMS(Config *config_): config(config_)
{
    maximumNumberOfGOMEAs   = config->maximumNumberOfGOMEAs;
    IMSsubgenerationFactor  = config->IMSsubgenerationFactor;
    basePopulationSize      = config->basePopulationSize;
	problemInstance 		= config->fitness;
    problemInstance->maximum_number_of_evaluations = config->maximumNumberOfEvaluations;
    problemInstance->maximum_number_of_seconds = config->maximumNumberOfSeconds;
	if( config->fix_seed )
		utils::initializeRandomNumberGenerator(config->randomSeed);
	else
		utils::initializeRandomNumberGenerator();

}

gomeaIMS::~gomeaIMS()
{}

void gomeaIMS::ezilaitini()
{
	for (int i = 0; i < GOMEAs.size(); ++i)
		delete GOMEAs[i];
	GOMEAs.clear();

	// delete problemInstance;
	if( isInitialized )
		delete sharedInformationInstance;
	isInitialized = false;
}

void gomeaIMS::initialize()
{
	utils::initStartTime();
    problemInstance->initializeRun();
    output = output_statistics_t();

	if( config->AnalyzeFOS )
	{
		prepareFolder(config->folder);
	}
    //initElitistFile(config->folder);

    sharedInformationInstance = new sharedInformation(config->maxArchiveSize);
	isInitialized = true;
	hasTerminated = false;
}

void gomeaIMS::run()
{
	initialize();

	try{
		while(!checkTermination())
		{
			if (numberOfGOMEAs < maximumNumberOfGOMEAs)
				initializeNewGOMEA();

			generationalStepAllGOMEAs();

			numberOfGenerationsIMS++;
		}
	}
	catch( utils::customException const& ){}
	hasTerminated = true;
	writeStatistics(numberOfGOMEAs - 1);
	ezilaitini();
}

void gomeaIMS::runGeneration()
{
	if( !isInitialized )
		initialize();

	if( checkTermination() )
		return;

	try{
		if (currentGOMEAIndex >= numberOfGOMEAs )
		{
			if( numberOfGOMEAs < maximumNumberOfGOMEAs)
				initializeNewGOMEA();
			else
				currentGOMEAIndex = numberOfGOMEAs-1;
		}

		if(!GOMEAs[currentGOMEAIndex]->terminated)
			GOMEAs[currentGOMEAIndex]->terminated = checkTerminationGOMEA(currentGOMEAIndex);
		
		if(!GOMEAs[currentGOMEAIndex]->terminated)
			runGeneration( currentGOMEAIndex );

		if( GOMEAs[currentGOMEAIndex]->numberOfGenerations % IMSsubgenerationFactor == 0 )
			currentGOMEAIndex++;
		else
			currentGOMEAIndex = minimumGOMEAIndex;
	}
	catch( utils::customException const& )
	{
		hasTerminated = true;
		writeStatistics( currentGOMEAIndex );
		ezilaitini();
	}
}

void gomeaIMS::runGeneration( int GOMEAIndex )
{
	GOMEAs[GOMEAIndex]->calculateAverageFitness();

	GOMEAs[GOMEAIndex]->makeOffspring();

	GOMEAs[GOMEAIndex]->copyOffspringToPopulation();

	GOMEAs[GOMEAIndex]->calculateAverageFitness();

	GOMEAs[GOMEAIndex]->numberOfGenerations++;
		
	writeStatistics( GOMEAIndex );
}

bool gomeaIMS::checkTermination()
{
    int i;
	
	if( checkEvaluationLimitTerminationCriterion() )
		hasTerminated = true;

	if( checkTimeLimitTerminationCriterion() )
		hasTerminated = true;

    if (numberOfGOMEAs == maximumNumberOfGOMEAs)
    {
        for (i = 0; i < maximumNumberOfGOMEAs; i++)
        {
            if (!GOMEAs[i]->terminated)
                return false;
        }

		hasTerminated = true;
    }
    
    return hasTerminated;
}

bool gomeaIMS::checkEvaluationLimitTerminationCriterion()
{
	if( !isInitialized )
		return( false );
	if( config->maximumNumberOfEvaluations > 0 && problemInstance->number_of_evaluations > config->maximumNumberOfEvaluations )
		hasTerminated = true;
	return hasTerminated; 
}

bool gomeaIMS::checkTimeLimitTerminationCriterion()
{
	if( !isInitialized )
		return( false );
	if( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSinceStartSeconds() > config->maximumNumberOfSeconds )
		hasTerminated = true;
	return hasTerminated; 
}

double gomeaIMS::getProgressUntilTermination()
{
	double overall_progress = -1.0;

	if( !isInitialized )
		return( -1.0 );

	if( config->maximumNumberOfSeconds > 0 )
	{
		double time_progress = 100.0*utils::getElapsedTimeSinceStartSeconds()/(config->maximumNumberOfSeconds);
		overall_progress = fmax( overall_progress, time_progress );
	}

	if (numberOfGOMEAs == maximumNumberOfGOMEAs && config->maximumNumberOfGenerations > 0 )
	{
		double generational_progress = 100.0*GOMEAs[maximumNumberOfGOMEAs-1]->numberOfGenerations / config->maximumNumberOfGenerations;
		overall_progress = fmax( overall_progress, generational_progress );
	}

	if( config->maximumNumberOfEvaluations > 0 )
	{
		double evaluation_progress = 100.0*problemInstance->number_of_evaluations / config->maximumNumberOfEvaluations;
		overall_progress = fmax( overall_progress, evaluation_progress );
	}

	overall_progress = fmin( overall_progress, 100.0 );

	return overall_progress;
}

void gomeaIMS::initializeNewGOMEA()
{
    #ifdef DEBUG
        cout << "Current number Of GOMEAs is " << numberOfGOMEAs << " | Creating New GOMEA!\n";
    #endif

    Population *newPopulation = NULL;

    if (numberOfGOMEAs == 0)
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGOMEAs, basePopulationSize);
    else
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGOMEAs, 2 * GOMEAs[numberOfGOMEAs-1]->populationSize, GOMEAs[0]->FOSInstance );
    
    GOMEAs.push_back(newPopulation);
    numberOfGOMEAs++;
}

void gomeaIMS::generationalStepAllGOMEAs()
{
    int GOMEAIndexSmallest, GOMEAIndexBiggest;

    GOMEAIndexBiggest    = numberOfGOMEAs - 1;
    GOMEAIndexSmallest = 0;
    while(GOMEAIndexSmallest <= GOMEAIndexBiggest)
    {
        if (!GOMEAs[GOMEAIndexSmallest]->terminated)
            break;

        GOMEAIndexSmallest++;
    }

    GOMEAGenerationalStepAllGOMEAsRecursiveFold(GOMEAIndexSmallest, GOMEAIndexBiggest);
}

void gomeaIMS::GOMEAGenerationalStepAllGOMEAsRecursiveFold(int GOMEAIndexSmallest, int GOMEAIndexBiggest)
{
    int i, GOMEAIndex;

    for(i = 0; i < IMSsubgenerationFactor-1; i++)
    {
        for(GOMEAIndex = GOMEAIndexSmallest; GOMEAIndex <= GOMEAIndexBiggest; GOMEAIndex++)
        {
            if(!GOMEAs[GOMEAIndex]->terminated)
                GOMEAs[GOMEAIndex]->terminated = checkTerminationGOMEA(GOMEAIndex);

            if((!GOMEAs[GOMEAIndex]->terminated) && (GOMEAIndex >= minimumGOMEAIndex))
			{
				runGeneration( GOMEAIndex );
			}
        }

        for(GOMEAIndex = GOMEAIndexSmallest; GOMEAIndex < GOMEAIndexBiggest; GOMEAIndex++)
            GOMEAGenerationalStepAllGOMEAsRecursiveFold(GOMEAIndexSmallest, GOMEAIndex);
    }
}

bool gomeaIMS::checkTerminationGOMEA(int GOMEAIndex)
{
	if( checkTermination() )
		return true;

	if( config->maximumNumberOfGenerations > 0 && (int) GOMEAs[GOMEAIndex]->numberOfGenerations >= config->maximumNumberOfGenerations )
	{
        if( GOMEAIndex == minimumGOMEAIndex )
			minimumGOMEAIndex = GOMEAIndex+1;
		return true;
	}

	if( numberOfGOMEAs > 1 )
	{
		for (int i = GOMEAIndex+1; i < numberOfGOMEAs; i++)
		{        
			if (GOMEAs[i]->averageFitness > GOMEAs[GOMEAIndex]->averageFitness)
			{
				minimumGOMEAIndex = GOMEAIndex+1;
				return true;
			}
		}
	}

	if (!GOMEAs[GOMEAIndex]->allSolutionsAreEqual())
		return false;

	if( GOMEAIndex == minimumGOMEAIndex )
		minimumGOMEAIndex = GOMEAIndex+1;
    return true;
}

void gomeaIMS::writeStatistics( int population_index )
{
	assert( sharedInformationInstance != NULL );
	int key = numberOfStatisticsWrites;
    double evals = problemInstance->number_of_evaluations;
    //double elitist_evals = sharedInformationInstance->elitistSolutionHittingTimeEvaluations;
    double time_s = sharedInformationInstance->elitistSolutionHittingTimeMilliseconds/1000.0;
	double best_fitness = sharedInformationInstance->elitistFitness;
    output.addMetricValue("evaluations",key,evals);
    //output.addMetricValue("elitist_hitting_evaluations",key,elitist_evals);
    output.addMetricValue("time",key,time_s);
    output.addMetricValue("best_obj_val",key,best_fitness);
    output.addMetricValue("population_index",key,population_index);
    output.addMetricValue("generation",key,(int)GOMEAs[population_index]->numberOfGenerations);
    output.addMetricValue("pop_size",key,(int)GOMEAs[population_index]->populationSize);
	numberOfStatisticsWrites++;
}


}}