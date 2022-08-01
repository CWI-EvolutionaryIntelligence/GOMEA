#include <iostream>
#include <fstream>
using namespace std;

#include "gomeaIMS.hpp"
#include "utils.hpp"

gomeaIMS::gomeaIMS()
{
	return;
}

gomeaIMS::gomeaIMS(Config *config_): config(config_)
{
    maximumNumberOfGOMEAs   = config->maximumNumberOfGOMEAs;
    IMSsubgenerationFactor  = config->IMSsubgenerationFactor;
    basePopulationSize      = config->basePopulationSize;
}

gomeaIMS::~gomeaIMS()
{
	if( isInitialized )
	{
		for (int i = 0; i < numberOfGOMEAs; ++i)
			delete GOMEAs[i];

		delete problemInstance;
		delete sharedInformationInstance;
	}
}

void gomeaIMS::initialize()
{
	startTimer();

	prepareFolder(config->folder);
    initElitistFile(config->folder);
    
	createProblemInstance(config->problemIndex, config->numberOfVariables, config, &problemInstance, config->problemInstancePath);
    #ifdef DEBUG
        cout << "Problem Instance created! Problem number is " << config->problemIndex << endl;
    #endif

    sharedInformationInstance = new sharedInformation(config->maxArchiveSize);
    #ifdef DEBUG
        cout << "Shared Information instance created!\n";
    #endif

	isInitialized = true;
}

void gomeaIMS::run()
{
	if( !isInitialized )
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
	catch( customException const& )
	{}
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
	catch( customException const& )
	{}
}

void gomeaIMS::runGeneration( int GOMEAIndex )
{
	//printf("GOMEA[%d] - pop size %d - generation %d\n",GOMEAIndex,GOMEAs[GOMEAIndex]->populationSize,GOMEAs[GOMEAIndex]->numberOfGenerations);
	
	GOMEAs[GOMEAIndex]->calculateAverageFitness();

	GOMEAs[GOMEAIndex]->makeOffspring();

	GOMEAs[GOMEAIndex]->copyOffspringToPopulation();

	GOMEAs[GOMEAIndex]->calculateAverageFitness();

	GOMEAs[GOMEAIndex]->numberOfGenerations++;
}

bool gomeaIMS::checkTermination()
{
    int i;

	if( checkTimeLimitTerminationCriterion() )
		return true;

    if (numberOfGOMEAs == maximumNumberOfGOMEAs)
    {
        for (i = 0; i < maximumNumberOfGOMEAs; i++)
        {
            if (!GOMEAs[i]->terminated)
                return false;
        }

        return true;
    }
    
    return false;
}

bool gomeaIMS::checkTimeLimitTerminationCriterion()
{
	if( !isInitialized )
		return( false );
	if( config->maximumNumberOfSeconds > 0 && getElapsedTime() > config->maximumNumberOfSeconds )
		return true;
	return false;
}

double gomeaIMS::getProgressUntilTermination()
{
	double overall_progress = -1.0;

	if( !isInitialized )
		return( -1.0 );

	if( config->maximumNumberOfSeconds > 0 )
	{
		double time_progress = 100.0*getElapsedTime()/(config->maximumNumberOfSeconds);
		overall_progress = fmax( overall_progress, time_progress );
	}

	if (numberOfGOMEAs == maximumNumberOfGOMEAs && config->maximumNumberOfGenerations > 0 )
	{
		double generational_progress = 100.0*GOMEAs[maximumNumberOfGOMEAs-1]->numberOfGenerations / config->maximumNumberOfGenerations;
		overall_progress = fmax( overall_progress, generational_progress );
	}

	if( config->maximumNumberOfEvaluations > 0 )
	{
		double evaluation_progress = 100.0*sharedInformationInstance->numberOfEvaluations / config->maximumNumberOfEvaluations;
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
				runGeneration( GOMEAIndex );
        }

        for(GOMEAIndex = GOMEAIndexSmallest; GOMEAIndex < GOMEAIndexBiggest; GOMEAIndex++)
            GOMEAGenerationalStepAllGOMEAsRecursiveFold(GOMEAIndexSmallest, GOMEAIndex);
    }
}

bool gomeaIMS::checkTerminationGOMEA(int GOMEAIndex)
{
	if( config->maximumNumberOfGenerations > 0 && GOMEAs[GOMEAIndex]->numberOfGenerations >= config->maximumNumberOfGenerations )
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

    for (size_t i = 1; i < GOMEAs[GOMEAIndex]->populationSize; i++)
    {
        for (size_t j = 0; j < config->numberOfVariables; j++)
        {
            if (GOMEAs[GOMEAIndex]->population[i]->genotype[j] != GOMEAs[GOMEAIndex]->population[0]->genotype[j])
                return false;
        }
    }

	if( GOMEAIndex == minimumGOMEAIndex )
		minimumGOMEAIndex = GOMEAIndex+1;
    return true;
}
