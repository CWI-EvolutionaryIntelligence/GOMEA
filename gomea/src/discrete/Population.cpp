#include "gomea/src/discrete/Population.hpp"

Population::Population(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, FOS_t FOSInstance_ ): 
        config(config_), 
        problemInstance(problemInstance_),
        sharedInformationPointer(sharedInformationPointer_),
        GOMEAIndex(GOMEAIndex_), 
        populationSize(populationSize_)
{
        terminated = false;
        numberOfGenerations = 0;
        averageFitness = 0.0;
        
        population.resize(populationSize);
        offspringPopulation.resize(populationSize);
        noImprovementStretches.resize(populationSize);
        
        vector<int> allGenes(config->numberOfVariables);
        iota(allGenes.begin(), allGenes.end(), 0);

        for (size_t i = 0; i < populationSize; ++i)
        {
            noImprovementStretches[i] = 0;

            population[i] = new Individual(config->numberOfVariables, config->alphabetSize);
            population[i]->randomInit(&config->rng);
            evaluateSolution(population[i], NULL, allGenes, 0.0);
            
            offspringPopulation[i] = new Individual(config->numberOfVariables, config->alphabetSize);
            *offspringPopulation[i] = *population[i];            
        }
			
		if( FOSInstance_ == NULL )
		{
			FOSInstance = createFOSInstance(config->FOSIndex, config->numberOfVariables, config->alphabetSize, config->similarityMeasure, config->maximumFOSSetSize);
		}
		else FOSInstance = FOSInstance_;
        
        #ifdef DEBUG
            cout << "New Population created! Population #" << GOMEAIndex << " PopulationSize:" << populationSize << endl;
            cout << this;
        #endif
}

Population::~Population()
{
    for (size_t i = 0; i < populationSize; ++i)
    {
        delete population[i];
        delete offspringPopulation[i];
    }
}

ostream & operator << (ostream &out, const Population &populationInstance)
{
    for (size_t i = 0; i < populationInstance.populationSize; ++i)
        out << populationInstance.population[i] << " ";
    out << endl;
    return out;
}


void Population::calculateAverageFitness()
{
    averageFitness = 0.0;
    for (size_t i = 0; i < populationSize; ++i)
        averageFitness += population[i]->fitness;
    averageFitness /= populationSize;
}

void Population::copyOffspringToPopulation()
{
    for(size_t i = 0; i < populationSize; i++)
    {
        *population[i] = *offspringPopulation[i];
    }
}

void Population::makeOffspring()
{
	if( config->similarityMeasure == 2 )
	{
		if( FOSInstance->FOSSize() == 0 )
    		FOSInstance->learnFOS(problemInstance->getMIMatrix(), &config->rng);
	}
	else
    	FOSInstance->learnFOS(population, &config->rng);
    
    FOSInstance->setCountersToZero();
    if (config->AnalyzeFOS)
    {
        FOSInstance->writeToFileFOS(config->folder, GOMEAIndex, numberOfGenerations);
    }
 
    generateOffspring();

    if (config->AnalyzeFOS)
        FOSInstance->writeFOSStatistics(config->folder, GOMEAIndex, numberOfGenerations);

}

void Population::generateOffspring()
{
    vector<vector<int> > neighbors;
   
	assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
    vector<int> FOSOrder;
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(FOSOrder, &config->rng); 
	else if( config->useParallelFOSOrder )
		FOSInstance->determineParallelFOSOrder(FOSOrder, problemInstance->getVIG(), &config->rng);

    for(size_t i = 0; i < populationSize; i++)
    {
			if( !config->useParallelFOSOrder && !config->fixFOSOrderForPopulation )
				FOSInstance->shuffleFOS(FOSOrder, &config->rng); 

            Individual backup = *population[i];    
            
            bool solutionHasChanged;
            solutionHasChanged = GOM(i, &backup, FOSOrder );
            
            /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
            if (config->useForcedImprovements)
            {
                if ((!solutionHasChanged) || (noImprovementStretches[i] > (1+(log(populationSize)/log(10)))))
                    FI(i, &backup, FOSOrder);
            }

        if(!(offspringPopulation[i]->fitness > population[i]->fitness))
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
    }
}

bool Population::GOM(size_t offspringIndex, Individual *backup, vector<int> FOSOrder )
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
    
    *offspringPopulation[offspringIndex] = *population[offspringIndex];
            
    vector<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
    {
        int ind = FOSOrder[i];

        if (FOSInstance->FOSElementSize(ind) == 0 || FOSInstance->FOSElementSize(ind) == config->numberOfVariables)
            continue;

        bool donorEqualToOffspring = true;
        size_t indicesTried = 0;

        while (donorEqualToOffspring && indicesTried < donorIndices.size())
        {
            int j = config->rng() % (donorIndices.size() - indicesTried);
            swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
            donorIndex = donorIndices[indicesTried];
            indicesTried++;

            if (donorIndex == offspringIndex)
                continue;

            vector<int> touchedGenes;
            for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];            
                offspringPopulation[offspringIndex]->genotype[variableFromFOS] = population[donorIndex]->genotype[variableFromFOS];
                touchedGenes.push_back(variableFromFOS);

                if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
                    donorEqualToOffspring = false;            
            }

            if (!donorEqualToOffspring)
            {
                evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

                // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
                // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
                if ((!thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness >= backup->fitness)) || 
                        (thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness > backup->fitness)))     
                {             
                    *backup = *offspringPopulation[offspringIndex];
                    
                    solutionHasChanged = true;

                    FOSInstance->improvementCounters[ind]++;
                }
                else
                {
                    *offspringPopulation[offspringIndex] = *backup;
                }

                FOSInstance->usageCounters[ind]++;

            }

            break;
        }
    }
    return solutionHasChanged;
}


bool Population::FI(size_t offspringIndex, Individual *backup, vector<int> FOSOrder )
{
		if( !config->useParallelFOSOrder && !config->fixFOSOrderForPopulation )
			FOSInstance->shuffleFOS(FOSOrder, &config->rng); 

        bool solutionHasChanged = 0;

        for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
        {
            int ind = FOSOrder[i];
            vector<int> touchedGenes;            
            bool donorEqualToOffspring = true;
            for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];
                offspringPopulation[offspringIndex]->genotype[variableFromFOS] = sharedInformationPointer->elitist.genotype[variableFromFOS];
                touchedGenes.push_back(variableFromFOS);
                if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
                    donorEqualToOffspring = false;
            }

            if (!donorEqualToOffspring)
            {
                evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

                if (offspringPopulation[offspringIndex]->fitness > backup->fitness)
                {
                    *backup = *offspringPopulation[offspringIndex];
                    solutionHasChanged = true;
                }
                else
                {
                    *offspringPopulation[offspringIndex] = *backup;
                }
            }
            if (solutionHasChanged)
                break;
        }

        if (!solutionHasChanged)
        {
            *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
        }

        return solutionHasChanged;
}


void Population::evaluateSolution(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{    
    checkTimeLimit();

    /* Do the actual evaluation */
    archiveRecord searchResult;
    
    if (config->saveEvaluations)
        sharedInformationPointer->evaluatedSolutions->checkAlreadyEvaluated(solution->genotype, &searchResult);
    
    if (searchResult.isFound)
        solution->fitness = searchResult.value;
    else
    { 
        //cout << "before eval" << solution -> fitness << endl;
        if (config->usePartialEvaluations && solutionBefore != NULL)
        {
            problemInstance->calculateFitnessPartialEvaluations(solution, solutionBefore, touchedGenes, fitnessBefore);
            sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / config->numberOfVariables;
        }
        else
        {
            problemInstance->calculateFitness(solution);
            sharedInformationPointer->numberOfEvaluations += 1;
        }

        if (config->saveEvaluations)
            sharedInformationPointer->evaluatedSolutions->insertSolution(solution->genotype, solution->fitness);
    }

    updateElitistAndCheckVTR(solution);
}

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && getTime(sharedInformationPointer->startTimeMilliseconds)/1000.0 > config->maximumNumberOfSeconds)
    {
        throw customException("time");
    }
}

void Population::updateElitistAndCheckVTR(Individual *solution)
{
    /* Update elitist solution */
    if (sharedInformationPointer->firstEvaluationEver || (solution->fitness > sharedInformationPointer->elitist.fitness))
    {
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = getTime(sharedInformationPointer->startTimeMilliseconds);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = sharedInformationPointer->numberOfEvaluations;

        sharedInformationPointer->elitist = *solution;
        
        /* Check the VTR */
        if (solution->fitness >= config->vtr)
        {
            writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            cout << "VTR HIT!\n";
            throw customException("vtr");
        }
    
        writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
        if( config->writeElitists )
			writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
    }

    sharedInformationPointer->firstEvaluationEver = false;
}

