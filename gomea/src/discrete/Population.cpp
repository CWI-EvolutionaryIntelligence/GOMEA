#include "gomea/src/discrete/Population.hpp"

namespace gomea{
namespace discrete{

Population::Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ ): 
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
        
        vec_t<int> allGenes(problemInstance->number_of_variables);
        iota(allGenes.begin(), allGenes.end(), 0);

        for (size_t i = 0; i < populationSize; ++i)
        {
            noImprovementStretches[i] = 0;

            population[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
            population[i]->randomInit(&gomea::utils::rng);
            problemInstance->evaluate(population[i]);
            updateElitistAndCheckVTR(population[i]);
            
            offspringPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
            *offspringPopulation[i] = *population[i];
        }
			
		if( config->linkage_config != NULL )
		{
            printf("%s\n",linkage_model_t::getTypeName(config->linkage_config->type).c_str());
			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
		}
		else if( FOSInstance_ == NULL )
		{
			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
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

bool Population::allSolutionsAreEqual()
{
    for (size_t i = 1; i < populationSize; i++)
    {
        for (int j = 0; j < problemInstance->number_of_variables; j++)
        {
            if (population[i]->variables[j] != population[0]->variables[j])
                return false;
        }
    }
    return true;
}

void Population::calculateAverageFitness()
{
    averageFitness = 0.0;
    for (size_t i = 0; i < populationSize; ++i)
        averageFitness += population[i]->getObjectiveValue();
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
    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (config->similarityMeasure == 2)
        {
            if (FOSInstance->size() == 0)
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix() );
        }
        else
            FOSInstance->learnLinkageTreeFOS(population, config->alphabetSize );
    }

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
    vec_t<vec_t<int> > neighbors;
   
	assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(); 
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }

    for(size_t i = 0; i < populationSize; i++)
    {
			if( !config->useParallelFOSOrder && !config->fixFOSOrderForPopulation )
				FOSInstance->shuffleFOS(); 

            solution_t<char> backup = *population[i];
            
            bool solutionHasChanged;
            solutionHasChanged = GOM(i);
            
            /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
            if (config->useForcedImprovements)
            {
                if ((!solutionHasChanged) || (noImprovementStretches[i] > (1+(log(populationSize)/log(10)))))
                    FI(i);
            }

        if(!(offspringPopulation[i]->getObjectiveValue() > population[i]->getObjectiveValue()))
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
    }
}

bool Population::GOM(size_t offspringIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
    
    *offspringPopulation[offspringIndex] = *population[offspringIndex];
            
    vec_t<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];

        if (FOSInstance->elementSize(ind) == 0 || (int) FOSInstance->elementSize(ind) == problemInstance->number_of_variables)
            continue;

        bool donorEqualToOffspring = true;
        size_t indicesTried = 0;

        while (donorEqualToOffspring && indicesTried < donorIndices.size())
        {
            int j = gomea::utils::rng() % (donorIndices.size() - indicesTried);
            swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
            donorIndex = donorIndices[indicesTried];
            indicesTried++;

            if (donorIndex == offspringIndex)
                continue;

            vec_t<char> donorGenes;
            for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];
                //offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                donorGenes.push_back(population[donorIndex]->variables[variableFromFOS]);
                if (donorGenes[j] != offspringPopulation[offspringIndex]->variables[variableFromFOS])
                    donorEqualToOffspring = false;
            }
            partial_solution_t<char> *partial_offspring = new partial_solution_t<char>(donorGenes, FOSInstance->FOSStructure[ind]);

            if (!donorEqualToOffspring)
            {
                //evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->getObjectiveValue());
                problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

                // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
                // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
                if ((!thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() >= offspringPopulation[offspringIndex]->getObjectiveValue())) || 
                        (thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() > offspringPopulation[offspringIndex]->getObjectiveValue())))     
                {
                    offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                    // offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                    //*backup = *offspringPopulation[offspringIndex];
                    
                    solutionHasChanged = true;
                    updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);

                    FOSInstance->improvementCounters[ind]++;
                }

                FOSInstance->usageCounters[ind]++;

            }

            break;
        }
    }
    return solutionHasChanged;
}


bool Population::FI(size_t offspringIndex)
{
    if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
        FOSInstance->shuffleFOS();

    bool solutionHasChanged = 0;

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];
        vec_t<char> touchedGenes = vec_t<char>(FOSInstance->elementSize(ind));
        bool donorEqualToOffspring = true;
        for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            touchedGenes[j] = sharedInformationPointer->elitist.variables[variableFromFOS];
            if (population[offspringIndex]->variables[variableFromFOS] != touchedGenes[j])
                donorEqualToOffspring = false;
        }
        gomea::partial_solution_t<char> *partial_offspring = new gomea::partial_solution_t<char>(touchedGenes, FOSInstance->FOSStructure[ind]);

        if (!donorEqualToOffspring)
        {
            problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

            if (partial_offspring->getObjectiveValue() > offspringPopulation[offspringIndex]->getObjectiveValue() ) 
            {
                offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);
                solutionHasChanged = true;
            }
        }
        delete partial_offspring;
        if (solutionHasChanged)
            break;
    }

    if (!solutionHasChanged)
    {
        *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

/*void Population::evaluateSolution(solution_t<char> *parent, gomea::partial_solution_t<char> *solution ) 
{
    checkTimeLimit();

    //cout << "before eval" << solution -> fitness << endl;
    if (config->usePartialEvaluations && solution != NULL)
    {
        problemInstance->calculateFitnessPartialEvaluations(solution, solutionBefore, touchedGenes, fitnessBefore);
        sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
    }
    else
    {
        problemInstance->calculateFitness(solution);
        sharedInformationPointer->numberOfEvaluations += 1;
    }

    updateElitistAndCheckVTR(solution);
}*/

/*void Population::evaluateSolution(solution_t<char> *solution, solution_t<char> *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore)
{
    checkTimeLimit();

    // Do the actual evaluation
    archiveRecord searchResult;
    
    if (config->saveEvaluations)
        sharedInformationPointer->evaluatedSolutions->checkAlreadyEvaluated(solution->variables, &searchResult);
    
    if (searchResult.isFound)
        solution->setObjectiveValue(searchResult.value);
    else
    { 
        //cout << "before eval" << solution -> fitness << endl;
        if (config->usePartialEvaluations && solutionBefore != NULL)
        {
            assert(0);
            // TODO
            //problemInstance->evaluatePartialSolution(solution, solutionBefore, touchedGenes, fitnessBefore);
            sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
        }
        else
        {
            assert(0);
            // TODO
            //problemInstance->evaluate(solution);
            sharedInformationPointer->numberOfEvaluations += 1;
        }

        if (config->saveEvaluations)
            sharedInformationPointer->evaluatedSolutions->insertSolution(solution->variables, solution->getObjectiveValue());
    }

    updateElitistAndCheckVTR(solution);
}*/

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSeconds(sharedInformationPointer->startTime) > config->maximumNumberOfSeconds)
    {
        terminated = true;
        throw customException("time");
    }
}

void Population::updateElitistAndCheckVTR(solution_t<char> *solution)
{
    /* Update elitist solution */
    if (sharedInformationPointer->firstEvaluationEver || (solution->getObjectiveValue() > sharedInformationPointer->elitist.getObjectiveValue()))
    {
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = utils::getElapsedTimeMilliseconds(sharedInformationPointer->startTime);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = problemInstance->number_of_evaluations;

        sharedInformationPointer->elitist = *solution;
        
        /* Check the VTR */
        if (problemInstance->use_vtr && solution->getObjectiveValue() >= problemInstance->vtr)
        {
            writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            //cout << "VTR HIT!\n";
            terminated = true;
            throw customException("vtr");
        }
    
        writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
        if( config->writeElitists )
			writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
    }

    sharedInformationPointer->firstEvaluationEver = false;
}

}}