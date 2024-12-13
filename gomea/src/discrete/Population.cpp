#include "gomea/src/discrete/Population.hpp"

namespace gomea{
namespace discrete{

Population::Population(Config *config_, output_statistics_t *output_, fitness_t<char> *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ ): 
        config(config_),
        output(output_),
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

            population[i] = new solution_t<char>(problemInstance->number_of_variables, problemInstance->alphabet_size);
            population[i]->randomInit(&gomea::utils::rng);
            problemInstance->evaluate(population[i]);
            
            offspringPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, problemInstance->alphabet_size);
            *offspringPopulation[i] = *population[i];
        }
			
		if( config->linkage_config != NULL )
		{
			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
            FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
		}
		else if( FOSInstance_ == NULL )
		{
			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
		}
		else FOSInstance = FOSInstance_;
        
        #ifdef DEBUG
            std::cout << "New Population created! Population #" << GOMEAIndex << " PopulationSize:" << populationSize << endl;
            std::cout << this;
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

std::ostream & operator << (std::ostream &out, const Population &populationInstance)
{
    out << "Generation " << populationInstance.numberOfGenerations << ":" << std::endl;
    for (size_t i = 0; i < populationInstance.populationSize; ++i)
        out << *populationInstance.population[i] << std::endl;
    out << std::endl;
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

double Population::getFitnessMean()
{
	double objective_avg = 0.0;
	for(int i = 0; i < populationSize; i++ )
		objective_avg  += population[i]->getObjectiveValue();
	objective_avg = objective_avg / ((double) populationSize);
	return( objective_avg );
}

double Population::getFitnessVariance()
{
	double objective_avg = getFitnessMean();
	double objective_var = 0.0;
	for(int i = 0; i < populationSize; i++ )
		objective_var  += (population[i]->getObjectiveValue()-objective_avg)*(population[i]->getObjectiveValue()-objective_avg);
	objective_var = objective_var / ((double) populationSize);

	if( objective_var <= 0.0 )
		objective_var = 0.0;
	return( objective_var );
}

double Population::getConstraintValueMean()
{
	double constraint_avg = 0.0;
	for(int i = 0; i < populationSize; i++ )
		constraint_avg  += population[i]->getConstraintValue();
	constraint_avg = constraint_avg / ((double) populationSize);

	return( constraint_avg );
}

double Population::getConstraintValueVariance()
{
	double constraint_avg = getConstraintValueMean();

	double constraint_var = 0.0;
	for(int i = 0; i < populationSize; i++ )
		constraint_var  += (population[i]->getConstraintValue()-constraint_avg)*(population[i]->getConstraintValue()-constraint_avg);
	constraint_var = constraint_var / ((double) populationSize);

	if( constraint_var <= 0.0 )
		constraint_var = 0.0;
	return( constraint_var );
}

solution_t<char> *Population::getBestSolution()
{
	int index_best = 0;
	for(int j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population[j]->getObjectiveValue(), population[j]->getConstraintValue(), population[index_best]->getObjectiveValue(), population[index_best]->getConstraintValue()) )
		{
			index_best = j;
        }
    }
	return( population[index_best] );
}

solution_t<char> *Population::getWorstSolution()
{
	int index_worst = 0;
	for(int j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population[index_worst]->getObjectiveValue(), population[index_worst]->getConstraintValue(), population[j]->getObjectiveValue(), population[j]->getConstraintValue()) )
		{
			index_worst = j;
        }
    }
	return( population[index_worst] );
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
    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < populationSize; ++i)
            updateElitistAndCheckVTR(population[i]);
    }

    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix(FOSInstance->getSimilarityMeasure()), false);
                FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
            }
        }
        else
        {
            FOSInstance->learnLinkageTreeFOS(population, problemInstance->alphabet_size );
            FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);
        }
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

    for (size_t i = 0; i < populationSize; i++)
    {
        if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();

        solution_t<char> backup = *population[i];

        bool solutionHasChanged;
        solutionHasChanged = GOM(i);

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
                FI(i);
        }

        if (!(offspringPopulation[i]->getObjectiveValue() > population[i]->getObjectiveValue()))
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
            std::swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
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
                //problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring, FOSInstance->getDependentSubfunctions(ind) );
                problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );
                if( problemInstance->output_frequency == NEW_ELITE && !problemInstance->elitist_was_written )
                    writeStatistics();

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
            delete partial_offspring;

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
            if( problemInstance->output_frequency == NEW_ELITE && !problemInstance->elitist_was_written )
                writeStatistics();

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

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSeconds(utils::start_time) > config->maximumNumberOfSeconds)
    {
        terminated = true;
        throw utils::terminationException("time");
    }
}

void Population::updateElitistAndCheckVTR(solution_t<char> *solution)
{
    /* Update elitist solution */
    //if (sharedInformationPointer->firstEvaluationEver || (solution->getObjectiveValue() > sharedInformationPointer->elitist.getObjectiveValue()))
    if (sharedInformationPointer->firstEvaluationEver || problemInstance->betterFitness(solution,&sharedInformationPointer->elitist) )
    {
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = utils::getElapsedTimeMilliseconds(utils::start_time);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = problemInstance->number_of_evaluations;

        sharedInformationPointer->elitist = *solution;
		sharedInformationPointer->elitistFitness = solution->getObjectiveValue();
		sharedInformationPointer->elitistConstraintValue = solution->getConstraintValue();
        
        /* Check the VTR */
        if (problemInstance->use_vtr && solution->getObjectiveValue() >= problemInstance->vtr)
        {
            //writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            //writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            //std::cout << "VTR HIT!\n";
            terminated = true;
            throw utils::terminationException("vtr");
        }
    
        //writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
        //if( config->writeElitists )
			//writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
    }

    sharedInformationPointer->firstEvaluationEver = false;
}

void Population::writeStatistics( bool is_final )
{
    assert( sharedInformationPointer != NULL );
	int key = output->number_of_writes;
    double evals = problemInstance->number_of_evaluations;
	double elapsed_time = utils::getElapsedTimeSinceStartSeconds();
	double eval_time = utils::getTimer("eval_time");
    //double elitist_evals = sharedInformationPointer->elitistSolutionHittingTimeEvaluations;
    //double time_s = sharedInformationPointer->elitistSolutionHittingTimeMilliseconds/1000.0;
	//double best_fitness = sharedInformationPointer->elitistFitness;
    output->addGenerationalMetricValue("generation",key,(int)numberOfGenerations);
    output->addGenerationalMetricValue("evaluations",key,evals);
    //output->addMetricValue("elitist_hitting_evaluations",key,elitist_evals);
    output->addGenerationalMetricValue("time",key,elapsed_time);
    output->addGenerationalMetricValue("eval_time",key,eval_time);
    output->addGenerationalMetricValue("population_index",key,(int)GOMEAIndex);
    output->addGenerationalMetricValue("population_size",key,(int)populationSize);
    output->addGenerationalMetricValue("best_obj_val",key,problemInstance->elitist_objective_value);
    output->addGenerationalMetricValue("best_cons_val",key,problemInstance->elitist_constraint_value);
	if( config->generational_solution )
		output->addGenerationalMetricValue("best_solution",key,sharedInformationPointer->elitist.variablesToString());
    //output->addMetricValue("obj_val_avg",key,population_objective_avg);
    //output->addMetricValue("obj_val_var",key,population_objective_var);

	if( is_final ){
		output->addFinalMetricValue("evaluations",evals);
		output->addFinalMetricValue("time",elapsed_time);
		output->addFinalMetricValue("eval_time",eval_time);
		output->addFinalMetricValue("best_solution",sharedInformationPointer->elitist.variablesToString());
		output->addFinalMetricValue("best_obj_val",problemInstance->elitist_objective_value);
		output->addFinalMetricValue("best_cons_val",problemInstance->elitist_constraint_value);
	}

    problemInstance->elitist_was_written = true;
	output->number_of_writes++;
}


}}