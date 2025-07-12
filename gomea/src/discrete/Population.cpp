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

        initializeAndEvaluatePopulation();

        for (size_t i = 0; i < populationSize; ++i)
        {            
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
        if( config->gene_invariant ){
            FOSInstance->MI_truncation_factor = 0.5;
        }
        
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

void Population::initializeAndEvaluatePopulation()
{
    for (size_t i = 0; i < populationSize; ++i)
    {
        noImprovementStretches[i] = 0;
        population[i] = new solution_t<char>(problemInstance->number_of_variables, problemInstance->alphabet_size);
    }

    if( config->gene_invariant ){
        initializePopulationProbabilisticallyComplete();
    }
    else{
        initializePopulationRandomUniform();
    }

    // Evaluate the initial population
    for (size_t i = 0; i < populationSize; ++i){
        problemInstance->evaluate(population[i]);
    }
}

// Initialize the population with a probabilistically complete set of solutions, ensuring that each gene has an equal frequency in the population.
void Population::initializePopulationProbabilisticallyComplete()
{
    // Find nearest (smaller) multiple of alphabet_size
    int permutation_size = problemInstance->alphabet_size * (populationSize / problemInstance->alphabet_size);
    for( int j = 0; j < problemInstance->number_of_variables; j++ )
    {
        // Fill the first part of the population with an equal number of each gene
        vec_t<int> perm = gomea::utils::randomPermutation(permutation_size);
        for( int i = 0; i < permutation_size; i++ ){
            population[i]->variables[j] = perm[i] % problemInstance->alphabet_size;
        }
        // Fill the rest of the population with random (unique) values
        if(permutation_size != populationSize){
            vec_t<int> perm_last_chunk = gomea::utils::randomPermutation(problemInstance->alphabet_size);
            for( int i = permutation_size; i < populationSize; i++ ){
                population[i]->variables[j] = perm_last_chunk[i]; // already in range [0, alphabet_size)
            }
        }
    }
}

void Population::initializePopulationRandomUniform()
{
    // Initialize the population with random uniform solutions
    for (size_t i = 0; i < populationSize; ++i){
        population[i]->randomInit(&gomea::utils::rng);
    }
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

    if( FOSInstance->type == linkage::linkage_model_type::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance, population);
                FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
            }
        }
        else
        {
            FOSInstance->learnLinkageTreeFOS(problemInstance, population);
            FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);
        }
    }

    FOSInstance->setCountersToZero();
    if (config->AnalyzeFOS)
    {
        FOSInstance->writeToFileFOS(config->folder, GOMEAIndex, numberOfGenerations);
    }

    // Initialize the offspring population with copies of the current population
    for( size_t offspringIndex = 0; offspringIndex < populationSize; offspringIndex++ )
        *offspringPopulation[offspringIndex] = *population[offspringIndex];
 
    generateOffspring();

    if (config->AnalyzeFOS)
        FOSInstance->writeFOSStatistics(config->folder, GOMEAIndex, numberOfGenerations);

}

void Population::generateOffspring()
{
    assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS();
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }
    
    if(config->gene_invariant){
        std::vector<std::pair<int, int>> GOM_pairs;
        // Generate random pairs of (individual_index, FOS_index)
        for (int i = 0; i < populationSize; ++i) {
            for (int j = 0; j < FOSInstance->size(); ++j) {
                GOM_pairs.emplace_back(i, j);
            }
        }
        std::shuffle(GOM_pairs.begin(), GOM_pairs.end(), gomea::utils::rng);

        for (auto &[individual_index, FOS_index] : GOM_pairs) {
            GeneInvariantGOM(individual_index, FOS_index);
        }
    }
    else{
        bool *solutionHasChanged = new bool[populationSize];
        bool *isTheElitistSolution = new bool[populationSize];
        for (size_t i = 0; i < populationSize; i++){
            solutionHasChanged[i] = false;
            isTheElitistSolution[i] = (*offspringPopulation[i] == sharedInformationPointer->elitist);
        }

        /* Phase 1: optimal mixing with random donors */
        for (size_t i = 0; i < populationSize; i++){
            if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
                FOSInstance->shuffleFOS();
    
            for (size_t j = 0; j < FOSInstance->size(); j++){
                solutionHasChanged[i] |= GOM(i, FOSInstance->FOSorder[j], isTheElitistSolution[i]);
            }
        }

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            for (size_t i = 0; i < populationSize; i++)
            {
                if ((!solutionHasChanged[i]) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
                    FI(i);
            }
        }
        
        delete[] isTheElitistSolution;
        delete[] solutionHasChanged;
    }

    /* Update or reset no-improvement stretch */
    for (size_t i = 0; i < populationSize; i++)
    {
        if (!problemInstance->betterFitness(offspringPopulation[i], population[i]))
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
    }
}

bool Population::GOM(size_t offspringIndex, int FOS_index, bool isElitistSolution)
{
    int ind = FOS_index; //FOSInstance->FOSorder[i];
    bool solutionHasChanged = false;

    if (FOSInstance->elementSize(ind) == 0 || (int) FOSInstance->elementSize(ind) == problemInstance->number_of_variables){
        return solutionHasChanged;
    }
           
    vec_t<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    bool donorEqualToOffspring = true;
    size_t indicesTried = 0;
    while (donorEqualToOffspring && indicesTried < donorIndices.size())
    {
        int j = gomea::utils::rng() % (donorIndices.size() - indicesTried);
        std::swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
        size_t donorIndex = donorIndices[indicesTried];
        indicesTried++;

        if (donorIndex == offspringIndex)
            continue;

        vec_t<char> donorGenes;
        for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            donorGenes.push_back(population[donorIndex]->variables[variableFromFOS]);
            if (donorGenes[j] != offspringPopulation[offspringIndex]->variables[variableFromFOS])
                donorEqualToOffspring = false;
        }
        partial_solution_t<char> *partial_offspring = new partial_solution_t<char>(donorGenes, FOSInstance->FOSStructure[ind]);

        if (!donorEqualToOffspring)
        {
            problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );
            if( problemInstance->output_frequency == NEW_ELITE && !problemInstance->elitist_was_written )
                writeStatistics();

            // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
            // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
            if ((!isElitistSolution && (partial_offspring->getObjectiveValue() >= offspringPopulation[offspringIndex]->getObjectiveValue())) || 
                    (isElitistSolution && (partial_offspring->getObjectiveValue() > offspringPopulation[offspringIndex]->getObjectiveValue())))     
            {
                offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                
                solutionHasChanged = true;
                updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);

                FOSInstance->improvementCounters[ind]++;
            }

            FOSInstance->usageCounters[ind]++;

        }
        delete partial_offspring;

        break;
    }
    return solutionHasChanged;
}

void Population::GeneInvariantGOM(size_t parent_index, int FOS_index)
{
    if (FOSInstance->elementSize(FOS_index) == 0 || (int) FOSInstance->elementSize(FOS_index) == problemInstance->number_of_variables){
        return;
    }
           
    vec_t<int> mateIndices(populationSize);
    iota(mateIndices.begin(), mateIndices.end(), 0);

    size_t mate_index = 1 - parent_index; // population size == 2
    if( populationSize > 2 ){ // Choose 2 random mates, then pick best
        std::vector<size_t> mates;
        size_t indicesTried = 0;
        while (mates.size() < 2)
        {
            int j = gomea::utils::rng() % (mateIndices.size() - indicesTried);
            std::swap(mateIndices[indicesTried], mateIndices[indicesTried + j]);
            size_t ind = mateIndices[indicesTried];
            indicesTried++;
            if (ind != parent_index)
                mates.push_back(ind);
        }
        assert(mates.size() == 2);
        mate_index = mates[1];
        if( problemInstance->betterFitness(offspringPopulation[mates[0]], offspringPopulation[mates[1]]) )
            mate_index = mates[0];
    }

    // Swap parent and mate indices if mate is at least as good as parent
    if( !problemInstance->betterFitness(offspringPopulation[parent_index], offspringPopulation[mate_index]) )
        std::swap(parent_index, mate_index);

    bool parent_and_mate_equal = true;
    vec_t<char> donor_genes_parent, donor_genes_mate;
    for(size_t j = 0; j < FOSInstance->elementSize(FOS_index); j++)
    {
        int variableFromFOS = FOSInstance->FOSStructure[FOS_index][j];
        donor_genes_parent.push_back(offspringPopulation[parent_index]->variables[variableFromFOS]);
        donor_genes_mate.push_back(offspringPopulation[mate_index]->variables[variableFromFOS]);
        if (donor_genes_parent[j] != donor_genes_mate[j])
            parent_and_mate_equal = false;
    }

    if (!parent_and_mate_equal)
    {
        // Evaluate parent with donor genes from mate
        partial_solution_t<char> *partial_offspring_parent = new partial_solution_t<char>(donor_genes_mate, FOSInstance->FOSStructure[FOS_index]);
        problemInstance->evaluatePartialSolution(offspringPopulation[parent_index], partial_offspring_parent );
        
        if( problemInstance->output_frequency == NEW_ELITE && !problemInstance->elitist_was_written )
            writeStatistics();

        // Accept the change if change to parent is at least equally good (allows random walk in neutral fitness landscape)
        if ( !problemInstance->betterFitness(offspringPopulation[parent_index], partial_offspring_parent) )  // accept
        {
            // Also evaluate the mate with donor genes from parent
            partial_solution_t<char> *partial_offspring_mate = new partial_solution_t<char>(donor_genes_parent, FOSInstance->FOSStructure[FOS_index]);
            problemInstance->evaluatePartialSolution(offspringPopulation[mate_index], partial_offspring_mate );
        
            // Then accept changes to both parent and mate
            offspringPopulation[parent_index]->insertPartialSolution(partial_offspring_parent);
            offspringPopulation[mate_index]->insertPartialSolution(partial_offspring_mate);
            
            updateElitistAndCheckVTR(offspringPopulation[parent_index]);
            updateElitistAndCheckVTR(offspringPopulation[mate_index]);

            FOSInstance->improvementCounters[FOS_index]++;
     
            delete partial_offspring_mate;
        }
        // Do nothing if change is rejected

        FOSInstance->usageCounters[FOS_index]++;

        delete partial_offspring_parent;
    }

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