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

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/real_valued/population.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

population_t::population_t( fitness_t *fitness, int population_size, double lower_init, double upper_init )
{
	this->population_size = population_size;
	this->fitness = fitness;
	initializeDefaultParameters();
	initializeParameterRangeBounds( lower_init, upper_init );
}

population_t::~population_t()
{
	for(int j = 0; j < population_size; j++ )
		delete( individuals[j] );
	free( individuals );
	
	free( selection );
	free( ranks );

	free( mean_shift_vector );
	free( prev_mean_vector );

	free( individual_NIS );

	free( lower_init_ranges );
	free( upper_init_ranges );
	
	for(int j = 0; j < linkage_model->size(); j++ )
		free( sampled_solutions[j] );
	free( sampled_solutions );
}

void population_t::initialize()
{
	initializeNewPopulationMemory();

	initializePopulationAndFitnessValues();
}

void population_t::runGeneration()
{
	if( population_terminated )
		return;

	makeSelection();

	updateElitist();

	estimateDistribution();

	copyBestSolutionsToPopulation();

	generateAndEvaluateNewSolutions();

	number_of_generations++;
}

void population_t::updateElitist()
{
	solution_t<double> *best_so_far = individuals[0];
	for( int i = 1; i < population_size; i++ )
	{
		if( fitness->betterFitness( individuals[i], best_so_far ) )
		   best_so_far = individuals[i];	
	}
	objective_value_elitist = best_so_far->getObjectiveValue();
	constraint_value_elitist = best_so_far->getConstraintValue();
}

void population_t::makeSelection()
{
	computeRanks();
	int *sorted = mergeSort( ranks, population_size );

	if( ranks[sorted[selection_size-1]] == 0 )
		makeSelectionUsingDiversityOnRank0();
	else
	{
		for(int i = 0; i < selection_size; i++ )
			selection[i] = individuals[sorted[i]];
	}

	free( sorted );
}
		
void population_t::makeSelectionUsingDiversityOnRank0()
{
	int number_of_rank0_solutions = 0;
	for(int i = 0; i < population_size; i++ )
	{
		if( ranks[i] == 0 )
			number_of_rank0_solutions++;
	}

	int *preselection_indices = (int *) Malloc( number_of_rank0_solutions*sizeof( int ) );
	int k = 0;
	for(int i = 0; i < population_size; i++ )
	{
		if( ranks[i] == 0 )
		{
			preselection_indices[k] = i;
			k++;
		}
	}

	int index_of_farthest = 0;
	double distance_of_farthest = individuals[preselection_indices[0]]->getObjectiveValue();
	for(int i = 1; i < number_of_rank0_solutions; i++ )
	{
		if( individuals[preselection_indices[i]]->getObjectiveValue() > distance_of_farthest )
		{
			index_of_farthest    = i;
			distance_of_farthest = individuals[preselection_indices[i]]->getObjectiveValue();
		}
	}

	int number_selected_so_far = 0;
	int *selection_indices = (int *) Malloc( selection_size*sizeof( int ) );
	selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
	preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
	number_of_rank0_solutions--;
	number_selected_so_far++;

	double *nn_distances = (double *) Malloc( number_of_rank0_solutions*sizeof( double ) );
	for(int i = 0; i < number_of_rank0_solutions; i++ )
		nn_distances[i] = distanceEuclidean( individuals[preselection_indices[i]]->variables, individuals[selection_indices[number_selected_so_far-1]]->variables ); 

	while( number_selected_so_far < selection_size )
	{
		index_of_farthest    = 0;
		distance_of_farthest = nn_distances[0];
		for(int i = 1; i < number_of_rank0_solutions; i++ )
		{
			if( nn_distances[i] > distance_of_farthest )
			{
				index_of_farthest    = i;
				distance_of_farthest = nn_distances[i];
			}
		}

		selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
		preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
		nn_distances[index_of_farthest]           = nn_distances[number_of_rank0_solutions-1];
		number_of_rank0_solutions--;
		number_selected_so_far++;

		for(int i = 0; i < number_of_rank0_solutions; i++ )
		{
			double value = distanceEuclidean( individuals[preselection_indices[i]]->variables, individuals[selection_indices[number_selected_so_far-1]]->variables );
			if( value < nn_distances[i] )
				nn_distances[i] = value;
		}
	}

	for(int i = 0; i < selection_size; i++ )
		selection[i] = individuals[selection_indices[i]];

	free( nn_distances );
	free( selection_indices );
	free( preselection_indices );
}

void population_t::computeRanks()
{
	vec_t<int> sorted(population_size);
	for(int i = 0; i < population_size; i++ )
		sorted[i] = i;
	std::sort( sorted.begin(), sorted.end(), [&](int x, int y){return fitness->betterFitness(individuals[x],individuals[y]);}); 

	int rank = 0;
	ranks[sorted[0]] = rank;
	for(int i = 1; i < population_size; i++ )
	{
		if( individuals[sorted[i]]->getObjectiveValue() != individuals[sorted[i-1]]->getObjectiveValue() )
			rank++;

		ranks[sorted[i]] = rank;
	}
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void population_t::estimateDistribution()
{
	if( !linkage_model->is_conditional )
	{
		/*if( fitness->hasVariableInteractionGraph() )
			linkage_model->randomizeOrder(fitness->variable_interaction_graph);
		else*/
		linkage_model->shuffleFOS();
	}
	for( int i = 0; i < linkage_model->size(); i++ )
		estimateDistribution(i);
	updateAMSMeans();
}

void population_t::estimateDistribution( int FOS_index )
{
	linkage_model->estimateDistribution( FOS_index, selection, selection_size );
}

double population_t::estimateMean( int var )
{
	double mean = 0.0;
	for(int j = 0; j < selection_size; j++ )
		mean += selection[j]->variables[var];
	mean /= (double) selection_size;
	return( mean );
}

void population_t::updateAMSMeans()
{
	for(int i = 0; i < fitness->number_of_variables; i++ )
	{
		double new_mean = estimateMean(i);
		if( number_of_generations > 0 )
			mean_shift_vector[i] = new_mean - prev_mean_vector[i];

		prev_mean_vector[i] = new_mean;
	}
}

void population_t::copyBestSolutionsToPopulation()
{
	assert( num_elitists_to_copy == 1 ); // elitists to be copied should first be copied to avoid overwriting them beforehand
	for( int i = 0; i < num_elitists_to_copy; i++ )
	{
		for(int k = 0; k < fitness->number_of_variables; k++ )
			individuals[i]->variables[k] = selection[i]->variables[k];

		individuals[i]->setObjectiveValue(selection[i]->getObjectiveValue());
		individuals[i]->setConstraintValue(selection[i]->getConstraintValue());
		individuals[i]->setFitnessBuffers(selection[i]->fitness_buffers);
	}
}

void population_t::getBestInPopulation( int *individual_index )
{
	*individual_index = 0;
	for(int i = 0; i < population_size; i++ )
		if( fitness->betterFitness(individuals[i]->getObjectiveValue(), individuals[i]->getConstraintValue(), individuals[*individual_index]->getObjectiveValue(), individuals[*individual_index]->getConstraintValue()))
			*individual_index = i;
}

void population_t::evaluateCompletePopulation()
{
	for(int j = 0; j < population_size; j++ )
		fitness->evaluate( individuals[j] );        
}

void population_t::generateAndEvaluateNewSolutions()
{
	if( !fitness->black_box_optimization && (number_of_generations+1) % 50 == 0 )
		evaluateCompletePopulation();

	short *individual_improved = (short *) Malloc( population_size*sizeof( short ) );
	for(int k = num_elitists_to_copy; k < population_size; k++ )
		individual_improved[k] = 0;

	double alpha_AMS = 0.5*tau*(((double) population_size)/((double) (population_size-1)));
	int number_of_AMS_solutions = (int) (alpha_AMS*(population_size-1));
	/*for(int j = 0; j < linkage_model->size(); j++ )
	{
		samples_drawn_from_normal[j] = 0;
		out_of_bounds_draws[j]       = 0;
	}*/ // BLA - resets every generation in distribution class?

	linkage_model->shuffleFOS();
	for(int g = 0; g < linkage_model->size(); g++ )
	{
		int FOS_index = linkage_model->FOSorder[g];

		if( selection_during_gom )
		{
			makeSelection();
			estimateDistribution(FOS_index);
		}
		if( update_elitist_during_gom )
			updateElitist();

		for(int k = num_elitists_to_copy; k < population_size; k++ )
			sampled_solutions[FOS_index][k] = linkage_model->generatePartialSolution( FOS_index, individuals[k] );
		
		if( number_of_generations > 0 )
		{
			for(int k = num_elitists_to_copy; k <= number_of_AMS_solutions; k++ )
				applyPartialAMS( sampled_solutions[FOS_index][k], linkage_model->getDistributionMultiplier(FOS_index) );
		}
		
		for(int k = num_elitists_to_copy; k < population_size; k++ )
			fitness->evaluatePartialSolution( individuals[k], sampled_solutions[FOS_index][k] );
		
		int num_improvements = 0;
		short *accept_improvement = (short*) Malloc( population_size * sizeof(short) );
		for(int k = num_elitists_to_copy; k < population_size; k++ )
		{
			accept_improvement[k] = checkForImprovement( individuals[k], sampled_solutions[FOS_index][k] );
			individual_improved[k] = accept_improvement[k];
			if( accept_improvement[k] ) num_improvements++;
		}

		for(int k = num_elitists_to_copy; k < population_size; k++ )
		{
			if( accept_improvement[k] || randomRealUniform01() < linkage_model->getAcceptanceRate() )
			{
				sampled_solutions[FOS_index][k]->is_accepted = 1;
				insertImprovement( individuals[k], sampled_solutions[FOS_index][k] );
			}
			else
			{
				for(int i = 0; i < sampled_solutions[FOS_index][k]->getNumberOfTouchedVariables(); i++ )
				{
					int ind = sampled_solutions[FOS_index][k]->touched_indices[i];
					sampled_solutions[FOS_index][k]->touched_variables[i] = individuals[k]->variables[ind];
				}
				sampled_solutions[FOS_index][k]->setObjectiveValue(individuals[k]->getObjectiveValue());
				sampled_solutions[FOS_index][k]->setConstraintValue(individuals[k]->getConstraintValue());
				sampled_solutions[FOS_index][k]->setFitnessBuffers(individuals[k]->fitness_buffers);
			}	

			if( fitness->betterFitness( sampled_solutions[FOS_index][k]->getObjectiveValue(), sampled_solutions[FOS_index][k]->getConstraintValue(), objective_value_elitist, constraint_value_elitist ) )
				sampled_solutions[FOS_index][k]->improves_elitist = 1;
		}
		free( accept_improvement );

		linkage_model->adaptDistributionMultiplier( FOS_index, &sampled_solutions[FOS_index][num_elitists_to_copy], population_size-num_elitists_to_copy );
	}

	for(int g = 0; g < linkage_model->size(); g++ )
		for(int k = num_elitists_to_copy; k < population_size; k++ )
			delete( sampled_solutions[g][k] );
	
	if( number_of_generations > 0 )
	{
		for(int k = num_elitists_to_copy; k <= number_of_AMS_solutions; k++ )
			individual_improved[k] |= applyAMS(k);
	}
		
	short generational_improvement = 0;
	for(int i = num_elitists_to_copy; i < population_size; i++ )
	{
		if( !individual_improved[i] )
			individual_NIS[i]++;
		else
		{
			individual_NIS[i] = 0;
			generational_improvement = 1;
		}
	}

	int best_individual_index;
	getBestInPopulation( &best_individual_index );
	for(int k = num_elitists_to_copy; k < population_size; k++ )
		if( individual_NIS[k] > maximum_no_improvement_stretch )
			applyForcedImprovements( k, best_individual_index );

	if( generational_improvement )
		linkage_model->no_improvement_stretch = 0;
	else
	{
		short all_multipliers_leq_one = 1;
		for(int j = 0; j < linkage_model->size(); j++ )
			if( linkage_model->getDistributionMultiplier(j) > 1.0 )
			{
				all_multipliers_leq_one = 0;
				break;
			}

		if( all_multipliers_leq_one )
			linkage_model->no_improvement_stretch++;
	}

	free( individual_improved );
}

void population_t::applyPartialAMS( partial_solution_t<double> *solution, double cmul )
{
	short out_of_range = 1;
	double shrink_factor = 2;
	double *result = (double*) Malloc( solution->getNumberOfTouchedVariables() * sizeof(double) );
	while( (out_of_range == 1) && (shrink_factor > 1e-10) )
	{
		shrink_factor *= 0.5;
		out_of_range   = 0;
		for(int m = 0; m < solution->getNumberOfTouchedVariables(); m++ )
		{
			int im = solution->touched_indices[m];
			result[m] = solution->touched_variables[m] + shrink_factor * delta_AMS * cmul * (mean_shift_vector[im]);
			if( !fitness->isParameterInRangeBounds( result[m], im ) )
			{
				out_of_range = 1;
				break;
			}
		}
	}
	if( !out_of_range )
	{
		for(int m = 0; m < solution->getNumberOfTouchedVariables(); m++ )
		{
			solution->touched_variables[m] = result[m];
		}
	}
	free( result );
}
		
short population_t::checkForImprovement( solution_t<double> *solution, partial_solution_t<double> *part )
{
	return( fitness->betterFitness( part->getObjectiveValue(), part->getConstraintValue(), solution->getObjectiveValue(), solution->getConstraintValue() ) );
}

void population_t::insertImprovement( solution_t<double> *solution, partial_solution_t<double> *part )
{
	for(int j = 0; j < part->getNumberOfTouchedVariables(); j++ )
	{
		int ind = part->touched_indices[j];
		solution->variables[ind] = part->touched_variables[j];
	}
	//solution->buffer += part->buffer;
	solution->setObjectiveValue(part->getObjectiveValue());
	solution->setConstraintValue(part->getConstraintValue());
	solution->setFitnessBuffers(part->fitness_buffers);
}

short population_t::applyAMS( int individual_index )
{
	short out_of_range  = 1;
	short improvement   = 0;
	double delta_AMS     = 2;
	double shrink_factor = 2;
	solution_t<double> *solution_AMS = new solution_t<double>(fitness->number_of_variables);
	while( (out_of_range == 1) && (shrink_factor > 1e-10) )
	{
		shrink_factor *= 0.5;
		out_of_range   = 0;
		for(int m = 0; m < fitness->number_of_variables; m++ )
		{
			solution_AMS->variables[m] = individuals[individual_index]->variables[m] + shrink_factor*delta_AMS*(mean_shift_vector[m]);
			if( !fitness->isParameterInRangeBounds( solution_AMS->variables[m], m ) )
			{
				out_of_range = 1;
				break;
			}
		}
	}
	if( !out_of_range )
	{
		short improvement;
		fitness->evaluate( solution_AMS );
		improvement = fitness->betterFitness(solution_AMS->getObjectiveValue(), solution_AMS->getConstraintValue(), individuals[individual_index]->getObjectiveValue(), individuals[individual_index]->getConstraintValue()); 
		//if( improvement )
		if( randomRealUniform01() < linkage_model->getAcceptanceRate() || improvement ) // BLA
		{
			individuals[individual_index]->setObjectiveValue(solution_AMS->getObjectiveValue());
			individuals[individual_index]->setConstraintValue(solution_AMS->getConstraintValue());
			individuals[individual_index]->setFitnessBuffers(solution_AMS->fitness_buffers);
			for(int m = 0; m < fitness->number_of_variables; m++ )
				individuals[individual_index]->variables[m] = solution_AMS->variables[m];
			improvement = 1;
		}
	}
	delete( solution_AMS );

	return( improvement );
}

void population_t::applyForcedImprovements( int individual_index, int donor_index )
{
	short improvement = 0;
	double alpha = 1.0;

	while( alpha >= 0.01 )
	{
		alpha *= 0.5;
		for(int io = 0; io < linkage_model->size(); io++ )
		{
			int i = linkage_model->FOSorder[io];
			vec_t<int> touched_indices = linkage_model->sets[i];
			int num_touched_indices = linkage_model->elementSize(i);

			vec_t<double> FI_vars = vec_t<double>(num_touched_indices);
			for(int j = 0; j < num_touched_indices; j++ )
			{
				FI_vars[j] = alpha*individuals[individual_index]->variables[touched_indices[j]] + (1-alpha)*individuals[donor_index]->variables[touched_indices[j]];
			}
			partial_solution_t<double> *FI_solution = new partial_solution_t<double>( FI_vars, touched_indices );
			fitness->evaluatePartialSolution( individuals[individual_index], FI_solution );
			improvement = fitness->betterFitness( FI_solution->getObjectiveValue(), FI_solution->getConstraintValue(), individuals[individual_index]->getObjectiveValue(), individuals[individual_index]->getConstraintValue() );

			if( improvement )
			{
				for(int j = 0; j < num_touched_indices; j++ )
					individuals[individual_index]->variables[touched_indices[j]] = FI_solution->touched_variables[j];
				//individuals[individual_index]->buffer += FI_solution->buffer;
				individuals[individual_index]->setObjectiveValue(FI_solution->getObjectiveValue());
				individuals[individual_index]->setConstraintValue(FI_solution->getConstraintValue());
				individuals[individual_index]->setFitnessBuffers(FI_solution->fitness_buffers);
			}
			delete FI_solution;

			if( improvement )
				break;
		}
		if( improvement )
			break;
	}

	if( !improvement )
	{
		for(int i = 0; i < fitness->number_of_variables; i++ )
			individuals[individual_index]->variables[i] = individuals[donor_index]->variables[i];
		individuals[individual_index]->setObjectiveValue(individuals[donor_index]->getObjectiveValue());
		individuals[individual_index]->setConstraintValue(individuals[donor_index]->getConstraintValue());
		individuals[individual_index]->setFitnessBuffers(individuals[donor_index]->fitness_buffers);
	}
}

double population_t::getFitnessMean()
{
	double objective_avg = 0.0;
	for(int i = 0; i < population_size; i++ )
		objective_avg  += individuals[i]->getObjectiveValue();
	objective_avg = objective_avg / ((double) population_size);

	return( objective_avg );
}

double population_t::getFitnessVariance()
{
	double objective_avg = getFitnessMean();

	double objective_var = 0.0;
	for(int i = 0; i < population_size; i++ )
		objective_var  += (individuals[i]->getObjectiveValue()-objective_avg)*(individuals[i]->getObjectiveValue()-objective_avg);
	objective_var = objective_var / ((double) population_size);

	if( objective_var <= 0.0 )
		objective_var = 0.0;
	return( objective_var );
}

void population_t::initializeDefaultParameters()
{
	eta_cov = 1.0;
	tau = 0.35;
	st_dev_ratio_threshold = 1.0;
	distribution_multiplier_decrease = 0.9;
	distribution_multiplier_increase = 1.0/distribution_multiplier_decrease;
	maximum_no_improvement_stretch = 100;
	delta_AMS = 2.0;
	selection_size = (int) (tau * population_size);
}

void population_t::initializeNewPopulationMemory()
{
	individuals = (solution_t<double>**) Malloc( population_size*sizeof( solution_t<double>* ) );
	for(int j = 0; j < population_size; j++ )
		individuals[j] = new solution_t<double>(fitness->number_of_variables);

	ranks = (double *) Malloc( population_size*sizeof( double ) );

	selection = (solution_t<double>**) Malloc( selection_size*sizeof( solution_t<double> * ) );

	mean_shift_vector = (double *) Malloc( fitness->number_of_variables*sizeof( double ) );
	prev_mean_vector = (double *) Malloc( fitness->number_of_variables*sizeof( double ) );

	individual_NIS = (int*) Malloc( population_size*sizeof(int));

	initializeFOS();

	population_terminated = 0;

	number_of_generations = 0;
}

/**
 * Initializes the linkage tree
 */
void population_t::initializeFOS()
{
	linkage_model_rv_pt new_FOS = NULL;

	int max_clique_size;
	bool include_cliques_as_fos_elements, include_full_fos_element; 
	if( FOS_element_size > 0 )
	{
		//new_FOS = new linkage_model_rv_t( fitness->number_of_variables, FOS_element_size );
		new_FOS = linkage_model_rv_t::marginal_product_model( fitness->number_of_variables, FOS_element_size );
	}
	else if( FOS_element_size == -3 )
	{
		vec_t<vec_t<int>> FOS;
		vec_t<int> full_fos_element;
		for( int i = 0; i < fitness->number_of_variables; i++ )
		{
			vec_t<int> e;
			e.push_back(i);
			full_fos_element.push_back(i);
			FOS.push_back(e);
		}
		FOS.push_back(full_fos_element);
		//new_FOS = new linkage_model_rv_t( fitness->number_of_variables );
		new_FOS = linkage_model_rv_t::custom_fos( fitness->number_of_variables, FOS );
	}
	else if( FOS_element_size == -4 )
	{
		static_linkage_tree = 1;
		FOS_element_ub = 100;
		double **cov = NULL;
		new_FOS = linkage_model_rv_t::linkage_tree( fitness->number_of_variables );
		// TODO
	}
	else if( FOS_element_size <= -10 )
	{
		int id = -1 * FOS_element_size;
		id /= 10;
		include_full_fos_element = (id%10) == 1;
		id /= 10;
		include_cliques_as_fos_elements = (id%10) == 1;
		id /= 10;
		max_clique_size = id;
		new_FOS = linkage_model_rv_t::conditional(fitness->number_of_variables,fitness->variable_interaction_graph,max_clique_size,include_cliques_as_fos_elements,include_full_fos_element);
	}
	else
	{
		assert(0);
	}

	//new_FOS->print();
	linkage_model = new_FOS;
}

/**
 * Initializes the parameter range bounds.
 */
void population_t::initializeParameterRangeBounds( double lower_user_range, double upper_user_range )
{
	lower_init_ranges  = (double *) Malloc( fitness->number_of_variables*sizeof( double ) );
	upper_init_ranges  = (double *) Malloc( fitness->number_of_variables*sizeof( double ) );

	for(int i = 0; i < fitness->number_of_variables; i++ )
	{
		lower_init_ranges[i] = lower_user_range;
		if( lower_user_range < fitness->getLowerRangeBound(i) )
			lower_init_ranges[i] = fitness->getLowerRangeBound(i);
		if( lower_user_range > fitness->getUpperRangeBound(i) )
			lower_init_ranges[i] = fitness->getLowerRangeBound(i);

		upper_init_ranges[i] = upper_user_range;
		if( upper_user_range > fitness->getUpperRangeBound(i) )
			upper_init_ranges[i] = fitness->getUpperRangeBound(i);
		if( upper_user_range < fitness->getLowerRangeBound(i) )
			upper_init_ranges[i] = fitness->getUpperRangeBound(i);
	}
}


void population_t::initializeProblem( int problem_index )
{
	this->problem_index = problem_index;
	switch( problem_index )
	{
		default: break;
	}
}

void population_t::initializePopulationAndFitnessValues()
{
	for(int j = 0; j < population_size; j++ )
	{
		individual_NIS[j] = 0;
		for(int k = 0; k < fitness->number_of_variables; k++ )
			individuals[j]->variables[k] = lower_init_ranges[k] + (upper_init_ranges[k] - lower_init_ranges[k])*randomRealUniform01();

		fitness->evaluate( individuals[j] );
	}

	sampled_solutions = (partial_solution_t<double>***) Malloc( linkage_model->size() * sizeof(partial_solution_t<double>**) );
	for(int j = 0; j < linkage_model->size(); j++ )
		sampled_solutions[j] = (partial_solution_t<double>**) Malloc( population_size * sizeof(partial_solution_t<double>*) );
}

}}
