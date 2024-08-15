#include "gomea/src/discrete/sampler.hpp"

namespace gomea{

vec_t<char> sampler_Dt::sampleSolution( linkage_model_pt linkage_model, int FOS_index, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index )
{
	if( linkage_model->FOSStructure[FOS_index].size() == fitness_function->number_of_variables && linkage_model->is_conditional)
	{
		return sampleFullSolutionConditional(linkage_model->factorization, parent, population, parent_index);
	}
	else
	{
		assert( FOS_index < linkage_model->factorization->size() );
		assert( linkage_model->factorization->factors[FOS_index].size() == linkage_model->elementSize(FOS_index) );
		for(int i = 0; i < linkage_model->elementSize(FOS_index); i++ )
			assert( linkage_model->factorization->factors[FOS_index].variables[i] == linkage_model->FOSStructure[FOS_index][i] );
		return samplePartialSolution(linkage_model->factorization->factors[FOS_index], parent, population, parent_index);
	}
}

vec_t<char> sampler_Dt::sampleFullSolutionConditional( factorization_t *factorization, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	/*if(factorization.size() == 0)
	{
		initializeFactorization();
	}*/

	// Offspring sample starts as copy of parent
	vec_t<char> sample(number_of_variables);
	for( int i = 0; i < number_of_variables; i++ )
		sample[i] = parent->variables[i];

	for( int i = 0; i < factorization->size(); i++ )
	{
		cond_factor_t *factor = factorization->factors[factorization->order[i]];
		vec_t<char> fact_sample = samplePartialSolution(factor,sample,population,parent_index);
		// Insert partial sample into offspring sample
		for( int j = 0; j < factor->variables.size(); j++ )
		{
			int var_index = factor->variables[j];
			sample[var_index] = fact_sample[j];
		}
	}
	return sample;
}

vec_t<char> sampler_Dt::samplePartialSolution( cond_factor_t *factor, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	return samplePartialSolution(factor, parent->variables, population, parent_index);
}

vec_t<char> sampler_Dt::samplePartialSolution( cond_factor_t *factor, const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index )
{
	vec_t<char> sample; // In the same order as 'variables'
	vec_t<int> donorIndices(population.size());
    std::iota(donorIndices.begin(), donorIndices.end(), 0);

	std::shuffle(donorIndices.begin(),donorIndices.end(),utils::rng);
	// Find a donor that is different from the parent and is a suitable conditioned donor for the parent
	int donors_tried = 0;
	int donor_index = -1; 
	while( donors_tried < population.size() )
	{
		donor_index = donorIndices[donors_tried];
		if( donor_index != parent_index &&
			!population[donor_index]->hasPartiallyEqualGenotype(parent, factor->variables) &&
			isConditionedDonor(factor, population[donor_index], parent) ) 
		{
			break; // Suitable donor found
		}
		donors_tried++;
	}
	if( donors_tried == population.size() ) // No suitable donor was found - take random solution from population
	{
		donor_index = donorIndices[0]; // First index is a random solution, as donorIndices was shuffled
	}
		
	// Insert into parent
	for( int x : factor->variables )
	{
		sample.push_back(population[donor_index]->variables[x]);
		assert(x == factor->variables[sample.size()-1] ); // check if the order of the indices in 'variables' is the same as the 'variable_groups' when concatenated
	}
	
	assert( sample.size() == factor->variables.size() );

	return sample;
}

bool sampler_Dt::isConditionedDonor( cond_factor_t *factor, solution_t<char> *donor_candidate, solution_t<char> *parent )
{
	return isConditionedDonor(factor, donor_candidate, parent->variables);
}

bool sampler_Dt::isConditionedDonor( cond_factor_t *factor, solution_t<char> *donor_candidate, const vec_t<char> &parent )
{
	for (int gene_ind : factor->variables_conditioned_on)
	{
		if (donor_candidate->variables[gene_ind] != parent[gene_ind])
		{
			return false;
		}
	}
	return true;
}

}