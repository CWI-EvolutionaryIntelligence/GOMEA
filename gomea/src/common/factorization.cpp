#include "gomea/src/common/factorization.hpp"

namespace gomea{

size_t factorization_t::size()
{
	return factors.size();
}

void factorization_t::addGroup( int var_index )
{
	std::vector<int> vec;
	vec.push_back(var_index);
	addGroup(vec);
}

void factorization_t::addGroup( const std::set<int> &group ) 
{
	std::vector<int> vec;
	for( int x : group )
		vec.push_back(x);
	addGroup(vec);
}

void factorization_t::addGroup( vec_t<int> variables, std::set<int> conditioned_variables )
{
	std::vector<int> cond;
	for( int x : conditioned_variables )
		cond.push_back(x);
	factors.push_back(new cond_factor_t(variables,cond));
}

void factorization_t::initializeCondFactorInteractionGraph( const graph_t &variable_interaction_graph )
{
	assert( factors.size() > 0 );

	// Save clique/factor index that each variable is in for faster lookup
	vec_t<int> clique_index(factors.size());
	for( int i = 0; i < factors.size(); i++ )
	{
		for( int x : factors[i]->variables )
			clique_index[x] = i;
	}

	// Find dependent cliques
	for( int i = 0; i < factors.size(); i++ )
	{
		std::set<int> neighboring_cliques;
		for( int x : factors[i]->variables )
		{
			for( int neighbor : variable_interaction_graph.at(x) )
			{
				int clique_of_neighbor = clique_index[neighbor];
				neighboring_cliques.insert(clique_of_neighbor);
			}
		}
		factorization_interaction_graph.insert({i,neighboring_cliques});
	}
}

/*vec_t<T> factorization_t::samplePartialSolutionConditional( vec_t<int> variables_to_sample, vec_t<cond_factor_t*> factorization, solution_t<T> *parent, const vec_t<solution_t<T>*> &population, int parent_index ) 
{
	if(factorization.size() == 0)
	{
		initializeFactorization();
	}

	//assert( is_conditional );
	if( FOSStructure[FOS_index].size() == number_of_variables )
	{
		// Offspring sample starts as copy of parent
		vec_t<T> sample(number_of_variables);
		for( int i = 0; i < number_of_variables; i++ )
			sample[i] = parent->variables[i];

		for( int i = 0; i < factorization.size(); i++ )
		{
			cond_factor_t *fact = factorization[factorization_order[i]];
			vec_t<T> fact_sample = samplePartialSolutionConditional(sample,population,parent_index);
			// Insert partial sample into offspring sample
			for( int j = 0; j < fact->variables.size(); j++ )
			{
				int var_index = fact->variables[j];
				sample[var_index] = fact_sample[j];
				assert(FOSStructure[FOS_index][var_index] == var_index);
			}
		}
		return sample;
	}
	else
	{
		return factorization[FOS_index]->samplePartialSolutionConditional(parent->variables,population,parent_index);
	}
}

template<class T>
vec_t<T> factorization_t<T>::samplePartialSolutionConditional( solution_t<T> *parent, const vec_t<solution_t<T>*> &population, int parent_index ) 
{
	return samplePartialSolutionConditional(parent->variables, population, parent_index);
}

void sampler_Dt::initializeFactorization()
{
	for( int i = 0; i < variable_groups.size(); i++ )
	{
		cond_factor_Dt *fact = new cond_factor_Dt(variable_groups[i], variable_groups_conditioned_on[i]);
		factorization.push_back( fact );
	}
}

vec_t<char> sampler_Dt::samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index )
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
			!population[donor_index]->hasPartiallyEqualGenotype(parent, this->variables) &&
			isConditionedDonor(population[donor_index], parent) ) 
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
	for( int x : variables )
	{
		sample.push_back(population[donor_index]->variables[x]);
		assert(x == variables[sample.size()-1] ); // check if the order of the indices in 'variables' is the same as the 'variable_groups' when concatenated
	}
	
	assert( sample.size() == variables.size() );

	return sample;
}

bool sampler_Dt::isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent )
{
	return isConditionedDonor(donor_candidate, parent->variables);
}

bool sampler_Dt::isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent )
{
	for (int gene_ind : variables_conditioned_on)
	{
		if (donor_candidate->variables[gene_ind] != parent[gene_ind])
		{
			return false;
		}
	}
	return true;
}

void sampler_Rt::initializeFactorization()
{
	for( int i = 0; i < variable_groups.size(); i++ )
	{
		cond_factor_Rt *fact = new cond_factor_Rt(variable_groups[i], variable_groups_conditioned_on[i]);
		factorization.push_back( fact );
	}
}*/

}