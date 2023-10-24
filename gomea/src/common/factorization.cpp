#include "gomea/src/common/factorization.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{

factorization_t::factorization_t() 
{
}

factorization_t::factorization_t( const vec_t<int> &variables )
{
	addGroupOfVariables(variables);
}

factorization_t::factorization_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

factorization_t::factorization_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

void factorization_t::addGroupOfVariables( vec_t<int> indices, vec_t<int> indices_cond )
{
	std::sort(indices.begin(),indices.end());
	std::sort(indices_cond.begin(),indices_cond.end());
	for( int i : indices )
		variables.push_back(i);
	std::sort(variables.begin(),variables.end());
	for( int i : indices_cond )
		variables_conditioned_on.push_back(i);
	sampling_order.push_back(sampling_order.size());
}

void factorization_t::addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond )
{
	vec_t<int> cond;
	for( int i : indices_cond )
		cond.push_back(i);
	addGroupOfVariables(indices,cond);
}
			
void factorization_t::addGroupOfVariables( int index, const vec_t<int> &indices_cond )
{
	vec_t<int> indices;
	indices.push_back(index);
	addGroupOfVariables(indices, indices_cond);
}
			
void factorization_t::addGroupOfVariables( int index, int index_cond )
{
	vec_t<int> indices, indices_cond;
	indices.push_back(index);
	indices_cond.push_back(index_cond);
	addGroupOfVariables(indices, indices_cond);
}

void factorization_t::updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited ) 
{
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	
	variables_conditioned_on.clear();

	// Add FOS element of all nodes in clique, conditioned on dependent, already visited variables
	std::set<int> cond;
	for( int v : variables )
		visited[v] = IN_CLIQUE;
	for( int v : variables )
	{
		for( int x : variable_interaction_graph.at(v) )
		{
			if( visited[x] == IS_VISITED )
				cond.insert(x);
		}
	}
	for( int v : variables )
		visited[v] = IS_VISITED;
	
	for( int x : cond )
		variables_conditioned_on.push_back(x);
}

void factorization_t::setOrder( const vec_t<int> &sampling_order ) 
{
	this->sampling_order = sampling_order;
}

void factorization_t::initializeFrequencyTables( const vec_t<solution_t<char>*> &population ) 
{
	vec_t<vec_t<int>> donor_list;
	// TODO
	assert(0);
}

bool factorization_t::isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent )
{
	return isConditionedDonor(donor_candidate, parent->variables);
}

bool factorization_t::isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent )
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

vec_t<char> factorization_t::samplePartialSolutionConditional( solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	return samplePartialSolutionConditional(parent->variables, population, parent_index);
}

vec_t<char> factorization_t:: samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	vec_t<char> sample; // In the same order as 'variables'
	vec_t<int> donorIndices(population.size());
    std::iota(donorIndices.begin(), donorIndices.end(), 0);
	/*for( int i = 0; i < variables.size(); i++ )
		printf("%d ",variables[i]);
	printf("\n");*/

	//int og = sampling_order[i];
	std::shuffle(donorIndices.begin(),donorIndices.end(),utils::rng);
	int donors_tried = 0;
	// Find a donor that is different from the parent and is a suitable conditioned donor for the parent
	while( donors_tried < population.size() )
	{
		int donor_index = donorIndices[donors_tried];
		if( donor_index != parent_index &&
			!population[donor_index]->hasPartiallyEqualGenotype(parent, this->variables) &&
			isConditionedDonor(population[donor_index], parent) ) 
		{
			for( int x : variables )
				sample.push_back(population[donor_index]->variables[x]);
			//printf("Donor found.\n");
			break;
		}
		donors_tried++;
	}
	if( donors_tried == population.size() ) // No suitable donor was found - take random solution from population
	{
		/*if( config->verbose )
		{
			printf("Random donor.\n");
			printf("P: ");
			for( int j = 0; j < parent.size(); j++ )
			{
				bool is_donated = false;
				for( int x : variables )
					if( j == x ) is_donated = true;
				
				if( is_donated ) printf("*");
				printf("%d", parent[j]);
				if( is_donated ) printf("* ");
				else printf(" ");
			}
			printf("\n");
		}*/

		int donor_index = donorIndices[0]; // First index is a random solution, as donorIndices was shuffled
		// Insert into parent
		for( int x : variables )
		{
			sample.push_back(population[donor_index]->variables[x]);
			assert(x == variables[sample.size()-1] ); // check if the order of the indices in 'variables' is the same as the 'variable_groups' when concatenated
		}
	}
	
	assert( sample.size() == variables.size() );

	return sample;
}


void factorization_t::print()
{
	printf("[");
	for( int x : variables )
		printf(" %d",x);
	printf("]->[");
	for( int x : variables_conditioned_on )
		printf(" %d",x);
	printf("]\n");
}

}
