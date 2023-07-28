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
#include "gomea/src/real_valued/linkage_model.hpp"
#include <queue>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
int		  static_linkage_tree = 0,                  /* Whether the FOS is fixed throughout optimization. */
		  random_linkage_tree = 0;                  /* Whether the fixed linkage tree is learned based on a random distance measure. */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

linkage_model_rv_pt linkage_model_rv_t::createFOSInstance( const linkage_config_t &config, size_t numberOfVariables, const graph_t &VIG ) 
{
	if( config.type != linkage::FROM_FILE )
		assert( numberOfVariables > 0 );
	switch( config.type )
	{
		case linkage::UNIVARIATE: return univariate(numberOfVariables);
		case linkage::MPM: return marginal_product_model(numberOfVariables, config.mpm_block_size);
		case linkage::FULL: return full(numberOfVariables);
		case linkage::LINKAGE_TREE: return linkage_tree(numberOfVariables, config.lt_similarity_measure, config.lt_filtered, config.lt_maximum_set_size, config.lt_is_static );
		case linkage::CONDITIONAL: return conditional(numberOfVariables,VIG,config.cond_max_clique_size,config.cond_include_cliques_as_fos_elements,config.cond_include_full_fos_element);
		case linkage::CUSTOM_LM: return custom_fos(numberOfVariables,config.FOS);
		case linkage::FROM_FILE: return from_file(config.filename);
	}
	throw std::runtime_error("Unknown linkage model.\n");
}

// Copy a FOS
linkage_model_rv_t::linkage_model_rv_t( const linkage_model_rv_t &other ) : linkage_model_t(other.number_of_variables)
{
	for(size_t i = 0; i < other.FOSStructure.size(); i++ )
	{
		std::vector<int> vec;
		for(size_t j = 0; j < FOSStructure[i].size(); j++ )
			vec.push_back(other.FOSStructure[i][j]);
		addGroup(vec);
	}
	p_accept = other.p_accept;
	is_static = other.is_static;
}

linkage_model_rv_t::linkage_model_rv_t( size_t number_of_variables, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element ) : linkage_model_t(number_of_variables)
{
	this->include_cliques_as_fos_elements = include_cliques_as_fos_elements;
	this->include_full_fos_element = include_full_fos_element;
	this->is_conditional = true;
	this->is_static = true; 
	this->max_clique_size = max_clique_size;
	assert( include_cliques_as_fos_elements || include_full_fos_element );
	
	const int UNVISITED = 0;
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	const int IN_QUEUE = 3;
	std::vector<int> visited(number_of_variables,0);
	vec_t<int> var_order = gomea::utils::randomPermutation( number_of_variables );

	assert( variable_interaction_graph.size() == number_of_variables );
	if( variable_interaction_graph.size() != number_of_variables )
		throw std::runtime_error("Incorrectly initialized Variable Interaction Graph\n");

	std::vector<int> VIG_order;
	conditional_distribution_t *full_cond = NULL;
	if( include_full_fos_element )
		full_cond = new conditional_distribution_t();
	for( int i = 0; i < number_of_variables; i++ )
	{
		int ind = var_order[i];
		if( visited[ind] == IS_VISITED )
			continue;
		visited[ind] = IN_CLIQUE;
	
		std::queue<int> q;
		q.push(ind);

		while( !q.empty() )
		{
			ind = q.front();
			q.pop();

			if( visited[ind] == IS_VISITED )
				continue;
			visited[ind] = IS_VISITED;

			VIG_order.push_back(ind);

			std::vector<int> clique;
			std::set<int> cond;
			clique.push_back(ind);
			for( int x : variable_interaction_graph.at(ind) ) // neighbors of ind
			{
				if( visited[x] == IS_VISITED )
					cond.insert(x);
			}

			for( int x : variable_interaction_graph.at(ind) ) // neighbors of ind
			{
				if( visited[x] != IS_VISITED )
				{
					bool add_to_clique = true;
					std::set<int> neighbors = variable_interaction_graph.at(x);
					if( (int) clique.size() >= max_clique_size )
						add_to_clique = false;
					if( add_to_clique )
					{
						for( int y : clique )
						{
							if( neighbors.find(y) == neighbors.end() ) // edge (x,y) does not exist
							{
								add_to_clique = false;
								//printf("no E(%d,%d)\n",x,y);
								break;
							}
						}
					}
					if( add_to_clique )
					{
						for( int y : cond )
						{
							if( neighbors.find(y) == neighbors.end() ) // edge (x,y) does not exist
							{
								add_to_clique = false;
								//printf("no E(%d,%d)\n",x,y);
								break;
							}
						}
					}
					if( add_to_clique )
						clique.push_back(x);
				}
			}
			for( int x : clique )
			{
				visited[x] = IS_VISITED;
				for( int y : variable_interaction_graph.at(x) ) // neighbors of x
				{
					if( visited[y] == UNVISITED )
					{
						q.push(y);
						visited[y] = IN_QUEUE;
					}
				}
			}
			if( include_cliques_as_fos_elements )	
				addConditionedGroup( clique, cond );
			if( include_full_fos_element )
				full_cond->addGroupOfVariables( clique, cond );
		}
	}
	if( include_full_fos_element ) 
	{
		if( size() != 1 ) // if length == 1, only 1 clique was found, which must have been the full model; in that case, do not add it again	
			addGroup( full_cond );
		else
			delete full_cond;
	}
}

linkage_model_rv_t::~linkage_model_rv_t()
{
	for( auto d : distributions )
		delete( d );
}

linkage_model_rv_pt linkage_model_rv_t::univariate(size_t numberOfVariables_)
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,1));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::full(size_t numberOfVariables_)
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,numberOfVariables_));
	return( new_fos );
}

linkage_model_rv_pt linkage_model_rv_t::marginal_product_model( size_t numberOfVariables_, size_t block_size )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,block_size));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::conditional( size_t numberOfVariables_, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,variable_interaction_graph,max_clique_size,include_cliques_as_fos_elements,include_full_fos_element));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::custom_fos( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,FOS));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::linkage_tree(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static_ ) 
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_, similarityMeasure_, filtered_, maximumSetSize_, is_static_));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::from_file( std::string filename )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(filename));
	return( new_fos );
}

void linkage_model_rv_t::clearDistributions()
{
	distributions.clear();
}

void linkage_model_rv_t::initializeDistributions()
{
	assert( distributions.size() == 0 || distributions.size() == FOSStructure.size() );
	if( distributions.size() == 0 )
	{
		for( vec_t<int> group : FOSStructure )
			distributions.push_back( new normal_distribution_t(group) );
	}
}
			
double linkage_model_rv_t::getDistributionMultiplier( int element_index )
{
	return( distributions[element_index]->distribution_multiplier );
}

double linkage_model_rv_t::getAcceptanceRate() 
{
	return( p_accept );
}

void linkage_model_rv_t::addGroup( std::vector<int> group ) 
{
	std::sort(group.begin(),group.end());
	FOSStructure.push_back(group);
	distributions.push_back( new normal_distribution_t(group) );
}

void linkage_model_rv_t::addGroup( distribution_t *dist )
{
	//std::sort(dist->variables.begin(),dist->variables.end());
	FOSStructure.push_back(dist->variables);
	distributions.push_back( dist );
}

void linkage_model_rv_t::addConditionedGroup( std::vector<int> variables ) 
{
	std::set<int> cond;
	addConditionedGroup(variables,cond);
}

void linkage_model_rv_t::addConditionedGroup( std::vector<int> variables, std::set<int> conditioned_variables )
{
	std::sort(variables.begin(),variables.end());
	FOSStructure.push_back(variables);
	conditional_distribution_t *dist = new conditional_distribution_t(variables,conditioned_variables);
	distributions.push_back(dist);
}

// Learn a linkage tree
void linkage_model_rv_t::learnLinkageTreeFOS( matE covariance_matrix )
{
	assert( type == linkage::LINKAGE_TREE );
	assert( !is_static );

	/* Compute Mutual Information matrix */
	vec_t<vec_t<double>> MI_matrix = computeMIMatrix( covariance_matrix, number_of_variables );
	linkage_model_t::learnLinkageTreeFOS(MI_matrix,true);
	clearDistributions();
	initializeDistributions();
}

void linkage_model_rv_t::learnLinkageTreeFOS( vec_t<vec_t<double>> similarity_matrix, bool include_full_fos_element )
{
	assert( type == linkage::LINKAGE_TREE );
	assert( !is_static );
	linkage_model_t::learnLinkageTreeFOS(similarity_matrix,include_full_fos_element);
	clearDistributions();
	initializeDistributions();
}

void linkage_model_rv_t::randomizeOrder( const graph_t &variable_interaction_graph ) 
{
	std::vector<int> visited(number_of_variables,0);
	std::vector<int> VIG_order = getVIGOrderBreadthFirst(variable_interaction_graph);
	/*printf("VIG_ORDER: ");
	for(int i = 0; i < VIG_order.size(); i++ )
		printf("%d ",VIG_order[i]);
	printf("\n");*/

	int FOS_length = 0;
	FOSorder = vec_t<int>(size());
	if( include_cliques_as_fos_elements )
	{
		for(int i = 0; i < number_of_variables; i++ )
		{
			assert( FOSStructure[i][0] == i );
			assert( elementSize(i) == 1 );
			FOSorder[i] = VIG_order[i];
			distributions[FOSorder[i]]->updateConditionals(variable_interaction_graph,visited);
		}
		FOS_length = number_of_variables;
	}
	for( int i = 0; i < number_of_variables; i++ )
		visited[i] = 0;
	if( include_full_fos_element )
	{
		FOSorder[FOS_length] = number_of_variables;
		distributions[FOS_length]->setOrder(VIG_order);
		distributions[FOS_length]->updateConditionals(variable_interaction_graph,visited);
		//distributions[FOS_length]->print();
	}
}

std::vector<int> linkage_model_rv_t::getVIGOrderBreadthFirst( const graph_t &variable_interaction_graph ) 
{
	const int UNVISITED = 0;
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	const int IN_QUEUE = 3;
	std::vector<int> visited(number_of_variables,0);
	vec_t<int> var_order = gomea::utils::randomPermutation( number_of_variables );

	std::vector<int> VIG_order;
	for( int i = 0; i < number_of_variables; i++ )
	{
		int ind = var_order[i];
		if( visited[ind] == IS_VISITED )
			continue;
		visited[ind] = IN_CLIQUE;
	
		std::queue<int> q;
		q.push(ind);

		while( !q.empty() )
		{
			ind = q.front();
			q.pop();

			if( visited[ind] == IS_VISITED )
				continue;
			visited[ind] = IS_VISITED;

			VIG_order.push_back(ind);

			for( int x : variable_interaction_graph.at(ind) )
			{
				if( visited[x] == UNVISITED )
				{
					q.push(x);
					visited[x] = IN_QUEUE;
					//printf("Q[ %d ]\n",x);
				}
			}
		}
	}
	return( VIG_order );
}

vec_t<vec_t<double>> linkage_model_rv_t::computeMIMatrix( matE covariance_matrix, int n )
{
    vec_t<vec_t<double>> MI_matrix;
	MI_matrix.resize(numberOfVariables);        
    for (size_t i = 0; i < numberOfVariables; ++i)
        MI_matrix[i].resize(numberOfVariables);         
	for(int i = 0; i < n; i++ )
	{
		MI_matrix[i][i] = 1e20;
		for(int j = 0; j < i; j++ )
		{
			double si = sqrt(covariance_matrix(i,i));
			double sj = sqrt(covariance_matrix(j,j));
			double r = covariance_matrix(i,j)/(si*sj);
			MI_matrix[i][j] = log(sqrt(1/(1-r*r)));
			MI_matrix[j][i] = MI_matrix[i][j];
		}
	}
	return( MI_matrix );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void linkage_model_rv_t::inheritDistributionMultipliers( linkage_model_rv_t *other, double *multipliers )
{
	int      i, *permutation;
	double   *multipliers_copy;

	multipliers_copy = (double*) utils::Malloc(size()*sizeof(double));
	for( i = 0; i < size(); i++ )
		multipliers_copy[i] = multipliers[i];

	permutation = matchFOSElements( other );

	for( i = 0; i < size(); i++ )
		multipliers[permutation[i]] = multipliers_copy[i];

	free( multipliers_copy );
	free( permutation );
}

int *linkage_model_rv_t::matchFOSElements( linkage_model_rv_t *other )
{
	int *permutation = (int *) utils::Malloc( size()*sizeof(int));
	int **FOS_element_similarity_matrix = (int**) utils::Malloc((size()-number_of_variables)*sizeof(int*));
	for(int i = 0; i < size()-number_of_variables; i++ )
		FOS_element_similarity_matrix[i] = (int*) utils::Malloc((size()-number_of_variables)*sizeof(int));
	for(int i = 0; i < number_of_variables; i++ )
	{
		for(int j = 0; j < number_of_variables; j++ )
		{
			if( other->FOSStructure[i][0] == FOSStructure[j][0] )
			{
				permutation[i] = j;
				break;
			}
		}
	}
	for(int i = number_of_variables; i < size(); i++ )
	{
		for(int j = number_of_variables; j < size(); j++ )
		{
			size_t a = 0;
			size_t b = 0;
			int matches = 0;
			while( a < other->FOSStructure[i].size() && b < FOSStructure[j].size() )
			{
				if( other->FOSStructure[i][a] < FOSStructure[j][b] )
					a++;
				else if( other->FOSStructure[i][a] > FOSStructure[j][b] )
					b++;
				else
				{
					a++;
					b++;
					matches++;
				}
			}
			FOS_element_similarity_matrix[i-number_of_variables][j-number_of_variables] = (int) 10000*(2.0*matches/(other->FOSStructure[i].size()+FOSStructure[j].size()));
		}
	}

	int *hungarian_permutation = hungarianAlgorithm(FOS_element_similarity_matrix, size()-number_of_variables);
	for(int i = 0; i < size()-number_of_variables; i++ )
		permutation[i+number_of_variables] = hungarian_permutation[i]+number_of_variables;

	for(int i = 0; i < size()-number_of_variables; i++ )
		free( FOS_element_similarity_matrix[i] );
	free( FOS_element_similarity_matrix );
	free( hungarian_permutation );

	return( permutation );
}

int *linkage_model_rv_t::hungarianAlgorithm( int **similarity_matrix, int dim )
{
	int x,y,ty;

	int *lx = (int*) utils::Malloc(dim*sizeof(int));
	int *ly = (int*) utils::Malloc(dim*sizeof(int));
	int *xy = (int*) utils::Malloc(dim*sizeof(int));
	int *yx = (int*) utils::Malloc(dim*sizeof(int));
	int *slack = (int*) utils::Malloc(dim*sizeof(int));
	int *slackx = (int*) utils::Malloc(dim*sizeof(int));
	int *prev = (int*) utils::Malloc(dim*sizeof(int));
	bool *S = (bool*) utils::Malloc(dim*sizeof(bool));
	bool *T = (bool*) utils::Malloc(dim*sizeof(bool));

	int root = -1;
	int max_match = 0;
	for(int i = 0; i < dim; i++ )
	{
		lx[i] = 0;
		ly[i] = 0;
		xy[i] = -1;
		yx[i] = -1;
	}
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			if(similarity_matrix[i][j] > lx[i])
				lx[i] = similarity_matrix[i][j];

	bool terminated = false;
	while(!terminated)
	{
		if (max_match == dim) break;

		int wr = 0;
		int rd = 0;
		int *q = (int*) utils::Malloc(dim*sizeof(int));
		for(int i = 0; i < dim; i++ )
		{
			S[i] = false;
			T[i] = false;
			prev[i] = -1;
		}

		for (x = 0; x < dim; x++)
		{
			if (xy[x] == -1)
			{
				q[wr++] = root = x;
				prev[x] = -2;
				S[x] = true;
				break;
			}
		}

		for (y = 0; y < dim; y++)
		{
			slack[y] = lx[root] + ly[y] - similarity_matrix[root][y];
			slackx[y] = root;
		}

		while ( 1 )
		{
			while (rd < wr)
			{
				x = q[rd++];
				for (y = 0; y < dim; y++)
				{
					if (similarity_matrix[x][y] == lx[x] + ly[y] && !T[y])
					{
						if (yx[y] == -1) break;
						T[y] = true;
						q[wr++] = yx[y];
						hungarianAlgorithmAddToTree(yx[y], x, S, prev, slack, slackx, lx, ly, similarity_matrix, dim);
					}
				}
				if (y < dim) break;
			}
			if (y < dim) break;

			int delta = 100000000;
			for(y = 0; y < dim; y++)
				if(!T[y] && slack[y] < delta)
					delta = slack[y];
			for(x = 0; x < dim; x++)
				if(S[x])
					lx[x] -= delta;
			for(y = 0; y < dim; y++)
				if(T[y])
					ly[y] += delta;
			for(y = 0; y < dim; y++)
				if(!T[y])
					slack[y] -= delta;

			wr = 0;
			rd = 0;
			for (y = 0; y < dim; y++)
			{
				if (!T[y] && slack[y] == 0)
				{
					if (yx[y] == -1)
					{
						x = slackx[y];
						break;
					}
					else
					{
						T[y] = true;
						if (!S[yx[y]])
						{
							q[wr++] = yx[y];
							hungarianAlgorithmAddToTree(yx[y], slackx[y], S, prev, slack, slackx, lx, ly, similarity_matrix, dim);
						}
					}
				}
			}
			if (y < dim) break;
		}

		if (y < dim)
		{
			max_match++;
			for (int cx = x, cy = y; cx != -2; cx = prev[cx], cy = ty)
			{
				ty = xy[cx];
				yx[cy] = cx;
				xy[cx] = cy;
			}
		}
		else terminated = true;

		free( q );
	}

	free( lx );
	free( ly );
	free( yx );
	free( slack );
	free( slackx );
	free( prev );
	free( S );
	free( T );

	return xy;
}

void linkage_model_rv_t::hungarianAlgorithmAddToTree(int x, int prevx, bool *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim) 
{
	S[x] = true;
	prev[x] = prevx;
	for (int y = 0; y < dim; y++)
	{
		if (lx[x] + ly[y] - similarity_matrix[x][y] < slack[y])
		{
			slack[y] = lx[x] + ly[y] - similarity_matrix[x][y];
			slackx[y] = x;
		}
	}
}

void linkage_model_rv_t::print()
{
	printf("FOS: {");
	for(int i = 0; i < size(); i++ )
	{
		printf("[");
		for(int j = 0; j < (int) FOSStructure[i].size(); j++ )
		{
			printf("%d", FOSStructure[i][j]);
			if( j != ((int)FOSStructure[i].size())-1)
				printf(",");
		}
		printf("]");
		printf(",");
	}
	printf("}\n");
}

partial_solution_t<double> *linkage_model_rv_t::generatePartialSolution( int FOS_index, solution_t<double> *solution_conditioned_on, fitness::fitness_generic_t *fitness_function )
{
	return( distributions[FOS_index]->generatePartialSolution(solution_conditioned_on, fitness_function) );
}

void linkage_model_rv_t::estimateDistributions( solution_t<double> **selection, int selection_size )
{
	for( int i = 0; i < size(); i++ )
		estimateDistribution( i, selection, selection_size );
	FOSorder = gomea::utils::randomPermutation( size() );
}

void linkage_model_rv_t::estimateDistribution( int FOS_index, solution_t<double> **selection, int selection_size )
{
	distributions[FOS_index]->estimateDistribution(selection,selection_size);
}

void linkage_model_rv_t::adaptDistributionMultiplier( int FOS_index, partial_solution_t<double> **solutions, int num_solutions )
{
	if( no_improvement_stretch >= maximum_no_improvement_stretch )
		distributions[FOS_index]->adaptDistributionMultiplierMaximumStretch(solutions,num_solutions);
	else
		distributions[FOS_index]->adaptDistributionMultiplier(solutions,num_solutions);
}

}}
