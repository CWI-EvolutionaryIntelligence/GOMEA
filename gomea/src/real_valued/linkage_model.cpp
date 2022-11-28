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
int       FOS_element_ub = 0,                       /* Cut-off value for bounded fixed linkage tree (BFLT). */
		  use_univariate_FOS = 0,                   /* Whether a univariate FOS is used. */
		  learn_linkage_tree = 0,                   /* Whether the FOS is learned at the start of each generation. */
		  static_linkage_tree = 0,                  /* Whether the FOS is fixed throughout optimization. */
		  random_linkage_tree = 0;                  /* Whether the fixed linkage tree is learned based on a random distance measure. */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

linkage_model_rv_pt linkage_model_rv_t::createFOSInstance( const linkage_config_t &config, size_t numberOfVariables )
{
	if( config.type != linkage::FROM_FILE )
		assert( numberOfVariables > 0 );
	switch( config.type )
	{
		case linkage::UNIVARIATE: return univariate(numberOfVariables);
		case linkage::MPM: return marginal_product_model(numberOfVariables, config.mpm_block_size);
		case linkage::LINKAGE_TREE: return linkage_tree(numberOfVariables, config.lt_similarity_measure, config.lt_filtered, config.lt_maximum_set_size );
		case linkage::CUSTOM_LM: return custom_fos(numberOfVariables,config.FOS);
		case linkage::FROM_FILE: return from_file(config.filename);
	}
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
}

// Learn a linkage tree
linkage_model_rv_t::linkage_model_rv_t( int problem_index, double **covariance_matrix, int n ) : linkage_model_t(n)
{
	/* Compute Mutual Information matrix */
	double **MI_matrix = NULL;
	if( learn_linkage_tree )
		MI_matrix = computeMIMatrix( covariance_matrix, number_of_variables );

	/* Initialize MPM to the univariate factorization */
	int **mpm             = (int **) Malloc( number_of_variables*sizeof( int * ) );
	int *mpm_num_ind 			= (int *) Malloc( number_of_variables*sizeof( int ) );
	int mpm_length        = number_of_variables;
	int **mpm_new         = NULL;
	for(int i = 0; i < number_of_variables; i++ )
	{
		int *indices             = (int *) Malloc( 1*sizeof( int ) );
		indices[0]               = i;
		mpm[i]                   = indices;
		mpm_num_ind[i] = 1;
	}

	/* Initialize LT to the initial MPM */
	if( problem_index != 14 )
	{
		for(int i = 0; i < mpm_length; i++ )
		{
			std::vector<int> vec;
			vec.push_back(mpm[i][0]);
			addGroup(vec);
		}
	}

	/* Initialize similarity matrix */
	S_matrix = NULL;
	if( !random_linkage_tree ){
		S_matrix = (double **) Malloc( number_of_variables*sizeof( double * ) );
		for(int i = 0; i < number_of_variables; i++ )
			S_matrix[i] = (double *) Malloc( number_of_variables*sizeof( double ) );
	}

	if( learn_linkage_tree )
	{
		for(int i = 0; i < mpm_length; i++ )
			for(int j = 0; j < mpm_length; j++ )
				S_matrix[i][j] = MI_matrix[mpm[i][0]][mpm[j][0]];
		for(int i = 0; i < mpm_length; i++ )
			S_matrix[i][i] = 0;

		for(int i = 0; i < number_of_variables; i++ )
			free( MI_matrix[i] );
		free( MI_matrix );
	}
	else if( random_linkage_tree )
	{
		S_vector = (double *) Malloc( number_of_variables*sizeof(double));
		for(int i = 0; i < number_of_variables; i++ )
			S_vector[i] = randu<double>();
	}
	else if( static_linkage_tree )
	{
		if( problem_index == 0 )
		{
			random_linkage_tree = 1;
			S_vector = (double *) Malloc( number_of_variables*sizeof(double));
			for(int i = 0; i < number_of_variables; i++ )
				S_vector[i] = randu<double>();
		}
		else if( problem_index == 7 )
		{
			S_matrix[0][0] = 0.0;
			for(int i = 1; i < number_of_variables; i++ )
			{
				S_matrix[i][i] = 0.0;
				S_matrix[i-1][i] = 1e8 + randu<double>();
				S_matrix[i][i-1] = S_matrix[i-1][i];
				for(int j = i+1; j < number_of_variables; j++ )
				{
					S_matrix[j][i] = randu<double>();
					S_matrix[i][j] = S_matrix[j][i];
				}
			}
		}
		else if( problem_index == 105 || problem_index == 106 )
		{
			for(int i = 0; i < number_of_variables-1; i++ )
			{
				for(int j = i+1; j < number_of_variables; j++ )
				{
					S_matrix[i][j] = 1.0 / covariance_matrix[mpm[i][0]][mpm[j][0]];
					S_matrix[j][i] = S_matrix[i][j];
				}
				S_matrix[i][i] = 0.0;
			}
		}
		else if( problem_index == 14 )
		{
			for(int i = 0; i < number_of_variables-2; i+=2 )
			{
				S_matrix[i][i] = 0.0;
				S_matrix[i+1][i+1] = 0.0;
				S_matrix[i][i+1] = 100*number_of_variables;
				S_matrix[i+1][i] = 100*number_of_variables;
				for(int j = i+2; j < number_of_variables; j++ )
				{
					S_matrix[i][j] = 0.0;
					S_matrix[i+1][j] = 0.0;
					S_matrix[j][i] = S_matrix[i][j];
					S_matrix[j][i+1] = S_matrix[i+1][j];
				}
			}
		}
		else if( problem_index == 13 || problem_index > 1000 )
		{
			int id = problem_index;
			double rotation_angle = (id%10)*5;
			id/=10;
			double conditioning_number = id%10;
			id/=10;
			int overlap_size = id%10;
			id/=10;
			int block_size = id;
			if( problem_index == 13 )
			{
				block_size = 5;
				overlap_size = 0;
				conditioning_number = 6;
				rotation_angle = 45;
			}
			for(int i = 0; i+block_size <= number_of_variables; i+=(block_size-overlap_size) )
			{
				for( int j = 0; j < block_size; j++ )
				{
					for( int k = 0; k < j; k++ )
					{
						S_matrix[i+j][i+k] = 1e8 + randu<double>();
						S_matrix[i+k][i+j] = S_matrix[i+j][i+k];
					}
					for( int k = j; i+k < number_of_variables; k++ )
					{
						S_matrix[i+j][i+k] = randu<double>();
						S_matrix[i+k][i+j] = S_matrix[i+j][i+k];
					}
					S_matrix[i+j][i+j] = 0.0;
				}
			}
		}
		else
		{
			printf("Implement this.\n");
			exit( 0 );
		}
	}

	int *NN_chain       = (int *) Malloc( (number_of_variables+2)*sizeof( int ) );
	int NN_chain_length = 0;
	short done          = 0;
	while( !done )
	{
		if( NN_chain_length == 0 )
		{
			NN_chain[NN_chain_length] = utils::randomInt( mpm_length );
			NN_chain_length++;
		}

		if( NN_chain[NN_chain_length-1] >= mpm_length ) NN_chain[NN_chain_length-1] = mpm_length-1;

		while( NN_chain_length < 3 )
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_num_ind, mpm_length );
			NN_chain_length++;
		}

		while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_num_ind, mpm_length );
			if( ((getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length],mpm_num_ind) == getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length-2],mpm_num_ind)))
					&& (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
				NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
			NN_chain_length++;
			if( NN_chain_length > number_of_variables )
				break;
		}
		int r0 = NN_chain[NN_chain_length-2];
		int r1 = NN_chain[NN_chain_length-1];

		if( r1 >= mpm_length || r0 >= mpm_length || mpm_num_ind[r0]+mpm_num_ind[r1] > FOS_element_ub )
		{
			NN_chain_length = 1;
			NN_chain[0] = 0;
			if( FOS_element_ub < number_of_variables )
			{
				done = 1;
				for(int i = 1; i < mpm_length; i++ )
				{
					if( mpm_num_ind[i] + mpm_num_ind[NN_chain[0]] <= FOS_element_ub ) done = 0;
					if( mpm_num_ind[i] < mpm_num_ind[NN_chain[0]] ) NN_chain[0] = i;
				}
				if( done ) break;
			}
			continue;
		}

		if( r0 > r1 )
		{
			int rswap = r0;
			r0    = r1;
			r1    = rswap;
		}
		NN_chain_length -= 3;

		if( r1 < mpm_length && r1 != r0 ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
		{
			int *indices = (int *) Malloc( (mpm_num_ind[r0]+mpm_num_ind[r1])*sizeof( int ) );

			int k = 0;
			for(int j = 0; j < mpm_num_ind[r0]; j++ )
			{
				indices[k] = mpm[r0][j];
				k++;
			}
			for(int j = 0; j < mpm_num_ind[r1]; j++ )
			{
				indices[k] = mpm[r1][j];
				k++;
			}

			std::vector<int> vec;
			int *sorted = mergeSortInt(indices, mpm_num_ind[r0]+mpm_num_ind[r1]);
			for(int j = 0; j < mpm_num_ind[r0]+mpm_num_ind[r1]; j++ )
				vec.push_back(indices[sorted[j]]);
			addGroup(vec);

			free( sorted );
			free( indices );

			double mul0 = ((double) mpm_num_ind[r0])/((double) mpm_num_ind[r0]+mpm_num_ind[r1]);
			double mul1 = ((double) mpm_num_ind[r1])/((double) mpm_num_ind[r0]+mpm_num_ind[r1]);
			if( random_linkage_tree )
			{
				S_vector[r0] = mul0*S_vector[r0]+mul1*S_vector[r1];
			}
			else
			{
				for(int i = 0; i < mpm_length; i++ )
				{
					if( (i != r0) && (i != r1) )
					{
						S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
						S_matrix[r0][i] = S_matrix[i][r0];
					}
				}
			}

			int **mpm_new                   = (int **) Malloc( (mpm_length-1)*sizeof( int * ) );
			int *mpm_new_number_of_indices = (int *) Malloc( (mpm_length-1)*sizeof( int ) );
			int mpm_new_length            = mpm_length-1;
			for(int i = 0; i < mpm_new_length; i++ )
			{
				mpm_new[i]                   = mpm[i];
				mpm_new_number_of_indices[i] = mpm_num_ind[i];
			}

			mpm_new[r0] = (int*) Malloc( vec.size() * sizeof(int) );
			for(size_t i = 0; i < vec.size(); i++ )
				mpm_new[r0][i] = vec[i];
			mpm_new_number_of_indices[r0] = mpm_num_ind[r0]+mpm_num_ind[r1];
			if( r1 < mpm_length-1 )
			{
				mpm_new[r1]                   = mpm[mpm_length-1];
				mpm_new_number_of_indices[r1] = mpm_num_ind[mpm_length-1];

				if( random_linkage_tree )
				{
					S_vector[r1] = S_vector[mpm_length-1];
				}
				else
				{
					for(int i = 0; i < r1; i++ )
					{
						S_matrix[i][r1] = S_matrix[i][mpm_length-1];
						S_matrix[r1][i] = S_matrix[i][r1];
					}

					for(int j = r1+1; j < mpm_new_length; j++ )
					{
						S_matrix[r1][j] = S_matrix[j][mpm_length-1];
						S_matrix[j][r1] = S_matrix[r1][j];
					}
				}
			}

			for(int i = 0; i < NN_chain_length; i++ )
			{
				if( NN_chain[i] == mpm_length-1 )
				{
					NN_chain[i] = r1;
					break;
				}
			}

			free( mpm[r0] );
			free( mpm );
			free( mpm_num_ind );
			mpm         = mpm_new;
			mpm_num_ind = mpm_new_number_of_indices;
			mpm_length  = mpm_new_length;

			if( mpm_length == 1 )
				done = 1;
		}
	}
	free( NN_chain );

	free( mpm_new );
	free( mpm_num_ind );

	if( random_linkage_tree )
		free( S_vector );
	else
	{
		for(int i = 0; i < number_of_variables; i++ )
			free( S_matrix[i] );
		free( S_matrix );
	}
	shuffleFOS();	
}

linkage_model_rv_t::linkage_model_rv_t( size_t number_of_variables, const std::map<int,std::set<int>> &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element ) : linkage_model_t(number_of_variables)
{
	this->include_cliques_as_fos_elements = include_cliques_as_fos_elements;
	this->include_full_fos_element = include_full_fos_element;
	this->is_conditional = true;
	this->max_clique_size = max_clique_size;
	assert( include_cliques_as_fos_elements || include_full_fos_element );
	
	const int UNVISITED = 0;
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	const int IN_QUEUE = 3;
	int visited[number_of_variables]{};
	vec_t<int> var_order = gomea::utils::randomPermutation( number_of_variables );

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
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,1));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::full(size_t numberOfVariables_)
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,numberOfVariables_));
	return( new_fos );
}

linkage_model_rv_pt linkage_model_rv_t::marginal_product_model( size_t numberOfVariables_, size_t block_size )
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,block_size));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::conditional( size_t numberOfVariables_, const std::map<int,std::set<int>> &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element )
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,variable_interaction_graph,max_clique_size,include_cliques_as_fos_elements,include_full_fos_element));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::custom_fos( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS )
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_,FOS));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::linkage_tree(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_ )
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(numberOfVariables_, similarityMeasure_, filtered_, maximumSetSize_));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::from_file( std::string filename )
{
	linkage_model_rv_pt new_fos = shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(filename));
	return( new_fos );
}

void linkage_model_rv_t::initializeDistributions()
{
	distributions.clear();
	for( vec_t<int> group : FOSStructure )
		distributions.push_back( new normal_distribution_t(group) );
}
			
int linkage_model_rv_t::getDistributionMultiplier( int element_index )
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

double linkage_model_rv_t::getSimilarity( int a, int b, int *mpm_num_ind )
{
	if( FOS_element_ub < number_of_variables && mpm_num_ind[a] + mpm_num_ind[b] > FOS_element_ub ) return( 0 );
	if( random_linkage_tree ) return( 1.0-fabs(S_vector[a]-S_vector[b]) );
	return( S_matrix[a][b] );
}

int linkage_model_rv_t::determineNearestNeighbour( int index, double **S_matrix, int *mpm_num_ind, int mpm_length )
{
	int result = 0;
	if( result == index )
		result++;
	for(int i = 1; i < mpm_length; i++ )
	{
		if( ((getSimilarity(index,i,mpm_num_ind) > getSimilarity(index,result,mpm_num_ind)) || ((getSimilarity(index,i,mpm_num_ind) == getSimilarity(index,result,mpm_num_ind)) && (mpm_num_ind[i] < mpm_num_ind[result]))) && (i != index) )
			result = i;
	}

	return( result );
}

void linkage_model_rv_t::randomizeOrder( const std::map<int,std::set<int>> &variable_interaction_graph ) 
{
	int visited[number_of_variables]{};
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

std::vector<int> linkage_model_rv_t::getVIGOrderBreadthFirst( const std::map<int,std::set<int>> &variable_interaction_graph ) 
{
	const int UNVISITED = 0;
	const int IS_VISITED = 1;
	const int IN_CLIQUE = 2;
	const int IN_QUEUE = 3;
	int visited[number_of_variables]{};
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

double **linkage_model_rv_t::computeMIMatrix( double **covariance_matrix, int n )
{
	double **MI_matrix = (double **) Malloc( n*sizeof( double * ) );
	for(int j = 0; j < n; j++ )
		MI_matrix[j] = (double *) Malloc( n*sizeof( double ) );
	for(int i = 0; i < n; i++ )
	{
		MI_matrix[i][i] = 1e20;
		for(int j = 0; j < i; j++ )
		{
			double si = sqrt(covariance_matrix[i][i]);
			double sj = sqrt(covariance_matrix[j][j]);
			double r = covariance_matrix[i][j]/(si*sj);
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

	multipliers_copy = (double*) Malloc(size()*sizeof(double));
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
	int *permutation = (int *) Malloc( size()*sizeof(int));
	int **FOS_element_similarity_matrix = (int**) Malloc((size()-number_of_variables)*sizeof(int*));
	for(int i = 0; i < size()-number_of_variables; i++ )
		FOS_element_similarity_matrix[i] = (int*) Malloc((size()-number_of_variables)*sizeof(int));
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

	int *lx = (int*) Malloc(dim*sizeof(int));
	int *ly = (int*) Malloc(dim*sizeof(int));
	int *xy = (int*) Malloc(dim*sizeof(int));
	int *yx = (int*) Malloc(dim*sizeof(int));
	int *slack = (int*) Malloc(dim*sizeof(int));
	int *slackx = (int*) Malloc(dim*sizeof(int));
	int *prev = (int*) Malloc(dim*sizeof(int));
	short *S = (short*) Malloc(dim*sizeof(short));
	short *T = (short*) Malloc(dim*sizeof(short));

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

	short terminated = 0;
	while(!terminated)
	{
		if (max_match == dim) break;

		int wr = 0;
		int rd = 0;
		int *q = (int*) Malloc(dim*sizeof(int));
		for(int i = 0; i < dim; i++ )
		{
			S[i] = 0;
			T[i] = 0;
			prev[i] = -1;
		}

		for (x = 0; x < dim; x++)
		{
			if (xy[x] == -1)
			{
				q[wr++] = root = x;
				prev[x] = -2;
				S[x] = 1;
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
						T[y] = 1;
						q[wr++] = yx[y];
						hungarianAlgorithmAddToTree(yx[y], x, S, prev, slack, slackx, lx, ly, similarity_matrix, dim);
					}
				}
				if (y < dim) break;
			}
			if (y < dim) break;

			int delta = 100000000;
			for(y = 0; y < dim; y++)
				if(T[y] == 0 && slack[y] < delta)
					delta = slack[y];
			for(x = 0; x < dim; x++)
				if(S[x] == 1)
					lx[x] -= delta;
			for(y = 0; y < dim; y++)
				if(T[y] == 1)
					ly[y] += delta;
			for(y = 0; y < dim; y++)
				if(T[y] == 0)
					slack[y] -= delta;

			wr = 0;
			rd = 0;
			for (y = 0; y < dim; y++)
			{
				if (T[y] == 0 && slack[y] == 0)
				{
					if (yx[y] == -1)
					{
						x = slackx[y];
						break;
					}
					else
					{
						T[y] = 1;
						if (S[yx[y]] == 0)
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
		else terminated = 1;

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

void linkage_model_rv_t::hungarianAlgorithmAddToTree(int x, int prevx, short *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim) 
{
	S[x] = 1;
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

partial_solution_t<double> *linkage_model_rv_t::generatePartialSolution( int FOS_index, solution_t<double> *solution_conditioned_on )
{
	return( distributions[FOS_index]->generatePartialSolution(solution_conditioned_on) );
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
