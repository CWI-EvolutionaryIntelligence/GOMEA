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

linkage_model_rv_pt linkage_model_rv_t::createFOSInstance( const linkage_config_t &config, size_t number_of_variables, const graph_t &VIG ) 
{
	if( config.type != linkage::FROM_FILE )
		assert( number_of_variables > 0 );
	switch( config.type )
	{
		case linkage::UNIVARIATE: return univariate(number_of_variables);
		case linkage::MPM: return marginal_product_model(number_of_variables, config.mpm_block_size);
		case linkage::FULL: return full(number_of_variables);
		case linkage::LINKAGE_TREE: return linkage_tree(number_of_variables, config.lt_similarity_measure, config.lt_filtered, config.lt_maximum_set_size, config.lt_is_static );
		case linkage::CONDITIONAL: return conditional(number_of_variables,VIG,config.cond_max_clique_size,config.cond_include_cliques_as_fos_elements,config.cond_include_full_fos_element);
		case linkage::CUSTOM_LM: return custom_fos(number_of_variables,config.FOS);
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

linkage_model_rv_pt linkage_model_rv_t::univariate(size_t number_of_variables_)
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_,1));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::full(size_t number_of_variables_)
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_,number_of_variables_));
	return( new_fos );
}

linkage_model_rv_pt linkage_model_rv_t::marginal_product_model( size_t number_of_variables_, size_t block_size )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_,block_size));
	return( new_fos );
}
			
linkage_model_rv_pt linkage_model_rv_t::conditional( size_t number_of_variables_, const graph_t &variable_interaction_graph, int max_clique_size, bool include_cliques_as_fos_elements, bool include_full_fos_element )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_,variable_interaction_graph,max_clique_size,include_cliques_as_fos_elements,include_full_fos_element));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::custom_fos( size_t number_of_variables_, const vec_t<vec_t<int>> &FOS )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_,FOS));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::linkage_tree(size_t number_of_variables_, int similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static_ ) 
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(number_of_variables_, similarityMeasure_, filtered_, maximumSetSize_, is_static_));
	return( new_fos );
}
    
linkage_model_rv_pt linkage_model_rv_t::from_file( std::string filename )
{
	linkage_model_rv_pt new_fos = std::shared_ptr<linkage_model_rv_t>(new linkage_model_rv_t(filename));
	return( new_fos );
}

distribution_Rt *linkage_model_rv_t::getDistribution( int FOS_index ){
	return( dynamic_cast<distribution_Rt*>(factorization->factors[FOS_index]->distribution) );
}

double linkage_model_rv_t::getDistributionMultiplier( int FOS_index )
{
	return( distribution_multipliers[FOS_index] );
}

void linkage_model_rv_t::adaptDistributionMultiplierNoStretch( int FOS_index, partial_solution_t<double>** partial_solutions, int num_solutions )
{
	bool improvementForFOSElement = false;
	for(int i = 0; i < num_solutions; i++ )
		if( partial_solutions[i]->improves_elitist )
			improvementForFOSElement = true;

	distribution_Rt *dist = getDistribution(FOS_index);
	if( (((double) dist->out_of_bounds_draws)/((double) dist->samples_drawn)) > 0.9 )
		distribution_multipliers[FOS_index] *= 0.5;

	double st_dev_ratio = dist->getStandardDeviationRatio(partial_solutions, num_solutions);

	if( improvementForFOSElement )
	{
		if( distribution_multipliers[FOS_index] < 1.0 )
			distribution_multipliers[FOS_index] = 1.0;

		if( st_dev_ratio > st_dev_ratio_threshold )
			distribution_multipliers[FOS_index] *= distribution_multiplier_increase;
	}
	else
	{
		if( distribution_multipliers[FOS_index] > 1.0 )
			distribution_multipliers[FOS_index] *= distribution_multiplier_decrease;

		if( distribution_multipliers[FOS_index] < 1.0)
			distribution_multipliers[FOS_index] = 1.0;
	}
}

void linkage_model_rv_t::adaptDistributionMultiplierMaximumStretch( int FOS_index, partial_solution_t<double>** partial_solutions, int num_solutions )
{
	bool improvementForFOSElement = false;
	for(int i = 0; i < num_solutions; i++ )
		if( partial_solutions[i]->improves_elitist )
			improvementForFOSElement = true;
			
	distribution_Rt *dist = getDistribution(FOS_index);
	if( (((double) dist->out_of_bounds_draws)/((double) dist->samples_drawn)) > 0.9 )
		distribution_multipliers[FOS_index] *= 0.5;

	double st_dev_ratio = dist->getStandardDeviationRatio(partial_solutions, num_solutions);

	if( improvementForFOSElement )
	{
		if( distribution_multipliers[FOS_index] < 1.0 )
			distribution_multipliers[FOS_index] = 1.0;

		if( st_dev_ratio > st_dev_ratio_threshold )
			distribution_multipliers[FOS_index] *= distribution_multiplier_increase;
	}
	else
	{
		distribution_multipliers[FOS_index] *= distribution_multiplier_decrease;
	}
}

double linkage_model_rv_t::getAcceptanceRate() 
{
	return( p_accept );
}

// Learn a linkage tree
void linkage_model_rv_t::learnLinkageTreeFOS( matE covariance_matrix )
{
	assert( type == linkage::LINKAGE_TREE );
	assert( !is_static );

	/* Compute Mutual Information matrix */
	vec_t<vec_t<double>> MI_matrix = computeMIMatrix( covariance_matrix, number_of_variables );
	linkage_model_t::learnLinkageTreeFOS(MI_matrix,true);
	//clearDistributions();
	//initializeDistributions();
}

void linkage_model_rv_t::learnLinkageTreeFOS( vec_t<vec_t<double>> similarity_matrix, bool include_full_fos_element )
{
	assert( type == linkage::LINKAGE_TREE );
	assert( !is_static );
	linkage_model_t::learnLinkageTreeFOS(similarity_matrix,include_full_fos_element);
	//clearDistributions();
	//initializeDistributions();
}

vec_t<vec_t<double>> linkage_model_rv_t::computeMIMatrix( matE covariance_matrix, int n )
{
    vec_t<vec_t<double>> MI_matrix;
	MI_matrix.resize(number_of_variables);        
    for (size_t i = 0; i < number_of_variables; ++i)
        MI_matrix[i].resize(number_of_variables);         
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

	int *hungarian_permutation = utils::hungarianAlgorithm(FOS_element_similarity_matrix, size()-number_of_variables);
	for(int i = 0; i < size()-number_of_variables; i++ )
		permutation[i+number_of_variables] = hungarian_permutation[i]+number_of_variables;

	for(int i = 0; i < size()-number_of_variables; i++ )
		free( FOS_element_similarity_matrix[i] );
	free( FOS_element_similarity_matrix );
	free( hungarian_permutation );

	return( permutation );
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

void linkage_model_rv_t::estimateDistributionForFOSElement( int FOS_index, solution_t<double> **selection, int selection_size )
{
	if(distribution_multipliers.size() != size())
		distribution_multipliers.resize(size());
	distribution_multipliers[FOS_index] = 1.0;

	if(FOSStructure[FOS_index].size() == number_of_variables)
	{
		for( int i = 0; i < factorization->size(); i++ )
			estimateDistributionForFactor( i, selection, selection_size, distribution_multipliers[FOS_index] );	
	}
	else
	{
		int factor_index = FOS_index;
		assert(factorization->factors[factor_index]->size() == FOSStructure[FOS_index].size());
		assert(factorization->factors[factor_index]->variables[0] == FOSStructure[FOS_index][0]);
		estimateDistributionForFactor( factor_index, selection, selection_size, distribution_multipliers[FOS_index] );
	}
}

void linkage_model_rv_t::estimateDistributionForFactor( int factor_index, solution_t<double> **selection, int selection_size, double distribution_multiplier )
{
	getDistribution(factor_index)->estimateDistribution(selection,selection_size,distribution_multiplier);
}

void linkage_model_rv_t::adaptDistributionMultiplier( int FOS_index, partial_solution_t<double> **solutions, int num_solutions )
{
	if( no_improvement_stretch >= maximum_no_improvement_stretch )
		adaptDistributionMultiplierMaximumStretch(FOS_index,solutions,num_solutions);
	else
		adaptDistributionMultiplierNoStretch(FOS_index,solutions,num_solutions);
}

bool linkage_model_rv_t::generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio )
{
	*st_dev_ratio = 0.0;
	bool generational_improvement = false;
	int number_of_improvements  = 0;
	for(int i = 0; i < num_solutions; i++ )
	{
		if( partial_solutions[i]->improves_elitist )
		{
			generational_improvement = true;
			break;
		}
	}

	for( size_t k = 0; k < factorization->size(); k++ )
	{
		double SDR = getDistribution(k)->getStandardDeviationRatio(partial_solutions,num_solutions);
		*st_dev_ratio = std::max( *st_dev_ratio, SDR );
	}
	
	return( generational_improvement );
}

}}
