#include "gomea/src/common/cond_factor.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{

cond_factor_t::cond_factor_t() 
{
}

cond_factor_t::cond_factor_t( const vec_t<int> &variables )
{
	addGroupOfVariables(variables);
}

cond_factor_t::cond_factor_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

cond_factor_t::cond_factor_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

void cond_factor_t::addGroupOfVariables( const vec_t<int> &indices, const vec_t<int> &indices_cond )
{
	for( int i : indices )
		variables.push_back(i);
	for( int i : indices_cond )
		variables_conditioned_on.push_back(i);
	std::sort(variables.begin(),variables.end());
	std::sort(variables_conditioned_on.begin(),variables_conditioned_on.end());
}

void cond_factor_t::addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond )
{
	vec_t<int> cond;
	for( int i : indices_cond )
		cond.push_back(i);
	addGroupOfVariables(indices,cond);
}
			
void cond_factor_t::addGroupOfVariables( int index, const vec_t<int> &indices_cond )
{
	vec_t<int> indices;
	indices.push_back(index);
	addGroupOfVariables(indices, indices_cond);
}
			
void cond_factor_t::addGroupOfVariables( int index, int index_cond )
{
	vec_t<int> indices, indices_cond;
	indices.push_back(index);
	indices_cond.push_back(index_cond);
	addGroupOfVariables(indices, indices_cond);
}

void cond_factor_t::updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited ) 
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
	std::sort(variables_conditioned_on.begin(), variables_conditioned_on.end());
}

void cond_factor_t::initializeFrequencyTables( const vec_t<solution_t<char>*> &population ) 
{
	vec_t<vec_t<int>> donor_list;
	// TODO
	assert(0);
}

bool cond_factor_t::isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent )
{
	return isConditionedDonor(donor_candidate, parent->variables);
}

bool cond_factor_t::isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent )
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

vec_t<char> cond_factor_t::samplePartialSolutionConditional( solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	return samplePartialSolutionConditional(parent->variables, population, parent_index);
}

vec_t<char> cond_factor_t::samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index ) 
{
	vec_t<char> sample; // In the same order as 'variables'
	vec_t<int> donorIndices(population.size());
    std::iota(donorIndices.begin(), donorIndices.end(), 0);
	/*for( int i = 0; i < variables.size(); i++ )
		printf("%d ",variables[i]);
	printf("\n");*/

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


void cond_factor_t::print()
{
	printf("[");
	for( int x : variables )
		printf(" %d",x);
	printf("]->[");
	for( int x : variables_conditioned_on )
		printf(" %d",x);
	printf("]\n");
}

double cond_factor_Rt::estimateMean( int var, solution_t<double> **selection, int selection_size )
{
	double mean = 0.0;
	for(int j = 0; j < selection_size; j++ )
		mean += selection[j]->variables[var];
	mean /= (double) selection_size;
	return( mean );
}

double cond_factor_Rt::estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size )
{
	double cov = 0.0;
	double meana = estimateMean(vara,selection,selection_size);
	double meanb = estimateMean(varb,selection,selection_size);
	for(int m = 0; m < selection_size; m++ )
		cov += (selection[m]->variables[vara]-meana)*(selection[m]->variables[varb]-meanb);
	cov /= (double) selection_size;
	cov = std::max(0.0, cov);
	return( cov );
}

vec_t<double> cond_factor_Rt::estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
{
	vec_t<double> mean_vector = vec_t<double>(variables.size());
	for( size_t i = 0; i < variables.size(); i++ )
		mean_vector[i] = estimateMean( variables[i], selection, selection_size );
	return( mean_vector );
}

matE cond_factor_Rt::estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
{
	/* First do the maximum-likelihood estimate from data */
	//mat covariance_matrix(variables.size(),variables.size(),fill::none);
	matE covariance_matrix = matE(variables.size(),variables.size());
	for( size_t j = 0; j < variables.size(); j++ )
	{
		int vara = variables[j];
		for( size_t k = j; k < variables.size(); k++ )
		{
			int varb = variables[k];
			double cov = estimateCovariance(vara,varb,selection,selection_size);
			//covariance_matrix(j,k) = (1-eta_cov)*covariance_matrix(j,k)+ eta_cov*cov;
			covariance_matrix(j,k) = cov * distribution_multiplier;
			covariance_matrix(k,j) = covariance_matrix(j,k);
		}
	}
	return( covariance_matrix );
}

void cond_factor_Rt::estimateDistribution( solution_t<double> **selection, int selection_size )
{
	int n = variables.size();
	samples_drawn = 0;
	out_of_bounds_draws = 0;

	mean_vector = estimateMeanVectorML(variables,selection,selection_size);
	
	/* Change the focus of the search to the best solution */
	if( distribution_multiplier < 1.0 )
		for(int j = 0; j < n; j++)
			mean_vector[j] = selection[0]->variables[variables[j]];

	covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size);
		
	int n_cond = variables_conditioned_on.size(); 
	if( n_cond > 0 )
	{
		matE A12( n, n_cond );
		for(int j = 0; j < n; j++ )
			for(int k = 0; k < n_cond; k++ )
				A12(j,k) = estimateCovariance(variables[j],variables_conditioned_on[k],selection,selection_size) * distribution_multiplier;
		//matE A22 = estimateCovarianceMatrixML(vars_cond,selection,selection_size);
		matE A22 = estimateCovarianceMatrixML(variables_conditioned_on,selection,selection_size);
		matE A22inv = utils::pinv(A22);
	   	//if( pinv(A22inv,A22) )
		{
			rho_matrix = A12*A22inv;
			matE submat = A12*A22inv*A12.transpose();
			//matE zeros = A22*A22inv*A22 - A22;
			covariance_matrix -= submat;
			covariance_matrix = covariance_matrix.unaryExpr([](double x) { return std::max(0.0,x); });
		}
		//else
		{
			//printf("pseudo-inverse failed\n");
		}
	}
	cholesky_decomposition = utils::choleskyDecomposition( covariance_matrix );
}

partial_solution_t<double> *cond_factor_Rt::generatePartialSolution( solution_t<double> *solution_conditioned_on, fitness::fitness_generic_t *fitness_function )
{
	vec_t<double> result = vec_t<double>(variables.size());
	vec_t<double> means = vec_t<double>(variables.size());
	
	int times_not_in_bounds = -1;
	out_of_bounds_draws--;

	vecE sample_result;
	vecE sample_means = vecE(mean_vector.size());
	bool ready = false;
	do
	{
		times_not_in_bounds++;
		samples_drawn++;
		out_of_bounds_draws++;

		if( times_not_in_bounds >= 100 )
		{
			printf("[C] Sampled out of bounds too many times.\n");
			for( size_t i = 0; i < sample_result.size(); i++ )
				printf("%10.3e ",sample_result[i]);
			printf("\n");
			for( size_t i = 0; i < sample_result.size(); i++ )
				printf("%10.3e ",sample_means[i]);
			printf("\n");
			for (size_t i = 0; i < sample_result.size(); i++)
				printf("%10.3e ",fitness_function->getLowerRangeBound(variables[i]));
			printf("\n");
			for (size_t i = 0; i < sample_result.size(); i++)
				printf("%10.3e ",fitness_function->getUpperRangeBound(variables[i]));
			printf("\n");
			exit(1);
			/*vecE sample_result = vec(num_indices, fill::none);
			vecE sample_means = vecE(num_indices, fill::none);
			for(int i = 0; i < num_indices; i++ )
			{
				sample_result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randomRealUniform01();
				sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]]) * 0.5;
			}*/
		}
		else
		{
			for( int i = 0; i < mean_vector.size(); i++ )
				sample_means[i] = mean_vector[i];

			int num_indices_cond = variables_conditioned_on.size();
			if( num_indices_cond > 0 )
			{
				vecE cond = vecE(num_indices_cond);
				for(int i = 0; i < num_indices_cond; i++ )
					cond[i] = solution_conditioned_on->variables[variables_conditioned_on[i]] - mean_vector_conditioned_on[i];
				vecE sample_mean_inc = rho_matrix*cond;
				for(int i = 0; i < mean_vector.size(); i++ )
					sample_means[i] += sample_mean_inc[i];
			}
			vecE sample_zs = utils::random1DNormalUnitVector(variables.size());
			sample_result = sample_means + cholesky_decomposition * sample_zs;
		}

		ready = true;
		if( fitness_function != NULL )
		{
			for (size_t i = 0; i < sample_result.size(); i++)
			{
				if (!fitness_function->isParameterInRangeBounds(sample_result[i], variables[i]))
				{
					ready = false;
					break;
				}
			}
		}
	}
	while( !ready );
	
	for( int i = 0; i < variables.size(); i++ )
	{
		result[i] = sample_result[i];
		means[i] = sample_means[i];
	}

	partial_solution_t<double> *sol_res = new partial_solution_t<double>(result,variables);
	sol_res->setSampleMean(means);
	return(sol_res);
}

bool cond_factor_Rt::generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio )
{
	*st_dev_ratio = 0.0;
	bool generational_improvement = false;

	for( size_t k = 0; k < variable_groups.size(); k++ )
	{	
		vec_t<int> indices = variable_groups[k];
		int num_indices = variable_groups[k].size();
		int number_of_improvements  = 0;

		std::vector<double> average_z_of_improvements(num_indices,0.0);
		
		//matE cholinv = pseudoInverse( trimatl( cholesky_decompositions[k] ) );
		matE cholinv = utils::pinv(cholesky_decompositions[k].triangularView<Eigen::Lower>());
		vecE sample_means( num_indices );
		for(int i = 0; i < num_solutions; i++ )
		{
			if( partial_solutions[i]->improves_elitist )
			{
				number_of_improvements++;
				for(int j = 0; j < num_indices; j++ )
				{
					int ind = index_in_var_array[k][j];
					sample_means[j] = partial_solutions[i]->touched_variables[ind] - partial_solutions[i]->sample_means[ind];
				}
				vecE z = cholinv * sample_means; //(partial_solutions[i]->touched_variables - partial_solutions[i]->sample_means);
				for(int j = 0; j < num_indices; j++ )
					average_z_of_improvements[j] += z[j];
			}
		}
	
		// Determine st.dev. ratio
		if( number_of_improvements > 0 )
		{
			for(int i = 0; i < num_indices; i++ )
			{
				average_z_of_improvements[i] /= (double) number_of_improvements;
				*st_dev_ratio = std::max( *st_dev_ratio, std::abs(average_z_of_improvements[i]) );
			}
			generational_improvement = true;
		}
	}
	
	return( generational_improvement );
}

}
