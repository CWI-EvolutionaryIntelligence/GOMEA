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
#include "gomea/src/real_valued/distribution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace realvalued{

distribution_t::~distribution_t() 
{
}

void distribution_t::adaptDistributionMultiplier( partial_solution_t<double>** partial_solutions, int num_solutions )
{
	bool improvementForFOSElement = false;
	if( (((double) out_of_bounds_draws)/((double) samples_drawn)) > 0.9 )
		distribution_multiplier *= 0.5;

	double st_dev_ratio;
	improvementForFOSElement = generationalImprovementForOnePopulationForFOSElement( partial_solutions, num_solutions, &st_dev_ratio );

	if( improvementForFOSElement )
	{
		if( distribution_multiplier < 1.0 )
			distribution_multiplier = 1.0;

		if( st_dev_ratio > st_dev_ratio_threshold )
			distribution_multiplier *= distribution_multiplier_increase;
	}
	else
	{
		if( distribution_multiplier > 1.0 )
			distribution_multiplier *= distribution_multiplier_decrease;

		if( distribution_multiplier < 1.0)
			distribution_multiplier = 1.0;
	}
}

void distribution_t::adaptDistributionMultiplierMaximumStretch( partial_solution_t<double>** partial_solutions, int num_solutions )
{
	bool improvementForFOSElement = false;
	if( (((double) out_of_bounds_draws)/((double) samples_drawn)) > 0.9 )
		distribution_multiplier *= 0.5;

	double st_dev_ratio;
	improvementForFOSElement = generationalImprovementForOnePopulationForFOSElement( partial_solutions, num_solutions, &st_dev_ratio );

	if( improvementForFOSElement )
	{
		if( distribution_multiplier < 1.0 )
			distribution_multiplier = 1.0;

		if( st_dev_ratio > st_dev_ratio_threshold )
			distribution_multiplier *= distribution_multiplier_increase;
	}
	else
	{
		distribution_multiplier *= distribution_multiplier_decrease;
	}
}

void distribution_t::updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited )
{}
		
void distribution_t::setOrder( const vec_t<int> &order )
{
	exit(1);
}

void distribution_t::print(){}

double distribution_t::estimateMean( int var, solution_t<double> **selection, int selection_size )
{
	double mean = 0.0;
	for(int j = 0; j < selection_size; j++ )
		mean += selection[j]->variables[var];
	mean /= (double) selection_size;
	return( mean );
}

double distribution_t::estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size )
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

vec_t<double> distribution_t::estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
{
	vec_t<double> mean_vector = vec_t<double>(variables.size());
	for( size_t i = 0; i < variables.size(); i++ )
		mean_vector[i] = estimateMean( variables[i], selection, selection_size );
	return( mean_vector );
}

matE distribution_t::estimateUnivariateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
{
	/* First do the maximum-likelihood estimate from data */
	//matE covariance_matrix = matE(variables.size(),variables.size(),fill::zeros);
	matE covariance_matrix = matE::Zero(variables.size(),variables.size());
	for( size_t j = 0; j < variables.size(); j++ )
	{
		int vara = variables[j];
		double cov = estimateCovariance(vara,vara,selection,selection_size);
		//covariance_matrix(j,k) = (1-eta_cov)*covariance_matrix(j,k)+ eta_cov*cov;
		covariance_matrix(j,j) = cov * distribution_multiplier;
	}
	return( covariance_matrix );
}

matE distribution_t::estimateRegularCovarianceMatrixML( vec_t<int> &variables, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size )
{
	matE covariance_matrix;
	covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size);
	/*int n = variables.size();
	if( selection_size < n + 1 )
		covariance_matrix = estimateUnivariateCovarianceMatrixML(variables,selection,selection_size);
	else
	{
		covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size);
		if( selection_size < (0.5*n*(n+1))+1 )
			regularizeCovarianceMatrix(covariance_matrix,mean_vector,selection,selection_size);
	}*/
	return( covariance_matrix );
}

matE distribution_t::estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
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

matE distribution_t::estimateFullCovarianceMatrixML( solution_t<double> **selection, int selection_size )
{
	/* First do the maximum-likelihood estimate from data */
	//mat covariance_matrix(variables.size(),variables.size(),fill::none);
	int num_variables = selection[0]->getNumberOfVariables();
	matE covariance_matrix = matE(num_variables,num_variables);
	for( size_t j = 0; j < num_variables; j++ )
	{
		int vara = j;
		for( size_t k = j; k < num_variables; k++ )
		{
			int varb = k;
			double cov = estimateCovariance(vara,varb,selection,selection_size);
			//covariance_matrix(j,k) = (1-eta_cov)*covariance_matrix(j,k)+ eta_cov*cov;
			covariance_matrix(j,k) = cov;
			covariance_matrix(k,j) = covariance_matrix(j,k);
		}
	}
	return( covariance_matrix );
}

bool distribution_t::regularizeCovarianceMatrix( matE &cov_mat, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size )
{
	// regularization for small populations
	double number_of_samples = (double) selection_size;
	int n = variables.size();

	// either use the univariate matrix as a prior,
	// or a diagonal matrix with the mean variance on all diagonal entries
	bool use_univariate_as_prior = true;

	double meanvar = 0.0;
	if(!use_univariate_as_prior)
	{
		for(int i = 0; i < n; ++i) {
			meanvar += cov_mat(i,i);
		}
		meanvar /= (double) n;
	}
	
	double phi = 0.0;
	// y = x.^2
	// phiMat = y'*y/t-sample.^2
	// phi = sum(sum(phiMat))
	matE squared_cov = matE(n,n);
	double temp;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			squared_cov(i,j) = 0.0;
			for(int k = 0; k < selection_size; ++k)
			{
				temp = (selection[k]->variables[variables[i]]-mean_vector[i])*(selection[k]->variables[variables[j]]-mean_vector[j]);
				squared_cov(i,j) += temp*temp;
			}
			squared_cov(i,j) /= number_of_samples;
		}
	}

	// this can be implemented faster by considering only half this matrix,
	// and we dont need to store square_cov actually.
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			phi += squared_cov(i,j) - cov_mat(i,j) * cov_mat(i,j);

	// Frobenius norm, i.e.,
	// gamma = norm(sample - prior,'fro')^2;
	double gamma = 0.0;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			if(use_univariate_as_prior) {
				temp = std::abs(cov_mat(i,j) - (( i == j ) ? cov_mat(i,i) : 0.0));
			} else {
				temp = std::abs(cov_mat(i,j) - (( i == j ) ? meanvar : 0.0));
			}
			gamma += temp*temp;
		}
	}

	double kappa = phi/gamma;
	double shrinkage = std::max(0.0,std::min(1.0,kappa/number_of_samples));
	//std::cout << "Shrinkage with factor " << shrinkage << std::endl;
	//shrinkage = std::max(1e-10,shrinkage);

	if(shrinkage == 0.0) {
		return false;
	}

	// shrinking
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			if(use_univariate_as_prior) {
				cov_mat(i,j) = (1.0 - shrinkage) * cov_mat(i,j) + ((i==j) ? shrinkage*cov_mat(i,i) : 0.0);
			} else {
				cov_mat(i,j) = (1.0 - shrinkage) * cov_mat(i,j) + ((i==j) ? shrinkage*meanvar : 0.0);
			}
		}
	}

	return true;
}

normal_distribution_t::normal_distribution_t( vec_t<int> variables )
{
	this->variables = variables;
}

void normal_distribution_t::estimateDistribution( solution_t<double> **selection, int selection_size )
{
	samples_drawn = 0;
	out_of_bounds_draws = 0;
	
	mean_vector = estimateMeanVectorML(variables,selection,selection_size);
	
	/* Change the focus of the search to the best solution */
	if( distribution_multiplier < 1.0 )
		for(size_t j = 0; j < variables.size(); j++)
			mean_vector[j] = selection[0]->variables[variables[j]];

	covariance_matrix = estimateRegularCovarianceMatrixML(variables,mean_vector,selection,selection_size);
	cholesky_decomposition = utils::choleskyDecomposition( covariance_matrix );
}

partial_solution_t<double> *normal_distribution_t::generatePartialSolution( solution_t<double> *parent, fitness::fitness_generic_t *fitness_function )
{
	vec_t<int> indices = variables; 
	int num_indices = variables.size();
	vec_t<double> result(num_indices);

	int times_not_in_bounds = -1;
	out_of_bounds_draws--;

	bool ready = false;
	do
	{
		times_not_in_bounds++;
		samples_drawn++;
		out_of_bounds_draws++;

		if( times_not_in_bounds >= 100 )
		{
			printf("Sampled out of bounds too many times.\n");
			
			for( size_t i = 0; i < result.size(); i++ )
				printf("%10.3e ",mean_vector[i]);
			printf("\n");
			for (size_t i = 0; i < result.size(); i++)
				printf("%10.3e ",fitness_function->getLowerRangeBound(i));
			printf("\n");
			for (size_t i = 0; i < result.size(); i++)
				printf("%10.3e ",fitness_function->getUpperRangeBound(i));
			printf("\n");

			exit(1);
			// TODO
			/*result = vec(num_indices,fill::none);
			sample_means = vec(num_indices,fill::none);
			sample_zs = zeros<vec>(num_indices);
			for(int i = 0; i < num_indices; i++ )
			{
				result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randu<double>();
				sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*0.5;
			}*/
		}
		else
		{
			vecE sample = cholesky_decomposition * random1DNormalUnitVector(num_indices);
			for( size_t i = 0; i < sample.size(); i++ )
			{
				result[i] = mean_vector[i] + sample[i];
			}
		}
			
		ready = true;
		if (fitness_function != NULL)
		{
			for (size_t i = 0; i < result.size(); i++)
			{
				if (!fitness_function->isParameterInRangeBounds(result[i], indices[i]))
				{
					ready = false;
					break;
				}
			}
		}
	}
	while( !ready );

	partial_solution_t<double> *res_sol = new partial_solution_t<double>(result,indices);
	res_sol->setSampleMean( mean_vector );
	return( res_sol );
}

bool normal_distribution_t::generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio )
{
	bool generational_improvement = false;
	vec_t<int> indices = partial_solutions[0]->touched_indices; 
	int num_indices = indices.size();

	std::vector<double> average_z_of_improvements(num_indices,0.0);

	int number_of_improvements  = 0;
	//matE cholinv = pinv( trimatl( cholesky_decomposition ) );
	matE cholinv = utils::pinv(cholesky_decomposition.triangularView<Eigen::Lower>());
	for(int i = 0; i < num_solutions; i++ )
	{
		if( partial_solutions[i]->improves_elitist )
		{
			number_of_improvements++;
			vecE d = vecE(num_indices);
			for( int j = 0; j < num_indices; j++ )
			{
				d[j] = (partial_solutions[i]->touched_variables[j] - partial_solutions[i]->sample_means[j]);
			}
			vecE z = cholinv * d;
			for(int j = 0; j < num_indices; j++ )
			{
				average_z_of_improvements[j] += z[j];
			}
		}
	}

	// Determine st.dev. ratio
	*st_dev_ratio = 0.0;
	if( number_of_improvements > 0 )
	{
		for(int i = 0; i < num_indices; i++ )
		{
			average_z_of_improvements[i] /= (double) number_of_improvements;
			*st_dev_ratio = std::max( *st_dev_ratio, std::abs(average_z_of_improvements[i]) );
		}

		generational_improvement = true;
	}

	return( generational_improvement );
}

conditional_distribution_t::conditional_distribution_t() 
{
}

conditional_distribution_t::conditional_distribution_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

conditional_distribution_t::conditional_distribution_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables )
{
	addGroupOfVariables(variables,conditioned_variables);
}

void conditional_distribution_t::initializeMemory()
{
	mean_vectors.resize(variable_groups.size());
	mean_vectors_conditioned_on.resize(variable_groups.size());
	covariance_matrices.resize(variable_groups.size());
	rho_matrices.resize(variable_groups.size());
	cholesky_decompositions.resize(variable_groups.size());
	samples_drawn = 0;
	out_of_bounds_draws = 0;
}

void conditional_distribution_t::addGroupOfVariables( vec_t<int> indices, vec_t<int> indices_cond )
{
	std::sort(indices.begin(),indices.end());
	std::sort(indices_cond.begin(),indices_cond.end());
	vec_t<int> indices_map;
	for( int i : indices )
	{
		indices_map.push_back(variables.size());
		variables.push_back(i);
	}
	//std::sort(variables.begin(),variables.end());
	index_in_var_array.push_back(indices_map);
	variable_groups.push_back(indices);
	variables_conditioned_on.push_back(indices_cond);
	order.push_back(order.size());
}

void conditional_distribution_t::addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond )
{
	vec_t<int> cond;
	for( int i : indices_cond )
		cond.push_back(i);
	addGroupOfVariables(indices,cond);
}
			
void conditional_distribution_t::addGroupOfVariables( int index, const vec_t<int> &indices_cond )
{
	vec_t<int> indices;
	indices.push_back(index);
	addGroupOfVariables(indices, indices_cond);
}
			
void conditional_distribution_t::addGroupOfVariables( int index, int index_cond )
{
	vec_t<int> indices, indices_cond;
	indices.push_back(index);
	indices_cond.push_back(index_cond);
	addGroupOfVariables(indices, indices_cond);
}
			
void conditional_distribution_t::setOrder( const vec_t<int> &order ) 
{
	this->order = order;
}

void conditional_distribution_t::print()
{
	for(size_t i = 0; i < variable_groups.size(); i++ )
	{
		int og = order[i];
		printf("[");
		for( int x : variable_groups[og] )
			printf(" %d",x);
		printf("]->[");
		for( int x : variables_conditioned_on[og] )
			printf(" %d",x);
		printf("],");
	}
	printf("\n");
}

}}
