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

void distribution_Rt::print(){}

double distribution_Rt::estimateMean( int var, solution_t<double> **selection, int selection_size )
{
	double mean = 0.0;
	for(int j = 0; j < selection_size; j++ )
		mean += selection[j]->variables[var];
	mean /= (double) selection_size;
	return( mean );
}

double distribution_Rt::estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size )
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

vec_t<double> distribution_Rt::estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size )
{
	vec_t<double> mean_vector = vec_t<double>(variables.size());
	for( size_t i = 0; i < variables.size(); i++ )
		mean_vector[i] = estimateMean( variables[i], selection, selection_size );
	return( mean_vector );
}

matE distribution_Rt::estimateUnivariateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size, double distribution_multiplier )
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

matE distribution_Rt::estimateRegularCovarianceMatrixML( vec_t<int> &variables, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size, double distribution_multiplier )
{
	matE covariance_matrix;
	covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size,distribution_multiplier);
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

matE distribution_Rt::estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size, double distribution_multiplier )
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

matE distribution_Rt::estimateFullCovarianceMatrixML( solution_t<double> **selection, int selection_size, double distribution_multiplier )
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
			covariance_matrix(j,k) = cov * distribution_multiplier;
			covariance_matrix(k,j) = covariance_matrix(j,k);
		}
	}
	return( covariance_matrix );
}

bool distribution_Rt::regularizeCovarianceMatrix( matE &cov_mat, vec_t<double> &mean_vector, solution_t<double> **selection, int selection_size )
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

void distribution_Rt::estimateDistribution( solution_t<double> **selection, int selection_size, double distribution_multiplier )
{
	int n = variables.size();
	int n_cond = variables_conditioned_on.size();

	samples_drawn = 0;
	out_of_bounds_draws = 0;
	
	mean_vector = estimateMeanVectorML(variables,selection,selection_size);
	
	/* Change the focus of the search to the best solution */
	if( distribution_multiplier < 1.0 )
		for(size_t j = 0; j < variables.size(); j++)
			mean_vector[j] = selection[0]->variables[variables[j]];

	covariance_matrix = estimateRegularCovarianceMatrixML(variables,mean_vector,selection,selection_size,distribution_multiplier);
	if( n_cond > 0 )
	{
		mean_vector_conditioned_on = estimateMeanVectorML(variables_conditioned_on,selection,selection_size);
		
		matE A12( n, n_cond );
		for(int j = 0; j < n; j++ )
			for(int k = 0; k < n_cond; k++ )
				A12(j,k) = estimateCovariance(variables[j],variables_conditioned_on[k],selection,selection_size) * distribution_multiplier;
		//mat A22 = estimateCovarianceMatrixML(vars_cond,selection,selection_size);
		matE A22 = estimateRegularCovarianceMatrixML(variables_conditioned_on,mean_vector_conditioned_on,selection,selection_size,distribution_multiplier);
		matE A22inv = utils::pinv(A22);
	   	if( 1 ) //pinv(A22inv,A22) )
		{
			rho_matrix = A12*A22inv;
			matE submat = A12*A22inv*A12.transpose();
			//mat zeros = A22*A22inv*A22 - A22;
			covariance_matrix -= submat;
		}
		else
		{
			printf("pseudo-inverse failed\n");
		}
	}

	cholesky_decomposition = utils::choleskyDecomposition( covariance_matrix );
}

vecE distribution_Rt::sample()
{
	samples_drawn++;
	vecE sample = cholesky_decomposition * random1DNormalUnitVector(variables.size());
	for( size_t i = 0; i < sample.size(); i++ )
		sample[i] += mean_vector[i];
	return sample;
}

vecE distribution_Rt::sample(const vec_t<double> &sample_means)
{
	samples_drawn++;
	vecE sample = cholesky_decomposition * utils::random1DNormalUnitVector(variables.size());
	for( size_t i = 0; i < sample.size(); i++ )
	{
		sample[i] += sample_means[i];
	}
	return sample;
}

vec_t<double> distribution_Rt::getConditionalSampleMeans(vec_t<double> solution_conditioned_on)
{
	int num_indices_cond = variables_conditioned_on.size();
	assert( num_indices_cond > 0 );
	assert( solution_conditioned_on.size() == num_variables );

	vecE cond = vecE(num_indices_cond);
	for(int i = 0; i < num_indices_cond; i++ )
		cond[i] = solution_conditioned_on[variables_conditioned_on[i]] - mean_vector_conditioned_on[i];
	vecE sample_mean_inc = rho_matrix*cond;
	vec_t<double> sample_means( variables.size() );
	for(int i = 0; i < variables.size(); i++ )
		sample_means[i] = mean_vector[i] + sample_mean_inc[i];

	return sample_means;
}

vec_t<double> distribution_Rt::getConditionalSampleMeans(solution_t<double> *solution_conditioned_on)
{
	int num_indices_cond = variables_conditioned_on.size();
	assert( num_indices_cond > 0 );

	vecE cond = vecE(num_indices_cond);
	for(int i = 0; i < num_indices_cond; i++ )
		cond[i] = solution_conditioned_on->variables[variables_conditioned_on[i]] - mean_vector_conditioned_on[i];
	vecE sample_mean_inc = rho_matrix*cond;
	vec_t<double> sample_means( variables.size() );
	for(int i = 0; i < variables.size(); i++ )
		sample_means[i] = mean_vector[i] + sample_mean_inc[i];

	return sample_means;
}

double distribution_Rt::getStandardDeviationRatio( partial_solution_t<double> **partial_solutions, int num_solutions )
{
	double SDR = 0.0;
	int number_of_improvements = 0;
	std::vector<double> average_z_of_improvements(variables.size(),0.0);		
	//matE cholinv = pseudoInverse( trimatl( cholesky_decompositions[k] ) );
	matE cholinv = utils::pinv(cholesky_decomposition.triangularView<Eigen::Lower>());
	vecE sample_means( variables.size() );
	for(int i = 0; i < num_solutions; i++ )
	{
		assert(variables.size() == partial_solutions[i]->touched_variables.size());
		if( partial_solutions[i]->improves_elitist )
		{
			number_of_improvements++;
			for(int j = 0; j < variables.size(); j++ )
			{
				int ind = variables[j]; //index_in_var_array[k][j];
				assert(variables[j] == partial_solutions[i]->touched_variables[j]);
				sample_means[j] = partial_solutions[i]->touched_variables[ind] - partial_solutions[i]->sample_means[ind];
			}
			vecE z = cholinv * sample_means; //(partial_solutions[i]->touched_variables - partial_solutions[i]->sample_means);
			for(int j = 0; j < variables.size(); j++ )
				average_z_of_improvements[j] += z[j];
		}
	}

	// Determine st.dev. ratio
	if( number_of_improvements > 0 )
	{
		for(int i = 0; i < variables.size(); i++ )
		{
			average_z_of_improvements[i] /= (double) number_of_improvements;
			SDR = std::max( SDR, std::abs(average_z_of_improvements[i]) );
		}
	}
	return SDR;
}

}}
