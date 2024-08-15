#include "gomea/src/real_valued/sampler.hpp"

namespace gomea{
namespace realvalued{

partial_solution_t<double> *sampler_Rt::sampleSolution( linkage_model_pt linkage_model, int FOS_index, solution_t<double> *solution_conditioned_on )
{
	if( linkage_model->FOSStructure[FOS_index].size() == fitness_function->number_of_variables && linkage_model->is_conditional)
	{
		assert(solution_conditioned_on != nullptr);
		return sampleFullSolutionConditional(linkage_model->factorization, solution_conditioned_on->variables);
	}
	else
	{
		assert( FOS_index < linkage_model->factorization->size() );
		assert( linkage_model->factorization->factors[FOS_index].size() == linkage_model->elementSize(FOS_index) );
		for(int i = 0; i < linkage_model->elementSize(FOS_index); i++ )
			assert( linkage_model->factorization->factors[FOS_index].variables[i] == linkage_model->FOSStructure[FOS_index][i] );
		return samplePartialSolution(linkage_model->factorization->factors[FOS_index]);
	}
}

partial_solution_t<double> *sampler_Rt::sampleFullSolutionConditional( factorization_t *factorization, const vec_t<double> &solution_conditioned_on ) 
{
	// Offspring sample starts as copy of parent
	int number_of_variables = factorization->number_of_variables;
	vec_t<int> variables(number_of_variables);
	std::iota(variables.begin(), variables.end(), 0);
	partial_solution_t<double> *result = new partial_solution_t<double>(solution_conditioned_on, variables);

	for( int i = 0; i < factorization->size(); i++ )
	{
		cond_factor_t *factor = factorization->factors[factorization->order[i]];
		partial_solution_t<double> *fact_sample = samplePartialSolution(factor,result->touched_variables);
		// Insert partial sample into offspring sample
		for( int j = 0; j < factor->variables.size(); j++ )
		{
			int var_index = factor->variables[j];
			assert( result->touched_indices[j] == var_index );
			result[var_index] = fact_sample->touched_variables[j];
			result->sample_means[var_index] = fact_sample->sample_means[j];
		}
		delete fact_sample;
	}
	return result;
}

partial_solution_t<double> *sampler_Rt::samplePartialSolution( cond_factor_t *factor, const vec_t<double> &solution_conditioned_on )
{
	distribution_Rt *dist = static_cast<distribution_Rt*>(factor->distribution);

	vec_t<double> result = vec_t<double>(factor->variables.size());
	
	int times_not_in_bounds = 0;

	vecE sample_result;
	vec_t<double> sample_means;
	bool ready = false;
	do
	{
		if( times_not_in_bounds >= 100 )
		{
			exit(1);
			//vecE sample_result = vec(num_indices, fill::none);
			//vecE sample_means = vecE(num_indices, fill::none);
			//for(int i = 0; i < num_indices; i++ )
			//{
			//	sample_result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randomRealUniform01();
			//	sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]]) * 0.5;
			//}
		}
		else
		{
			int num_indices_cond = factor->variables_conditioned_on.size();
			if( num_indices_cond == 0 )
			{
				sample_means = dist->mean_vector;
				sample_result = dist->sample(sample_means);
			}
			else
			{
				assert(solution_conditioned_on.size() == factor->num_variables);
				sample_means = dist->getConditionalSampleMeans(solution_conditioned_on);
				sample_result = dist->sample(sample_means);
			}
		}

		ready = true;
		if( fitness_function != NULL )
		{
			for (size_t i = 0; i < sample_result.size(); i++)
			{
				if (!fitness_function->isParameterInRangeBounds(sample_result[i], factor->variables[i]))
				{
					ready = false;
					dist->out_of_bounds_draws++;
					times_not_in_bounds++;
					break;
				}
			}
		}
	}
	while( !ready );
	
	for( int i = 0; i < factor->variables.size(); i++ )
		result[i] = sample_result[i];

	partial_solution_t<double> *sol_res = new partial_solution_t<double>(result,factor->variables);
	sol_res->setSampleMean(sample_means);
	return(sol_res);
}

}}