#pragma once

#include "gomea/src/utils/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/cond_factor.hpp"
#include "gomea/src/common/factorization.hpp"
#include "gomea/src/common/distribution.hpp"
#include "gomea/src/common/sampler.hpp"
#include "gomea/src/real_valued/linkage_model.hpp"
#include "gomea/src/real_valued/distribution.hpp"

namespace gomea{

class sampler_Rt : public sampler_t<double>
{
    public:
        sampler_Rt( fitness::fitness_t<double> *fitness_function ) : sampler_t(fitness_function){};

        partial_solution_t<double> *sampleSolution( linkage_model_pt linkage_model, int FOS_index, solution_t<double> *solution_conditioned_on = nullptr );
        partial_solution_t<double> *samplePartialSolution( cond_factor_t *factor, const vec_t<double> &solution_conditioned_on = vec_t<double>() );
        partial_solution_t<double> *sampleFullSolutionConditional( factorization_t *factorization, const vec_t<double> &solution_conditioned_on );
};

}