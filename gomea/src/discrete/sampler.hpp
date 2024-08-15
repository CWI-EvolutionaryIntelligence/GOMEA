#pragma once

#include "gomea/src/utils/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/cond_factor.hpp"
#include "gomea/src/common/factorization.hpp"
#include "gomea/src/common/distribution.hpp"
#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/common/sampler.hpp"

namespace gomea{

class sampler_Dt : public sampler_t<char>
{
    public:
        sampler_Dt( int number_of_variables, int alphabet_size, fitness::fitness_t<char> *fitness_function ) : sampler_t(fitness_function){};

        vec_t<char> sampleSolution( linkage_model_pt linkage_model, int FOS_index, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index );
        vec_t<char> sampleFullSolutionConditional( factorization_t *factorization, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index );
        vec_t<char> samplePartialSolution( cond_factor_t *factor, solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index );
        vec_t<char> samplePartialSolution( cond_factor_t *factor, const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index );

        bool isConditionedDonor( cond_factor_t *factor, solution_t<char> *donor_candidate, solution_t<char> *parent );
        bool isConditionedDonor( cond_factor_t *factor, solution_t<char> *donor_candidate, const vec_t<char> &parent );
};

}