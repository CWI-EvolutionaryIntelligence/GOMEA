#pragma once

#include "gomea/src/utils/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/cond_factor.hpp"
#include "gomea/src/common/factorization.hpp"
#include "gomea/src/common/distribution.hpp"
#include "gomea/src/common/linkage_model.hpp"

namespace gomea{

template<class T>
class sampler_t
{
    public:
        sampler_t<T>( fitness::fitness_t<T> *fitness_function ) : fitness_function(fitness_function){};

    protected:
        fitness::fitness_t<T> *fitness_function;
};

}