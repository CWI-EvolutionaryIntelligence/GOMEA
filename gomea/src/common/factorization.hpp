#pragma once

#include "gomea/src/utils/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/cond_factor.hpp"

namespace gomea{

class factorization_t
{
    public:
        factorization_t( int number_of_variables ) : number_of_variables(number_of_variables){};

        size_t size();

        void addGroup( int var_index );
        void addGroup( const std::set<int> &group );
        void addGroup( vec_t<int> variables, std::set<int> conditioned_variables = std::set<int>() );

        void initializeCondFactorInteractionGraph( const graph_t &variable_interaction_graph );
        
        int number_of_variables;
        vec_t<cond_factor_t*> factors;
        vec_t<int> order;
        graph_t factorization_interaction_graph;

    protected:
        //vec_t<T> samplePartialSolutionConditional( solution_t<T> *parent, const vec_t<solution_t<T>*> &population, int parent_index );
        //virtual vec_t<T> samplePartialSolutionConditional( const vec_t<T> &parent, const vec_t<solution_t<T>*> &population, int parent_index ) = 0;

        //vec_t<vec_t<int>> variable_groups;
        //vec_t<vec_t<int>> variable_groups_conditioned_on;
};

/*class sampler_Dt : public factorization_t<char>
{
    public:
        sampler_Dt( int number_of_variables, int alphabet_size ) : factorization_t(number_of_variables){};

        void initializeFactorization() override;

        vec_t<char> samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index ) override;

        bool isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent );
        bool isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent );

        vec_t<cond_factor_Dt*> factorization;
};

class sampler_Rt : public factorization_t<double>
{
    public:
        sampler_Rt( int number_of_variables ) : factorization_t(number_of_variables){};

        void initializeFactorization() override;

        vec_t<double> samplePartialSolutionConditional( const vec_t<double> &parent, const vec_t<solution_t<double>*> &population, int parent_index ) override;

        vec_t<cond_factor_Rt*> factorization;
};*/

}