#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/utils/linalg.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{

class cond_factor_t{
	public:
		vec_t<int> variables;
		vec_t<int> variables_conditioned_on;

		void addGroupOfVariables( const vec_t<int> &indices, const vec_t<int> &indices_cond = vec_t<int>() );
		void addGroupOfVariables( const vec_t<int> &indices, const std::set<int> &indices_cond );
		void addGroupOfVariables( int index, const vec_t<int> &indices_cond );
		void addGroupOfVariables( int index, int index_cond );
		
		void updateConditionals( const std::map<int,std::set<int>> &variable_interaction_graph, std::vector<int> &visited );
		void setOrder( const vec_t<int> &order ); 
		void print();

	protected:
		virtual ~cond_factor_t(){};
		cond_factor_t();
		cond_factor_t( const vec_t<int> &variables ); 
		cond_factor_t( const vec_t<int> &variables, const vec_t<int> &conditioned_variables );
		cond_factor_t( const vec_t<int> &variables, const std::set<int> &conditioned_variables );
};

class cond_factor_Dt : public cond_factor_t{
	public:
		cond_factor_Dt(){};
		cond_factor_Dt( const vec_t<int> &variables ) : cond_factor_t(variables){};
		cond_factor_Dt( const vec_t<int> &variables, const vec_t<int> &conditioned_variables ) : cond_factor_t(variables,conditioned_variables){};
		cond_factor_Dt( const vec_t<int> &variables, const std::set<int> &conditioned_variables ) : cond_factor_t(variables,conditioned_variables){};

		vec_t<char> samplePartialSolutionConditional( solution_t<char> *parent, const vec_t<solution_t<char>*> &population, int parent_index = -1 );
		vec_t<char> samplePartialSolutionConditional( const vec_t<char> &parent, const vec_t<solution_t<char>*> &population, int parent_index = -1 );
		bool isConditionedDonor( solution_t<char> *donor_candidate, solution_t<char> *parent );
		bool isConditionedDonor( solution_t<char> *donor_candidate, const vec_t<char> &parent );

		void initializeFrequencyTables( const vec_t<solution_t<char>*> &population );

	private:
		vec_t<vec_t<int>> frequency_tables;
};

class cond_factor_Rt : public cond_factor_t{
	public:
		cond_factor_Rt(){};
		cond_factor_Rt( const vec_t<int> &variables ) : cond_factor_t(variables){};
		cond_factor_Rt( const vec_t<int> &variables, const vec_t<int> &conditioned_variables ) : cond_factor_t(variables,conditioned_variables){};
		cond_factor_Rt( const vec_t<int> &variables, const std::set<int> &conditioned_variables ) : cond_factor_t(variables,conditioned_variables){};

		double estimateMean( int var, solution_t<double> **selection, int selection_size );
		double estimateCovariance( int vara, int varb, solution_t<double> **selection, int selection_size );
		vec_t<double> estimateMeanVectorML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		matE estimateCovarianceMatrixML( vec_t<int> &variables, solution_t<double> **selection, int selection_size );
		void estimateDistribution( solution_t<double> **selection, int selection_size, double distribution_multiplier );

		partial_solution_t<double> *generatePartialSolution( solution_t<double> *solution_conditioned_on, fitness::fitness_generic_t *fitness_function );

		bool generationalImprovementForOnePopulationForFOSElement( partial_solution_t<double>** partial_solutions, int num_solutions, double *st_dev_ratio );

	private:
		int samples_drawn = 0;
		int out_of_bounds_draws = 0;

		vec_t<double> mean_vector;
		vec_t<double> mean_vector_conditioned_on;
		matE covariance_matrix;
		matE rho_matrix;
		matE cholesky_decomposition;
};

}
