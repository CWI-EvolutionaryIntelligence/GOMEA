#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace fitness{

// Specifies whether a function is subject to minimization or maximization
typedef enum{
	MIN,
	MAX
} opt_mode;

template<class T>
class fitness_t
{
	protected:
		fitness_t( int number_of_variables ); 
		fitness_t( int number_of_variables, double vtr );

	public:
		virtual ~fitness_t();
		
		// Properties
		std::string name;
		int number_of_variables;
		int number_of_objectives = 1;
		opt_mode optimization_mode;

		// Termination conditions
		double  maximum_number_of_evaluations = -1.0,
				maximum_number_of_seconds = -1.0;

		// Gray-box specific
		std::map<int,std::set<int>> variable_interaction_graph; // VIG[i] lists all variables dependent on variable x_i, excluding itself
		std::map<int,std::set<int>> subfunction_dependency_map; // map[i] lists all subfunctions dependent on variable x_i
		double rotation_angle = 0.0;
		int rotation_block_size = 0;

		// Options
		double vtr; // value-to-reach
		bool black_box_optimization = false;
		bool use_vtr;
		bool vtr_hit_status = false;
		
		// Optimization progress
		double number_of_evaluations = 0.0; // discounted in GBO
		int full_number_of_evaluations = 0; // not discounted in GBO
		double elitist_objective_value = 1e308;
		double elitist_constraint_value = 1e308;

		virtual int getNumberOfSubfunctions();
		virtual int getNumberOfFitnessBuffers();

		void evaluate( solution_t<T> *solution );
		void evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution );
		void evaluatePartialSolutionBlackBox( solution_t<T> *parent, partial_solution_t<T> *solution );
		
		static fitness_t *getFitnessClass( int problem_index, int number_of_variables, double vtr );
		bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y );
		bool betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );

		void initialize();	
		void initializeSubfunctionDependencyMap();
		
		virtual vec_t<int> inputsToSubfunction( int subfunction_index );
		virtual int getIndexOfFitnessBuffer( int subfunction_index );

		bool hasVariableInteractionGraph();
		virtual void initializeVariableInteractionGraph();
		void printVariableInteractionGraph();

		vec_t<vec_t<double>> getSimilarityMatrix( int similarity_measure_index );
		virtual double getSimilarityMeasure( size_t var_a, size_t var_b );
		
		bool isParameterInRangeBounds( double parameter, int dimension );
		virtual double getLowerRangeBound( int dimension );
		virtual double getUpperRangeBound( int dimension );
		
		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		void initializeObjectiveRotationMatrix( void );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		double *rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix );
		void ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size );

		void checkTermination();
		void checkEvaluationLimitTerminationCondition();
		void checkTimeLimitTerminationCondition();

	private:
		fitness_t( int number_of_variables, double vtr, bool use_vtr, opt_mode optimization_mode );
		virtual void evaluationFunction( solution_t<T> *solution ) = 0;
		virtual void partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution );

		vec_t<vec_t<double>> similarity_matrix;
};

template class fitness_t<char>;
template class fitness_t<double>;

}}
