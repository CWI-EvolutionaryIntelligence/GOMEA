#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/utils/time.hpp"
#include <map>

namespace gomea{
namespace fitness{

// Specifies whether a function is subject to minimization or maximization
typedef enum{
	MIN,
	MAX
} opt_mode;

class fitness_generic_t{
	public:
		bool isParameterInRangeBounds( double parameter, int dimension );
		virtual double getLowerRangeBound( int dimension );
		virtual double getUpperRangeBound( int dimension );
};

template<class T>
class fitness_t : public fitness_generic_t
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

        // Own copy of execution start time, for correct time condition checks
		time_t start_time;

		// Gray-box specific
		std::map<int,std::set<int>> variable_interaction_graph; // VIG[i] lists all variables dependent on variable x_i, excluding itself
		std::map<int,std::set<int>> subfunction_dependency_map; // map[i] lists all subfunctions dependent on variable x_i
		double rotation_angle = 0.0;
		int rotation_block_size = 0;
		double **rotation_matrix;

		// Options
		double vtr; // value-to-reach
		bool black_box_optimization = false;
		bool use_vtr;
		
		// Optimization progress
		double number_of_evaluations; // discounted in GBO
		int full_number_of_evaluations; // not discounted in GBO
		double elitist_objective_value;
		double elitist_constraint_value;
		bool vtr_hit_status;

		void evaluate( solution_t<T> *solution );
		void evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution );
		//void evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution, const std::set<int> &dependent_subfunctions );
		
		static fitness_t *getFitnessClass( int problem_index, int number_of_variables, double vtr );
		bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y );
		bool betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );

		virtual void initialize();
		void initializeRun();

		virtual void initializeVariableInteractionGraph();
		bool hasVariableInteractionGraph();
		void printVariableInteractionGraph();
		
		vec_t<vec_t<double>> getSimilarityMatrix( int similarity_measure_index );
		virtual double getSimilarityMeasure( size_t var_a, size_t var_b );
		
		
		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		double *rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix );
		void ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size );

		void checkTermination();
		void checkEvaluationLimitTerminationCondition();
		void checkTimeLimitTerminationCondition();

		int getNumberOfVariables();
		double getVTR();

	private:
		fitness_t( int number_of_variables, double vtr, bool use_vtr, opt_mode optimization_mode );
		virtual void evaluationFunction( solution_t<T> *solution ) = 0;
		//virtual void partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution, const std::set<int> &dependent_subfunctions );
		virtual void partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution );
		
		void evaluatePartialSolutionBlackBox( solution_t<T> *parent, partial_solution_t<T> *solution );

		vec_t<vec_t<double>> similarity_matrix;
};

}}
