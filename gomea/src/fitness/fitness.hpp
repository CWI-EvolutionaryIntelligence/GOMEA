#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"

namespace gomea{
namespace fitness{

class fitness_t
{
	public:
		fitness_t( int number_of_parameters, double vtr );

		// Properties
		std::string name;
		int number_of_parameters;
		double rotation_angle = 0.0;
		int rotation_block_size = 0;
		// Gray-box specific
		std::map<int,std::set<int>> variable_interaction_graph; // VIG[i] lists all variables dependent on variable x_i, excluding itself
		std::map<int,std::set<int>> subfunction_dependency_map; // map[i] lists all subfunctions dependent on variable x_i

		// Options
		double vtr; // value-to-reach
		short black_box_optimization = 0;
		short use_vtr = 0;
		short vtr_hit_status = 0;
		
		// Optimization progress
		double number_of_evaluations = 0.0; // discounted in GBO
		int full_number_of_evaluations = 0; // not discounted in GBO
		double elitist_objective_value = 1e308;
		double elitist_constraint_value = 1e308;

		virtual int getNumberOfSubfunctions();

		void evaluate( solution_t<double> *solution );
		void evaluatePartialSolution( solution_t<double> *parent, partial_solution_t<double> *solution );
		void evaluatePartialSolutionBlackBox( solution_t<double> *parent, partial_solution_t<double> *solution );
		
		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		void initializeObjectiveRotationMatrix( void );
		short isParameterInRangeBounds( double parameter, int dimension );
		virtual double getLowerRangeBound( int dimension );
		virtual double getUpperRangeBound( int dimension );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		double *rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix );

		double compute2DHyperVolume( const std::vector<double> &obj_f0, const std::vector<double> &obj_f1, std::vector<size_t> &sorted, double rx, double ry );
		double compute2DHyperVolume( double *obj_f0, double *obj_f1, int population_size );
		double compute2DUncrowdedHypervolume( double *obj_f0, double *obj_f1, int population_size );
		bool paretoDominates2D( double xf0, double xf1, double yf0, double yf1 );
		double distance_to_box(double ref_x, double ref_y, double p_x, double p_y);
		double distance_to_front(double p_x, double p_y, const std::vector<double>& obj_x, const std::vector<double>& obj_y, std::vector<size_t> &sorted_obj, double r_x, double r_y);

		void ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size );

		static fitness_t *getFitnessClass( int problem_index, int number_of_parameters, double vtr );
		static short betterFitness( solution_t<double> *sol_x, solution_t<double> *sol_y );
		static short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
		solution_t<double> *initializeSolution( int n );
		solution_t<double> *initializeSolution( double *variables );

		void initialize();	
		
		void initializeSubfunctionDependencyMap();
		virtual vec_t<int> inputsToSubfunction(int subfunction_index);

		bool hasVariableInteractionGraph();
		virtual void initializeVariableInteractionGraph();

	private:
		virtual void evaluationFunction( solution_t<double> *solution ) = 0;
		virtual void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );

		int evaluationEmbedded();
};

}}
