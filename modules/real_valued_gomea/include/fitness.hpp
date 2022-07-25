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

#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "tools.hpp"
#include "solution.hpp"
#include "partial_solution.hpp"
#include "fitness_buffer.hpp"
//#include "CECHeader.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

class fitness_t
{
	public:
		virtual ~fitness_t();

		// Properties
		std::string name;
		int number_of_parameters;
		double *lower_range_bound;
		double *upper_range_bound;
		double rotation_angle = 0.0;
		int rotation_block_size = 0;
		// Gray-box specific
		int number_of_subfunctions;		
		std::map<int,std::set<int>> variable_interaction_graph;

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

		static char *installedProblemName( int index );
		static int numberOfInstalledProblems( void );
		static void printAllInstalledProblems( void );

		void evaluate( solution_t *solution );
		void evaluatePartialSolution( solution_t *parent, partial_solution_t *solution );
		void evaluatePartialSolutionBlackBox( solution_t *parent, partial_solution_t *solution );
		
		void initializeFitnessFunction( void );
		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		void initializeObjectiveRotationMatrix( void );
		void initializeRangeBounds();
		short isParameterInRangeBounds( double parameter, int dimension );
		virtual double getLowerRangeBound( int dimension ) = 0;
		virtual double getUpperRangeBound( int dimension ) = 0;
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		double *rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix );

		double compute2DHyperVolume( const std::vector<double> &obj_f0, const std::vector<double> &obj_f1, std::vector<size_t> &sorted, double rx, double ry );
		double compute2DHyperVolume( double *obj_f0, double *obj_f1, int population_size );
		double compute2DUncrowdedHypervolume( double *obj_f0, double *obj_f1, int population_size );
		bool paretoDominates2D( double xf0, double xf1, double yf0, double yf1 );
		double distance_to_box(double ref_x, double ref_y, double p_x, double p_y);
		double distance_to_front(double p_x, double p_y, const std::vector<double>& obj_x, const std::vector<double>& obj_y, std::vector<size_t> &sorted_obj, double r_x, double r_y);

		void ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size );

		void sphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void sphereFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void ellipsoidFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void ellipsoidFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void cigarFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void tabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void cigarTabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void twoAxesFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void differentPowersFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void rosenbrockFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void parabolicRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void sharpRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void griewankFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void michalewiczFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void michalewiczFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void rastriginFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void rastriginFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void sumOfEllipsoidsFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void sumOfEllipsoidsFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
		void ciasBRFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
		void trapSphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );

		static fitness_t *getFitnessClass( int problem_index, int number_of_parameters, double vtr );
		static short betterFitness( solution_t *sol_x, solution_t *sol_y );
		static short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
		solution_t *initializeSolution( int n );
		solution_t *initializeSolution( double *variables );
		
		bool hasVariableInteractionGraph();
		virtual void initializeVariableInteractionGraph();

	private:
		virtual void evaluationFunction( solution_t *solution ) = 0;
		virtual void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );

		int evaluationEmbedded();
};

class sphereFunction_t : public fitness_t 
{
	public:
		sphereFunction_t( int number_of_parameters, double vtr );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double x );
};

class rosenbrockFunction_t : public fitness_t 
{
	public:
		rosenbrockFunction_t( int number_of_parameters, double vtr );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		void initializeVariableInteractionGraph();
		
	private:
		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		void univariatePartialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double x, double y );
};

class sorebFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		double conditioning_number, rotation_angle;
		int overlap_size;

		sorebFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, int block_size, int overlap_size );
		~sorebFunction_t();
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

		void initializeVariableInteractionGraph();

	private:
		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double *vars, int num_vars );

		int getStartingIndexOfBlock( int block_index );
		int getIndexOfFirstBlock( int var );
};

class osorebFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix_big, **rotation_matrix_small;
		int number_of_large_rotated_blocks;
		int number_of_small_rotated_blocks;

		osorebFunction_t( int number_of_parameters, double vtr );
		~osorebFunction_t();
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebChainFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		bool wrap_around;

		sorebChainFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around );
		~sorebChainFunction_t();
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebGridFunction_t : public fitness_t 
{
	public:
		bool wrap_around_x;
	   	bool wrap_around_y;
		int grid_width;
		std::map<int,double**> rotation_matrices;

		sorebGridFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y );
		~sorebGridFunction_t();
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		std::set<int> getNeighborsInGrid( int ind );

		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebCubeFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		bool wrap_around_x;
	   	bool wrap_around_y;
	   	bool wrap_around_z;
		int cube_width;
		std::map<int,double**> rotation_matrices;

		sorebCubeFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y, bool wrap_around_z );
		~sorebCubeFunction_t();
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		std::set<int> getNeighborsInGrid( int ind );

		void evaluationFunction( solution_t *solution );
		void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction( double *vars, int num_vars );
};

#ifdef CECLSGOFUNC
class CECLSGOFunctions_t: public fitness_t 
{
	public:
		CECLSGOFunctions_t( int id, int number_of_parameters, double vtr );
		Benchmarks *function;

		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		void evaluationFunction( solution_t *solution );
		//void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
};
#endif

class BD2FunctionHypervolume_t: public fitness_t 
{
	public:
		BD2FunctionHypervolume_t( int number_of_parameters, double vtr );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

		int front_size;
		int subfunction_size;
	private:
		void evaluationFunction( solution_t *solution );
		//void partialEvaluationFunction( solution_t *parent, partial_solution_t *solution );
		double subfunction_f0( double *x );
		double subfunction_f1( double *x );
};

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
double installedProblemLowerRangeBound( int index, int dimension );
double installedProblemUpperRangeBound( int index, int dimension );
void initializeParameterRangeBounds( void );
double sphereFunctionProblemLowerRangeBound( int dimension );
double sphereFunctionProblemUpperRangeBound( int dimension );
double ellipsoidFunctionLowerRangeBound( int dimension );
double ellipsoidFunctionUpperRangeBound( int dimension );
double cigarFunctionLowerRangeBound( int dimension );
double cigarFunctionUpperRangeBound( int dimension );
double tabletFunctionLowerRangeBound( int dimension );
double tabletFunctionUpperRangeBound( int dimension );
double cigarTabletFunctionLowerRangeBound( int dimension );
double cigarTabletFunctionUpperRangeBound( int dimension );
double twoAxesFunctionLowerRangeBound( int dimension );
double twoAxesFunctionUpperRangeBound( int dimension );
double differentPowersFunctionLowerRangeBound( int dimension );
double differentPowersFunctionUpperRangeBound( int dimension );
double rosenbrockFunctionLowerRangeBound( int dimension );
double rosenbrockFunctionUpperRangeBound( int dimension );
double parabolicRidgeFunctionLowerRangeBound( int dimension );
double parabolicRidgeFunctionUpperRangeBound( int dimension );
double sharpRidgeFunctionLowerRangeBound( int dimension );
double sharpRidgeFunctionUpperRangeBound( int dimension );
double griewankFunctionLowerRangeBound( int dimension );
double griewankFunctionUpperRangeBound( int dimension );
double michalewiczFunctionLowerRangeBound( int dimension );
double michalewiczFunctionUpperRangeBound( int dimension );
double rastriginFunctionLowerRangeBound( int dimension );
double rastriginFunctionUpperRangeBound( int dimension );
double sumOfEllipsoidsFunctionLowerRangeBound( int dimension );
double sumOfEllipsoidsFunctionUpperRangeBound( int dimension );
double ciasBRFunctionLowerRangeBound( int dimension );
double ciasBRFunctionUpperRangeBound( int dimension );
double trapSphereFunctionLowerRangeBound( int dimension );
double trapSphereFunctionUpperRangeBound( int dimension );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
