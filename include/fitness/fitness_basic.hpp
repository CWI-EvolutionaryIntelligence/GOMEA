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
#include "common/gomea_defs.hpp"
#include "common/solution.hpp"
#include "common/partial_solution.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

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

		void evaluate( solution_t<double> *solution );
		void evaluatePartialSolution( solution_t<double> *parent, partial_solution_t<double> *solution );
		void evaluatePartialSolutionBlackBox( solution_t<double> *parent, partial_solution_t<double> *solution );
		
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

		static fitness_t *getFitnessClass( int problem_index, int number_of_parameters, double vtr );
		static short betterFitness( solution_t<double> *sol_x, solution_t<double> *sol_y );
		static short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
		solution_t<double> *initializeSolution( int n );
		solution_t<double> *initializeSolution( double *variables );
		
		bool hasVariableInteractionGraph();
		virtual void initializeVariableInteractionGraph();

	private:
		virtual void evaluationFunction( solution_t<double> *solution ) = 0;
		virtual void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );

		int evaluationEmbedded();
};

}}
