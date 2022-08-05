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
#include "fitness/fitness_basic.hpp"
#include "common/solution.hpp"
#include "common/partial_solution.hpp"
#include "common/gomea_defs.hpp"
#include "utils/tools.hpp"
#include "real_valued_gomea/cython/RealValuedGOMEA.h"
#include "fitness/cython/Fitness.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class sphereFunction_t : public fitness_t 
{
	public:
		sphereFunction_t( int number_of_parameters, double vtr );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		void univariatePartialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
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

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );
};

class BD2FunctionHypervolume_t: public fitness_t 
{
	public:
		BD2FunctionHypervolume_t( int number_of_parameters, double vtr );
		
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );

		int front_size;
		int subfunction_size;
	private:
		void evaluationFunction( solution_t<double> *solution );
		//void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction_f0( double *x );
		double subfunction_f1( double *x );
};


}}
